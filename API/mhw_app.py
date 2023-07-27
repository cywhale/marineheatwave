import xarray as xr
import pandas as pd
import polars as pl
import numpy as np
from fastapi import FastAPI, HTTPException, Query, Response
from fastapi.responses import JSONResponse, FileResponse
from fastapi.openapi.docs import get_swagger_ui_html
from fastapi.openapi.utils import get_openapi
from typing import Optional, List
from pydantic import BaseModel
from datetime import date, datetime, timedelta
from tempfile import NamedTemporaryFile
import json, dask
from multiprocessing.pool import Pool
dask.config.set(pool=Pool(4))


def generate_custom_openapi():
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="Open API of Marine Heatwaves",
        version="1.0.0",
        description="Marine heatwaves (MHWs) evaluated from 0.25-degree gridded NOAA OISST v2.1, compiled by ODB, Taiwan. \n" +
                    "Reference: Jacox, Alexander, Bograd, and Scott (2020), Thermal Displacement by Marine Heatwaves, Nature, 584, 82–86, doi:10.1038/s41586-020-2534-z",
        routes=app.routes,
    )
    openapi_schema["servers"] = [
        {
            "url": "https://eco.odb.ntu.edu.tw"
        }
    ]
    # cause error: at responses.200/422.content.application/json.schema.$ref
    # openapi_schema["components"].pop("schemas", None)
    # if "schemas" in openapi_schema["components"]:
    #     openapi_schema["components"]["schemas"].pop("HTTPValidationError", None)
    #     openapi_schema["components"]["schemas"].pop("ValidationError", None)
    app.openapi_schema = openapi_schema
    return app.openapi_schema


app = FastAPI(docs_url=None)


@app.get("/api/swagger/openapi.json", include_in_schema=False)
async def custom_openapi():
    # app.openapi()) modify to customize openapi.json
    return JSONResponse(generate_custom_openapi())


@app.get("/api/swagger", include_in_schema=False)
async def custom_swagger_ui_html():
    return get_swagger_ui_html(
        openapi_url="/api/swagger/openapi.json",  # app.openapi_url
        title=app.title
    )

### Global variables ###
gridSz = 0.25
timeLimit = 365
LON_RANGE_LIMIT = 90
LAT_RANGE_LIMIT = 90
AREA_LIMIT = LON_RANGE_LIMIT * LAT_RANGE_LIMIT

@app.on_event("startup")
async def startup():
    global dz
    dz = xr.open_zarr('sst_anomaly.zarr', chunks='auto',
                      group='anomaly', decode_times=True)


def deg2str(value, isLon=True, oriNeg=False, roundTo=3):
    if (isLon):
      # if original lon is -180 -> 180 + 360 = 180 -> deg2str -> 180 cannot convert back to -180
      lon = value - 360 if value > 180 or oriNeg else value
      return str(round(lon, roundTo)).replace(".", "d")

    return str(round(value, roundTo)).replace(".", "d")

def to_nearest_grid_point(lon: float, lat: float) -> tuple:
    mlon = 180 if lon > 180 else (-180 if lon < -180 else lon)
    mlat = 90 if lat > 90 else (-90 if lat < -90 else lat)
    grid_lon = round(mlon * 4) / 4
    grid_lat = round(mlat * 4) / 4
    grid_lon = grid_lon + 360 if grid_lon < 0 else grid_lon
    return (grid_lon, grid_lat)


async def process_mhw_data(lon0: float, lat0: float, lon1: Optional[float], lat1: Optional[float], start: Optional[date], end: Optional[date], append: Optional[str]):
    if append is None:
        append = 'level'

    variables = [var.strip() for var in append.split(
        ',') if var.strip() in ['level', 'sst_anomaly', 'td']]
    if not variables:
        raise HTTPException(
            status_code=400, detail="Invalid variable(s). Allowed variables are 'sst_anomaly', 'level', and 'td'.")

    variables.sort() #in-place sort not return anything
    out_file = '-'.join(variables)

    if start is None:
        start_date = pd.to_datetime('1982-01-01')
    else:
        try:
            start_date = pd.to_datetime(start)
            start_date = start_date.replace(day=1)  # set day to 1
        except ValueError:
            raise HTTPException(
                status_code=400, detail="Invalid start date format")

    if end is None:
        end_date = pd.to_datetime(datetime.today())
    else:
        try:
            end_date = pd.to_datetime(end)
            end_date = end_date.replace(day=1)  # set day to 1
        except ValueError:
            raise HTTPException(
                status_code=400, detail="Invalid end date format")

    try:
        orig_lon0, orig_lon1 = lon0, lon1
        lon0, lat0 = to_nearest_grid_point(lon0, lat0)

        if lon1 is None or lat1 is None:
            # Only one point, no date range limitation
            out_file = (out_file + '_' + deg2str(lon0, True, orig_lon0<=-179.875) + '_' + deg2str(lat0, False) +
                        '_' + 'Point' + '_' + str(start_date.date()) + '_' + str(end_date.date()))
            data_subset = dz.sel(lon=slice(lon0, lon0+gridSz-0.01), lat=slice(
                lat0, lat0+gridSz-0.01), date=slice(start_date, end_date))

        else:
            # Bounding box, 1 month or 1 year date range limitation
            lon1, lat1 = to_nearest_grid_point(lon1, lat1)

            # Adjust temporal range limit based on spatial range
            lon_range = abs(orig_lon1 - orig_lon0)
            lat_range = abs(lat1 - lat0)
            area_range = lon_range * lat_range

            if (lon_range > LON_RANGE_LIMIT and lat_range > LAT_RANGE_LIMIT) or (area_range > AREA_LIMIT):
                # Large spatial range, limit to one month of data
                end_date = start_date + pd.DateOffset(months=1) - timedelta(days=1)
            elif (end_date - start_date).days > timeLimit:
                # Smaller spatial range, limit to one year of data
                end_date = start_date + timedelta(days=timeLimit)

            if lon0 > lon1 and np.sign(orig_lon0) == np.sign(orig_lon1):
                # Swap if lon0 > lon1 but the same sign
                lon0, lon1 = lon1, lon0
                orig_lon0, orig_lon1 = orig_lon1, orig_lon0

            if np.sign(orig_lon0) != np.sign(orig_lon1):
                # Requested area crosses the zero meridian
                if orig_lon1 < 0:
                    # Swap if orig_lon1 < 0 and now 180 < lon1 < 360
                    lon0, lon1 = lon1, lon0
                    orig_lon0, orig_lon1 = orig_lon1, orig_lon0

                subset1 = dz.sel(lon=slice(lon0, 360), lat=slice(lat0, lat1), date=slice(start_date, end_date))
                subset2 = dz.sel(lon=slice(0, lon1), lat=slice(lat0, lat1), date=slice(start_date, end_date))
                data_subset = xr.concat([subset1, subset2], dim='lon')
            else:
                # Requested area doesn't cross the zero meridian
                data_subset = dz.sel(lon=slice(lon0, lon1), lat=slice(lat0, lat1), date=slice(start_date, end_date))

            out_file = (out_file + '_' + deg2str(lon0, True, orig_lon0<=-179.875) + '_' + deg2str(lat0, False) +
                        '_' + deg2str(lon1, True, orig_lon1<=-179.875) + '_' + deg2str(lat1, False) +
                        '_' + 'BBOX' + '_' + str(start_date.date()) + '_' + str(end_date.date()))

        if data_subset.nbytes == 0:
            raise HTTPException(
                status_code=400, detail="No data available for the given parameters.")

        df = data_subset.to_dataframe().reset_index()
        mask = df['lon'] > 180
        df.loc[mask, 'lon'] = df.loc[mask, 'lon'] - 360
        df = df[['lon', 'lat', 'date'] +
                variables].dropna(how='all', subset=variables)

        if df.empty:
            raise HTTPException(
                status_code=400, detail="No data available after removing rows with NA values.")

        # if 'date' in df.columns:
        df['date'] = df['date'].apply(
            lambda x: x.isoformat() if not pd.isnull(x) else '')
        ##df.replace([np.inf, -np.inf], np.nan, inplace=True)  # Replace infinite values with NaN
        ##df.fillna("", inplace=True)  # Replace NaN values with None #not "null"
        # return df #Pandas version
        return out_file, pl.from_pandas(df)

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


class MHWResponse(BaseModel):
    lon: float
    lat: float
    date: date
    level: Optional[float]
    sst_anomaly: Optional[float]
    td: Optional[float]


@app.get("/api/mhw", response_model=List[MHWResponse], tags=["Marine Heatwave"], summary="Query MHW data as JSON")
async def read_mhw(
    lon0: float = Query(...,
                        description="Minimum longitude, range: [-180, 180]"),
    lat0: float = Query(..., description="Minimum latitude, range: [-90, 90]"),
    lon1: Optional[float] = Query(
        None, description="Maximum longitude, range: [-180, 180]"),
    lat1: Optional[float] = Query(
        None, description="Maximum latitude, range: [-90, 90]"),
    start: Optional[date] = Query(
        None, description="Start date of MHWs to query, minimum is 1982-02-01"),
    end: Optional[date] = Query(
        None, description="End date of MHWs to query, maximum is one month before the current date"),
    append: Optional[str] = Query(
        None, description="Data fields to append, separated by commas. Allowed fields: 'sst_anomaly', 'level', 'td'")
):
    """
    Query MHW data by longitude/latitude/date (in JSON).

    #### Usage
    * One-point MHWs without time-span limitation: e.g. /api/mhw?lon0=135&lat0=15
    * Bounding-box <= 90x90 in degrees: 1-year time-span limitation: e.g. /api/mhw?lon0=135&lon1&=140&lat0=15&lat1=30&start=2021-01-01
    * Bounding-box > 90x90 in degrees: 1-month time-span limitation: e.g. /api/mhw?lon0=-180&lon1&=180&lat0=-90&lat1=90&start=2021-01-01
    """

    _, df = await process_mhw_data(lon0, lat0, lon1, lat1, start, end, append)
    ##return JSONResponse(content=df.to_dict(orient='records')) #cannot handle NA when JSONResponse
    # return JSONResponse(content=json.loads(df.to_json(orient='records'))) # work version in Pandas
    return JSONResponse(content=df.to_dicts())

@app.get("/api/mhw/csv", response_class=Response, tags=["Marine Heatwave"], summary="Query MHW data as CSV")
async def read_mhw_csv(
    lon0: float = Query(...,
                        description="Minimum longitude, range: [-180, 180]"),
    lat0: float = Query(..., description="Minimum latitude, range: [-90, 90]"),
    lon1: Optional[float] = Query(
        None, description="Maximum longitude, range: [-180, 180]"),
    lat1: Optional[float] = Query(
        None, description="Maximum latitude, range: [-90, 90]"),
    start: Optional[date] = Query(
        None, description="Start date of MHWs to query, minimum is 1982-02-01"),
    end: Optional[date] = Query(
        None, description="End date of MHWs to query, maximum is one month before the current date"),
    append: Optional[str] = Query(
        None, description="Data fields to append, separated by commas. Allowed fields: 'sst_anomaly', 'level', 'td'")
):
    """
    Query MHW data by longitude/latitude/date (in csv).

    #### Usage
    * One-point MHWs without time-span limitation: e.g. /api/mhw?lon0=135&lat0=15
    * Bounding-box <= 90x90 in degrees: 1-year time-span limitation: e.g. /api/mhw?lon0=135&lon1&=140&lat0=15&lat1=30&start=2021-01-01
    * Bounding-box > 90x90 in degrees: 1-month time-span limitation: e.g. /api/mhw?lon0=-180&lon1&=180&lat0=-90&lat1=90&start=2021-01-01
    """

    out_file, df = await process_mhw_data(lon0, lat0, lon1, lat1, start, end, append)
    temp_file = NamedTemporaryFile(delete=False)
    # df.to_csv(temp_file.name, index=False) #Pandas solution, work
    df.write_csv(temp_file.name) #polars version

    return FileResponse(temp_file.name, media_type="text/csv", filename=out_file+".csv")
