import xarray as xr
import pandas as pd
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import JSONResponse, FileResponse
from fastapi.openapi.docs import get_swagger_ui_html
from fastapi.openapi.utils import get_openapi
from typing import Optional, List
from pydantic import BaseModel
from datetime import date, datetime, timedelta
from tempfile import NamedTemporaryFile
import dask
from multiprocessing.pool import Pool
dask.config.set(pool=Pool(4))


def generate_custom_openapi():
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="Open API of Marine Heatwaves, ODB",
        version="1.0.0",
        description="Marine heatwaves (MHWs) evaluated from 0.25-degree gridded NOAA OISST v2.1.\n" +
                    "Reference: Jacox, Alexander, Bograd, and Scott (2020), Thermal Displacement by Marine Heatwaves, Nature, 584, 82–86, doi:10.1038/s41586-020-2534-z",
        routes=app.routes,
    )
    openapi_schema["servers"] = [
        {
            "url": "https://eco.odb.ntu.edu.tw"
        }
    ]
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


@app.on_event("startup")
async def startup():
    global dz
    dz = xr.open_zarr('sst_anomaly.zarr', chunks='auto',
                      group='anomaly', decode_times=True)


def to_nearest_grid_point(lon: float, lat: float) -> tuple:
    grid_lon = round(lon * 4) / 4
    grid_lat = round(lat * 4) / 4
    return (grid_lon, grid_lat)


async def process_mhw_data(lon0: float, lat0: float, lon1: Optional[float], lat1: Optional[float], start: Optional[date], end: Optional[date], append: Optional[str]):
    if append is None:
        append = 'level'

    variables = [var.strip() for var in append.split(
        ',') if var.strip() in ['level', 'sst_anomaly', 'td']]
    if not variables:
        raise HTTPException(
            status_code=400, detail="Invalid variable(s). Allowed variables are 'sst_anomaly', 'level', and 'td'.")

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
        lon0, lat0 = to_nearest_grid_point(lon0, lat0)

        if lon1 is None or lat1 is None:
            # Only one point, no date range limitation
            data_subset = dz.sel(lon=slice(lon0+180, lon0+180+gridSz-0.01), lat=slice(
                lat0, lat0+gridSz-0.01), date=slice(start_date, end_date))

        else:
            # Bounding box, 1 year date range limitation
            lon1, lat1 = to_nearest_grid_point(lon1, lat1)
            if (end_date - start_date).days > timeLimit:
                end_date = start_date + timedelta(days=timeLimit)

            data_subset = dz.sel(lon=slice(
                lon0+180, lon1+180), lat=slice(lat0, lat1), date=slice(start_date, end_date))

        if data_subset.nbytes == 0:
            raise HTTPException(
                status_code=400, detail="No data available for the given parameters.")

        df = data_subset.to_dataframe().reset_index()
        df['lon'] = df['lon'] - 180
        df = df[['lon', 'lat', 'date'] +
                variables].dropna(how='all', subset=variables)

        if df.empty:
            raise HTTPException(
                status_code=400, detail="No data available after removing rows with NA values.")

        # if 'date' in df.columns:
        df['date'] = df['date'].apply(
            lambda x: x.isoformat() if not pd.isnull(x) else '')
        return df

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))


class MHWResponse(BaseModel):
    lon: float
    lat: float
    date: date
    level: Optional[float]
    sst_anomaly: Optional[float]
    td: Optional[float]


@app.get("/api/mhw", response_model=List[MHWResponse])
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
    df = await process_mhw_data(lon0, lat0, lon1, lat1, start, end, append)
    return JSONResponse(content=df.to_dict(orient='records'))


@app.get("/api/mhw/csv")
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
    df = await process_mhw_data(lon0, lat0, lon1, lat1, start, end, append)
    temp_file = NamedTemporaryFile(delete=False)
    df.to_csv(temp_file.name, index=False)
    startx = '1982-01-01' if start is None else str(start)
    return FileResponse(temp_file.name, media_type="text/csv", filename="mhw_" + startx + ".csv")
