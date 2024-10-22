import xarray as xr
import pandas as pd
# import polars as pl
# import numpy as np
from fastapi import FastAPI, HTTPException, Query, Response
from fastapi.responses import JSONResponse, FileResponse, StreamingResponse
from fastapi.openapi.docs import get_swagger_ui_html
from fastapi.openapi.utils import get_openapi
from contextlib import asynccontextmanager
from typing import Optional, List
from pydantic import BaseModel
from datetime import date, datetime  # , timedelta
from tempfile import NamedTemporaryFile
from src.mhw_utils import process_mhw_data
from src.mhw_plot import month_climatology, region_climatology, period2date
import src.config as config
# import json, dask
# This is a more low-level way of setting up Dask for parallel computations, and it doesn't provide some of the additional features of Dask's own distributed scheduler, such as advanced task prioritization and data locality awareness.
# from multiprocessing.pool import Pool
# dask.config.set(pool=Pool(4))
# from dask.distributed import Client
# client = Client('tcp://localhost:8786')
from src.dask_client_manager import get_dask_client
client = get_dask_client("mhwapi")

#app = FastAPI(docs_url=None)
#@app.on_event("startup")
#async def startup():
@asynccontextmanager
async def lifespan(app: FastAPI):
    # global dz # old use 'sst_anomaly.zarr', add sst -> mhw.zarr
    config.dz = xr.open_zarr('data/mhw.zarr', chunks='auto',
                             group='anomaly', decode_times=True)
    config.gridSz = 0.25
    config.timeLimit = 365
    config.LON_RANGE_LIMIT = 90
    config.LAT_RANGE_LIMIT = 90
    config.AREA_LIMIT = config.LON_RANGE_LIMIT * config.LAT_RANGE_LIMIT
    yield
    print("Application is shutting down!")
    config.dz.close()
    client.close()

#@app.on_event("shutdown")
#async def shutdown_event():
#    # Code to run when the application shuts down
#    print("Application is shutting down!")
#    client.close()


app = FastAPI(docs_url=None, lifespan=lifespan)


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


@app.get("/api/swagger/mhw/openapi.json", include_in_schema=False)
async def custom_openapi():
    # app.openapi()) modify to customize openapi.json
    return JSONResponse(generate_custom_openapi())


@app.get("/api/swagger/mhw", include_in_schema=False)
async def custom_swagger_ui_html():
    return get_swagger_ui_html(
        openapi_url="/api/swagger/mhw/openapi.json",  # app.openapi_url
        title=app.title
    )

### Global variables ###
### move to config.py ##


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
        None, description="Start date of MHWs to query, minimum is 1982-01-01"),
    end: Optional[date] = Query(
        None, description="End date of MHWs to query, maximum is one month before the current date"),
    append: Optional[str] = Query(
        None, description="Data fields to append, separated by commas. Allowed fields: 'sst', 'sst_anomaly', 'level', 'td'")
):
    """
    Query MHW data by longitude/latitude/date (in JSON).

    #### Usage
    * One-point MHWs without time-span limitation: e.g. /api/mhw?lon0=135&lat0=15
    * Bounding-box <= 90x90 in degrees: 1-year time-span limitation: e.g. /api/mhw?lon0=135&lon1&=140&lat0=15&lat1=30&start=2021-01-01
    * Bounding-box > 90x90 in degrees: 1-month time-span limitation: e.g. /api/mhw?lon0=-180&lon1&=180&lat0=-90&lat1=90&start=2021-01-01
    """

    try:
        _, df = await process_mhw_data(lon0=lon0, lat0=lat0, lon1=lon1, lat1=lat1, start=start, end=end, append=append)
        if df.is_empty():
            raise HTTPException(status_code=400, detail="No data available for the given parameters.")
        # return JSONResponse(content=df.to_dict(orient='records')) #cannot handle NA when JSONResponse
        # return JSONResponse(content=json.loads(df.to_json(orient='records'))) # work version in Pandas
        return JSONResponse(content=df.to_dicts())

    except HTTPException as herr:
        raise herr
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail="Internal server error. Please try it later or inform admin")


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
        None, description="Start date of MHWs to query, minimum is 1982-01-01"),
    end: Optional[date] = Query(
        None, description="End date of MHWs to query, maximum is one month before the current date"),
    append: Optional[str] = Query(
        None, description="Data fields to append, separated by commas. Allowed fields: 'sst', 'sst_anomaly', 'level', 'td'")
):
    """
    Query MHW data by longitude/latitude/date (in csv).

    #### Usage
    * One-point MHWs without time-span limitation: e.g. /api/mhw?lon0=135&lat0=15
    * Bounding-box <= 90x90 in degrees: 1-year time-span limitation: e.g. /api/mhw?lon0=135&lon1&=140&lat0=15&lat1=30&start=2021-01-01
    * Bounding-box > 90x90 in degrees: 1-month time-span limitation: e.g. /api/mhw?lon0=-180&lon1&=180&lat0=-90&lat1=90&start=2021-01-01
    """

    try:
        out_file, df = await process_mhw_data(lon0=lon0, lat0=lat0, lon1=lon1, lat1=lat1, start=start, end=end, append=append)
        if df.is_empty():
            raise HTTPException(status_code=400, detail="No data available for the given parameters.")

        temp_file = NamedTemporaryFile(delete=False)
        # df.to_csv(temp_file.name, index=False) #Pandas solution, work
        df.write_csv(temp_file.name)  # polars version

        return FileResponse(temp_file.name, media_type="text/csv", filename=out_file+".csv")

    except HTTPException as herr:
        raise herr
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail="Internal server error. Please try it later or inform admin")

@app.get("/api/mhw/plot", tags=["MHW_plot"], summary="Plot SST or SST Anomalies climatology as PNG", include_in_schema=False)
async def climatology(
    bbox: str = Query(..., description="Bounding-box (x0,y0,x1,y1),(...). Multi-bbox is allowed and x1, y1 can be skipped as point"),
    period: str = Query(..., description="Time period %Y%m%d-%Y%m%d,... Multi-period is allowed, and %m:month %d:day can be skipped"),
    sstvar: Optional[str] = Query(
        None, description="Plot climatology for sst(default): SST or sst_anomaly: SST Anomalies"),
    mode: Optional[str] = Query(
        None, description="Plot series(default): time-series or month: average by month of {sstvar} throughout {period}"),
    palette: Optional[str] = Query(
        None, description="Customize color plate for each {bbox} if {mode} is 'series' or for each {period} if {mode} is 'month'.")
):
    """
    Plot SST or SST Anomalies climatology (in PNG)

    #### Usage
    * Note that wider time/spatial span specified by {period} and {bbox} lead to much slower response.
    * {mode} = series: Time-series e.g. /api/mhw/plot?bbox=(-90,-10,-80,0),(-150,-5,-90,5)&period=2015-202306&sstvar=sst_anomaly
    * {mode} = month: Month climatology. /api/mhw/plot?bbox=(-90,-60,0,60)&sstvar=sst&period=1982-2011,2012-2021,2022,2023&mode=month
    """
    if not mode or mode.strip() not in ['month', 'series']:
        mode = 'series'

    if mode == 'month':
        bbox_list = bbox.strip('()').split('),(')
        bbox_coords = [float(coord) for coord in bbox_list[0].split(',')]
        if len(bbox_coords) == 2:
            lon0, lat0, lon1, lat1 = bbox_coords[0], bbox_coords[1], None, None
        elif len(bbox_coords) == 4:
            lon0, lat0, lon1, lat1 = bbox_coords
        else:
            raise HTTPException(
                status_code=400, detail="Invalid bounding box format. Must contain 2 or 4 coordinates.")

        buf = await month_climatology(lon0=lon0, lat0=lat0, lon1=lon1, lat1=lat1, period=period, sstvar=sstvar, palette=palette)

    else:
        if period == "all":
            start = '1982-01-01'
            end = str((datetime.today() - pd.DateOffset(months=1)
                       ).replace(day=1).date())
        else:
            start, end = period2date(period.strip().split(',')[0])

        buf = await region_climatology(bbox, start, end, sstvar, palette)

    return StreamingResponse(buf, media_type="image/png")
