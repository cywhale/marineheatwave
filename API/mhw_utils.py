import xarray as xr
import pandas as pd
import polars as pl
import numpy as np
from fastapi import HTTPException
from typing import Optional
from datetime import date, datetime, timedelta
import config


def deg2str(value, isLon=True, oriNeg=False, roundTo=3):
    if (isLon):
        # if original lon is -180 -> 180 + 360 = 180 -> deg2str -> 180 cannot convert back to -180
        lon = value - 360 if value > 180 or oriNeg else value
        return str(round(lon, roundTo)).replace(".", "d")

    return str(round(value, roundTo)).replace(".", "d")


def to_nearest_grid_point(lon: float, lat: float) -> tuple:
    mlon = 180 if lon > 180 else (-180 if lon < -180 else lon)
    mlat = 90 if lat > 90 else (-90 if lat < -90 else lat)
    mlon = mlon + 360 if mlon < 0 else mlon
    grid_lon = round(mlon * 4) / 4
    grid_lat = round(mlat * 4) / 4
    return (grid_lon, grid_lat)


def output_df(df, isoTime=False):
    if df.empty:
        raise HTTPException(
            status_code=400, detail="No data available after removing rows with NA values.")

    pl_df = pl.from_pandas(df)
    if (isoTime):
        pl_df = pl_df.with_columns(pl.col('date').dt.strftime('%Y-%m-%d'))
    # if 'date' in df.columns: #### its slow!
    # df['date'] = df['date'].apply(
    # lambda x: x.isoformat() if not pd.isnull(x) else '')
    # df.replace([np.inf, -np.inf], np.nan, inplace=True)  # Replace infinite values with NaN
    # df.fillna("", inplace=True)  # Replace NaN values with None #not "null"
    # return df #Pandas version
    return pl_df


async def process_mhw_data(lon0: float, lat0: float, lon1: Optional[float], lat1: Optional[float], start: Optional[date], end: Optional[date], append: Optional[str], mode: Optional[str] = None):
    if mode is None:
        mode = 'raw'

    if append is None:
        append = 'level'

    variables = list(set([var.strip() for var in append.split(
        ',') if var.strip() in ['level', 'sst_anomaly', 'td', 'sst']]))
    if not variables:
        raise HTTPException(
            status_code=400, detail="Invalid variable(s). Allowed variables are 'sst_anomaly', 'level', and 'td'.")

    variables.sort()  # in-place sort not return anything
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

        if lon1 is None or lat1 is None or (orig_lon0 == orig_lon1 and lat0 == lat1):
            # Only one point, no date range limitation
            out_file = (out_file + '_' + deg2str(lon0, True, orig_lon0 <= -179.875) + '_' + deg2str(lat0, False) +
                        '_' + 'Point' + '_' + str(start_date.date()) + '_' + str(end_date.date()))
            data_subset = config.dz.sel(lon=slice(lon0-0.5*config.gridSz, lon0+0.5*config.gridSz-0.001), lat=slice(
                lat0-0.5*config.gridSz, lat0+0.5*config.gridSz-0.001), date=slice(start_date, end_date))

        else:
            # Bounding box, 1 month or 1 year date range limitation
            lon1, lat1 = to_nearest_grid_point(lon1, lat1)

            if (mode != 'area_mean_sst' and mode != 'area_mean_sst_anomaly' and mode != 'month_mean'):
                # Adjust temporal range limit based on spatial range
                lon_range = abs(orig_lon1 - orig_lon0)
                lat_range = abs(lat1 - lat0)
                area_range = lon_range * lat_range

                if (lon_range > config.LON_RANGE_LIMIT and lat_range > config.LAT_RANGE_LIMIT) or (area_range > config.AREA_LIMIT):
                    # Large spatial range, limit to one month of data
                    end_date = start_date + \
                        pd.DateOffset(months=1) - timedelta(days=1)
                elif (end_date - start_date).days > config.timeLimit:
                    # Smaller spatial range, limit to one year of data
                    end_date = start_date + timedelta(days=config.timeLimit)

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

                subset1 = config.dz.sel(lon=slice(lon0, 360), lat=slice(
                    lat0, lat1), date=slice(start_date, end_date))
                subset2 = config.dz.sel(lon=slice(0, lon1), lat=slice(
                    lat0, lat1), date=slice(start_date, end_date))
                data_subset = xr.concat([subset1, subset2], dim='lon')
            else:
                # Requested area doesn't cross the zero meridian
                data_subset = config.dz.sel(lon=slice(lon0, lon1), lat=slice(
                    lat0, lat1), date=slice(start_date, end_date))

            out_file = (out_file + '_' + deg2str(lon0, True, orig_lon0 <= -179.875) + '_' + deg2str(lat0, False) +
                        '_' + deg2str(lon1, True, orig_lon1 <= -179.875) + '_' + deg2str(lat1, False) +
                        '_' + 'BBOX' + '_' + str(start_date.date()) + '_' + str(end_date.date()))

        if data_subset.nbytes == 0:
            raise HTTPException(
                status_code=400, detail="No data available for the given parameters.")

        data_subset.load()

        if mode in ['area_mean_sst', 'area_mean_sst_anomaly', 'month_mean']:
            # Perform the averaging using xarray
            area_mean = data_subset.mean(dim=['lon', 'lat'], skipna=True).dropna(
                dim='date', how='all').compute(timeout='30s')
            # Run the compute-intensive tasks in the Dask cluster
            # future = client.submit(process_data_subset, data_subset) #, append=append
            # area_mean = await client.gather(future)
            df = area_mean[['date', append]].to_dataframe().reset_index()
            # df.dropna(subset=[append], inplace=True)
            pl_df = output_df(df, False)

        else:
            df = data_subset.to_dataframe().reset_index()
            mask = df['lon'] > 180
            df.loc[mask, 'lon'] = df.loc[mask, 'lon'] - 360
            df = df[['lon', 'lat', 'date'] +
                    variables].dropna(how='all', subset=variables)
            pl_df = output_df(df, True)

        return out_file, pl_df

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
