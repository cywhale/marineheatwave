import pandas as pd
import xarray as xr
import numpy as np
import calendar, requests, os  # noqa: E401
from sqlalchemy import create_engine
from dotenv import load_dotenv
from urllib.parse import quote 

# Define the paths and connection strings
ZARR_PATH = 'data/mhw.zarr'
DATA_PATH = 'data/tmp/'
RECHUNK = False  #note, if re-chunk needed, it will overwrite whole datasets
TEST_Download = False
TEST_date= '2023-09-01'

load_dotenv()
DBUSER = os.getenv('DBUSER')
DBPASS = os.getenv('DBPASS')
DBHOST = os.getenv('DBHOST')
DBPORT = os.getenv('DBPORT')
DBNAME = os.getenv('DBNAME')
MHWTABLE = os.getenv('MHWTABLE')

PG_CONNECTION_STRING = 'postgresql://' + DBUSER + ':%s@' + DBHOST + ':' + DBPORT + '/' + DBNAME

def download_noaa_data(date, dest_dir):
    """
    Downloads NOAA OISST data for the given month and year.
    Args:
    - date (str): YYYY-MM-DD format
    - dest_dir (str): Destination directory to save the file
    Returns:
    - filename (str) if successful, None otherwise
    """
    year_month = date[:7].replace('-', '')
    year = int(year_month[:4])
    month = int(year_month[4:6])
    _, num_days = calendar.monthrange(year, month)

    downloaded_files = []
    
    for day in range(1, num_days + 1):
        filename = f"oisst-avhrr-v02r01.{year_month}{day:02d}.nc"
        file_path = os.path.join(dest_dir, filename)

        if os.path.exists(file_path):
            print(f"File {filename} already exists in the destination directory.")
            downloaded_files.append(file_path)
            continue

        url = f"https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/{year_month}/{filename}"

        response = requests.get(url, stream=True)
        
        if response.status_code == 200:
            with open(file_path, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192):
                    file.write(chunk)
            print(f"Downloaded {filename} to {dest_dir}")
            downloaded_files.append(file_path)
        else:
            print(f"Failed to download data from {url}. Error: {response.status_code}")

    return downloaded_files

def append_to_zarr():
    dz = xr.open_zarr(ZARR_PATH, group='anomaly', decode_times=True)

    # Check duplicates in the existing Zarr dataset
    #duplicate_dates_dz = dz['date'].to_index().duplicated()
    #if duplicate_dates_dz.any():
    #    print("Error! Duplicate dates in dz:", dz['date'][duplicate_dates_dz].values)
    #    return

    # Check the last date available in the Zarr dataset
    last_date_in_zarr = pd.to_datetime(dz['date'].values[-1])
    print("Get max date from zarr: ", last_date_in_zarr)

    # --- Append sst_anomaly, level, td data from Postgres ---
    engine = create_engine(PG_CONNECTION_STRING % quote(DBPASS))
    query = f"SELECT lat, lon, sst_anomaly, level, td, date FROM sst_anomaly_without_detrend WHERE date > '{last_date_in_zarr}';"
    df = pd.read_sql_query(query, engine)
    engine.dispose()

    if not df.empty:
        print("---- Merge data from DB start ----")
        df['date'] = pd.to_datetime(df['date'])
        dfs = [df[['date', 'lat', 'lon', var]].copy() for var in ['sst_anomaly', 'level', 'td']]

        # Merge dataframes on date, lat, lon
        df_final = dfs[0]
        for df_var in dfs[1:]:
            df_final = df_final.merge(df_var, on=['date', 'lat', 'lon'], how='outer')

        # Convert to xarray Dataset
        ds_db = df_final.set_index(['date', 'lat', 'lon']).to_xarray()

        # Merge this data with the existing Zarr dataset
        dz = xr.concat([dz, ds_db], dim='date')

        # --- Append sst data from the downloaded netcdf ---
        print("---- Merge data from Netcdf start ----")
        # Calculate next_month_date based on the last_date_in_zarr
        if last_date_in_zarr.month == 12:
            next_month_date = f"{last_date_in_zarr.year + 1}-01-01"
        else:
            next_month_date = f"{last_date_in_zarr.year}-{last_date_in_zarr.month + 1:02}-01"

        # Download the NOAA data for the next month
        filenames = download_noaa_data(next_month_date, DATA_PATH)

        if filenames is None or len(filenames) == 0:
            print("Error downloading NOAA data for date: ", next_month_date, ". Aborting update.")
            return
        
        print("Download nc for the date: ", next_month_date, " and get file: ", filenames)

        # ds_nc = xr.open_mfdataset(filenames, parallel=True, chunks={'time': '500MB'})
        try:
            with xr.open_mfdataset(filenames, parallel=False, chunks={'time': '500MB'}) as ds_nc:
                msst = ds_nc["sst"].resample(time='1MS').mean()
                ds_msst = msst.compute()
                ds_msst = ds_msst.squeeze('zlev').rename({'time': 'date'}).drop('zlev')
                print("Processing mean SST:", ds_nc)

                # Check if the 'next_month_date' already exists in dz['date']
                if np.datetime64(next_month_date) in dz['date'].values:
                    # Align the ds_msst Dataset with dz along 'lat' and 'lon'
                    # print(ds_msst)
                    # print("----debugging after align----")
                    ds_msst_aligned, _ = xr.align(ds_msst, dz['sst'], join='inner', exclude=['date'])
                    # print(ds_msst_aligned)

                    # Update the 'sst' values in dz for the specific date
                    dz['sst'].loc[dict(date=next_month_date)] = ds_msst_aligned.sel(date=next_month_date)
                else:
                    # If 'next_month_date' doesn't exist in dz['date'], simply concatenate as before
                    dz['sst'] = xr.concat([dz['sst'], ds_msst], dim='date')

        except Exception as e:
            print("Error processing NetCDF files: ", e)
            return
        
        if RECHUNK:
            # Re-chunk and save to Zarr
            chunk_size_date = dz['level'].chunks[0][0]
            chunk_size_lat = dz['level'].chunks[1][0]
            chunk_size_lon = dz['level'].chunks[2][0]
            dz = dz.chunk({'date': chunk_size_date, 'lat': chunk_size_lat, 'lon': chunk_size_lon})
            dz.to_zarr(ZARR_PATH, mode='w', append_dim='date', group='anomaly')
        else:
            new_data_slice = dz.sel(date=next_month_date)
            # Append only the new data slice to the Zarr store
            expanded_data_slice = new_data_slice.expand_dims('date')
            expanded_data_slice.to_zarr(ZARR_PATH, mode='a', append_dim='date', group='anomaly')

        print("All work append to zarr done: ", ZARR_PATH, ' with re-chunking:', RECHUNK)
    
    else:
        print("No new data got since: ", last_date_in_zarr)


def main():
    if TEST_Download:
        filenames = download_noaa_data(TEST_date, DATA_PATH)
        if filenames is None or len(filenames) == 0:
            print("Error downloading NOAA data for date: ", TEST_date, ". Aborting update.")
            return        
        print("Download nc for the date: ", TEST_date, " and get file: ", filenames)

    else:    
        append_to_zarr()


if __name__ == "__main__":
    main()
