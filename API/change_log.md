ODB API of Marine Heatwaves
#### ver 0.1.0 /api/mhw and /api/mhw/csv for basic query

#### ver 0.1.1 deployment on eco.odb.ntu.edu.tw 20230725

#### ver 0.1.2 rewrite zarr to solve nan merged in diff rows problem/use pypolars

#### ver 0.1.3 Breaking fix bug conversion to [0,360]/pm2 pre-stop kill gunicorn processes

#### ver 0.1.4 Breaking add SST/Climatology plot API/Use Dask-scheduler/worker

#### ver 0.1.5 Add automation for update data from Postgres/Netcdf

#### ver 0.1.6 small fix input lon/lat within a 0.25-grid that cause no data when slicing

#### ver 0.2.0 officially release as open data, remove polar-circle MHW grids/renew all data

    -- change to pipenv, modify pm2 management script

#### ver 0.2.1 re-strcture project and use pipenv, fix update_mhw download bugs

    -- move to FastAPI lifespan async method/small package upgrade
    -- various package upgrade/n2/fix lack of __init__.py

#### ver 0.2.2 fix API gridded function inconsistent with original PostGIS database bug

    -- fix xr.open_mfdataset parallel=True cause core dump problem/append script using Geoserver api

#### ver 0.2.3 fix API month_mean (should not) mixed with area_mean method bug
#### ver 0.2.4 Try to shared dask worker leak by unique naming scheme for Dask tasks/major package upgrade (numpy v2)
#### ver 0.2.5 Add LONG_TERM_RANGE (10x10) and LONG_TERM_LIMIT (10yr) criteria/package upgrade

    -- add extra range/time criteria (10 x 10 x 10yr)  used by Hidy2 to evaluate range average of anomalies
    -- Research example of detrended SST anomalies
