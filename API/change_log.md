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
    -- fix Basemap land mask wrongly affect data, need latlon=True in contourf

#### ver 0.2.6 Improve Swagger doc (small package upgrade/stable since last revision)

    -- improve Research example (detrended SST anomalies) code comments/n2
    -- add an ENSO timeseries example (moved from buoy project)
    -- try plot example to debug odbchat MHW cli patch

#### ver 0.2.7 Decouple shared Dask service, add resilient client/metadata consolidate/polars rtcompat

#### ver 0.2.8 Require TLS for shared Dask client; fix northern latitude boundary

    -- connect the MHW API to the shared Dask scheduler with its dedicated mutual-TLS client certificate
    -- map latitudes at or above 90 degrees to the dataset's northernmost grid point (89.875) to avoid empty point selections
    -- restrict the local-development Dask scheduler and worker in simu.sh to loopback, disable the worker dashboard, and fail closed if port 8786 is occupied or cannot be checked

#### ver 0.3.0 (draft, not tagged) Phase 1 core milestone: API performance and runtime environment

    -- read the zarr store directly with numpy-backed arrays (chunks=None) and only the requested variables; the API no longer connects to the shared Dask scheduler (src/dask_client_manager.py is kept as a no-op shim for old imports)
    -- CSV responses delete their temporary file after sending
    -- plots read only the plotted variable and draw on matplotlib Figure objects instead of global pyplot state
    -- heavy requests (plots, and JSON/CSV selections larger than MHW_HEAVY_BYTES, default 2 MiB) run one at a time per worker (MHW_HEAVY_SLOTS, default 1) in their own thread lane, so light queries are not starved; matplotlib drawing is serialised by a lock
    -- runtime: one uv venv shared by the API and the monthly updater, CPython 3.14.10, hash-locked dependencies in requirements/p1-2/ (zarr 3, xarray 2026.9, pandas 3.1, polars 2.2 with the compat runtime, SQLAlchemy 2.1 with psycopg 3, FastAPI 0.165 / Starlette 1.10, gunicorn 27 with uvicorn-worker)
    -- launchers: conf/p12_env.sh defines the venv and POLARS_FORCE_PKG=compat; conf/start_app.sh execs gunicorn (no --reload); new conf/run_update_mhw.sh runs update_mhw.py with the same venv; conf/ecosystem.config.js manages only mhwapi (explicit cwd and kill_timeout, pre_stop removed)
    -- acceptance on production (loopback): the 17 golden data cases are byte-identical; single-point median 0.121 s (about 3.0 s before Phase 1, 0.093 s with the P1-1 code on the old environment); the new environment is about 32 ms slower than the old one, within the provisional limits (median <= 0.3 s, paired difference CI upper bound <= +50 ms), which the owner has not yet formally adopted
    -- accepted output difference: the CSV Range 416 response now carries Content-Range "bytes */<size>" (Starlette 1.x); plot PNGs are checked by format and size, not by bytes
    -- not changed: response bodies of the data endpoints, update_mhw.py, and the OpenAPI info.version (1.0.0)
    -- known bug kept unchanged: a single bounding-box query crossing 180 degrees (e.g. lon0=170&lon1=-170) returns the wrong longitude range (the 340-degree complement) and a shortened period; work around it with two queries, one on each side of 180 degrees
    -- open: the first real monthly update on the new venv (2026-12-18), reboot recovery of the PM2 app, and the P1-3 monthly-update hardening; see docs/release_v0.3.0_phase1.md
