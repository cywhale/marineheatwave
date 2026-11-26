# Open API of Marine Heatwaves (MHWs), ODB

Marine heatwaves (MHWs) evaluated from 0.25-degree gridded NOAA OISST v2.1

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.8250030.svg)](https://doi.org/10.5281/zenodo.8250030)


#### ODB MHW open data official website

- https://eco.odb.ntu.edu.tw/pub/MHW

#### Reference

- Jacox, Alexander, Bograd, and Scott (2020), Thermal Displacement by Marine Heatwaves, Nature, 584, 82–86, doi:10.1038/s41586-020-2534-z",

#### Swagger API doc

- [ODB Marine Heatwaves (MHW) API manual/online try-out](https://api.odb.ntu.edu.tw/hub/swagger?node=odb_mhw_v1)

#### Usage: Endpoint /api/mhw

###### One-point MHWs without time-span limitation:

- e.g. https://eco.odb.ntu.edu.tw/api/mhw?lon0=135&lat0=15

###### Bonding-box <= 90 * 90 in degrees: 1-year time-span limitation for fetching MHWs

- e.g. https://eco.odb.ntu.edu.tw/api/mhw?lon0=135&lon1&=140&lat0=15&lat1=30&start=2021-01-01

###### Bonding-box > 90 * 90 in degrees: 1-month time-span limitation for fetching MHWs

- e.g. https://eco.odb.ntu.edu.tw/api/mhw?lon0=-180&lon1&=180&lat0=-90&lat1=90&start=2021-01-01


###### The parameter 'append' can append more data-variables (level, sst_anomaly, td) in data output

1. level: Severity of MHWs (ref: Hobday et al., 2018) 0 - 4: None, Moderate, Strong, Severe, Extreme

2. sst: Monthly SST

3. sst_anomaly: Monthly SST Anomalies

4. td: Thermal Displacement (ref: Jacox et al., 2020)
    
- e.g. https://eco.odb.ntu.edu.tw/api/mhw?lon0=135&lat0=15&append=sst_anomaly,level (default is level)

- The parameters lon1, lat1, start, end, and append are optional which can be none.

#### /api/mhw/csv

- The same usage as /api/mhw, but in CSV file response

#### Demo on <a href="https://eco.odb.ntu.edu.tw/pub/MHW/" target="_blank">ODB Marine Heatwaves</a>

<p align="center"><img src="https://eco.odb.ntu.edu.tw/pub/MHW/assets/202401_level.jpg" width=640 alt="Marine Heatwaves level 202401" /></p>

<p align="center"><em>Global marine heatwave levels at 0.25 degree grid for January 2024.</em></p>

#### Other resources

###### MHW evaluation source code

- https://github.com/yeh-tc/Python_project-Thermal_Displacement

###### QGIS tutorial for MHW WMS layer

- https://github.com/oceandatabank/MHW_QGIS

###### Introduction of MHW database by ODB (e-paper in pdf)

- https://www.odb.ntu.edu.tw/wp-content/uploads/2023/05/ODB_Newsletter_21.pdf

#### Development and deployment (maintainers)

- Release notes for the current milestone: [docs/release_v0.3.0_phase1.md](docs/release_v0.3.0_phase1.md); version history: [change_log.md](change_log.md)

###### Dependencies

- The **current** dependency set is the uv lock set in [requirements/p1-2/](requirements/p1-2/README.md): CPython 3.14.10, one venv shared by the API and the monthly updater.
    - `runtime.lock`: production (API + updater)
    - `dev.lock`: tests, installed together with `runtime.lock`
    - `tools.lock`: benchmark/check runners only
- `runtime.lock` contains `dask[array]` because the updater (`update_mhw.py`) reads data through xarray with Dask chunks. The API itself reads the zarr store without Dask and does **not** use a shared Dask scheduler.
- `requirements.txt` and `Pipfile` are the **legacy** pins used before 2026-11-25 (Python 3.11, zarr 2, dask/distributed 2025.3). They are kept unchanged for reference and are not used by production.

###### Production layout

- PM2 runs the app `mhwapi` from `conf/ecosystem.config.js`. `conf/start_app.sh` sources `conf/p12_env.sh` (venv path, `POLARS_FORCE_PKG=compat`) and execs gunicorn with two `uvicorn_worker.UvicornWorker` workers.
- The monthly cron job calls a wrapper script that lives on the production host, outside the repository. Its code is versioned as `conf/update_mhw_cron.sh.template`; the production copy differs only in the config block (host paths). The wrapper:
    - runs `conf/run_update_mhw.sh`, which runs `update_mhw.py` with the same venv;
    - restarts `mhwapi` only after a successful append;
    - exits non-zero on any failure (exit codes in its header).

###### Local development (legacy script)

- **Do not run `conf/simu.sh`.** It is the old local-development script, not migrated to the uv locks, and it is kept only as a record of the old commands. Run as a whole, on any machine including your own, it:
    - `pip install`s `requirements.txt` into whatever Python is current, changing that environment;
    - starts a local Dask scheduler and worker, which the current API no longer needs;
    - starts gunicorn with TLS files and `--reload`;
    - then `kill -9`s every process whose command line matches `mhw_app`, `dask scheduler` or `dask worker`, which can include unrelated processes;
    - and runs `pm2 start ./conf/ecosystem.config.js`, whose `cwd` is the production path.
- For a current local setup, create a venv from `requirements/p1-2/` (runtime + dev) and run `pytest` in `API/`.
- Moving `simu.sh` to the uv locks (without Dask) is proposed as a small separate change; not decided yet.
