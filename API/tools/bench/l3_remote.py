"""L3 of spec 16 r2 §4: direct zarr reads on vm24, no HTTP, no distributed Dask.

Meant to be piped to the host interpreter over ssh; it only reads the store and
prints JSON lines to stdout (nothing is written on the host):

    ssh vm24 'nice -n 10 /home/odbadmin/.pyenv/versions/py311/bin/python - <ZARR> <ROUNDS>' < tools/bench/l3_remote.py

Modes (same selections as 0.2.8's process_mhw_data):
  auto_all   open_zarr(chunks='auto') + synchronous dask, all 4 variables  (0.2.8 minus distributed)
  auto_req   same, only the requested variables
  none_all   open_dataset(engine='zarr', chunks=None, cache=False), all 4 variables
  none_req   same, only the requested variables                           (P1-1a read path)
  none_req_pipeline  none_req + to_dataframe/dropna/polars/to_dicts/json.dumps (P1-1a CPU, no HTTP)
"""
import json
import sys
import time
import warnings

warnings.simplefilter("ignore")

import dask  # noqa: E402
import xarray as xr  # noqa: E402

ZARR = sys.argv[1]
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 else 7
dask.config.set(scheduler="synchronous")

FULL = slice("1982-01-01", "2026-10-01")
QUERIES = {
    # J01 / Q1: to_nearest_grid_point(135, 15) -> slice(135.0, 135.125)
    "Q1_point_135_15": dict(sel=dict(lon=slice(135.0, 135.125), lat=slice(15.0, 15.125), date=FULL), req=["level"]),
    # Q3: to_nearest_grid_point(135.3, 15.3) -> slice(135.25, 135.375)
    "Q3_point_135.3_15.3": dict(sel=dict(lon=slice(135.25, 135.375), lat=slice(15.25, 15.375), date=FULL), req=["level"]),
    # J09 / Q4: 10x10 deg, 10 years, sst_anomaly
    "Q4_bbox_10x10_10y": dict(sel=dict(lon=slice(130.0, 140.0), lat=slice(10.0, 20.0),
                                       date=slice("2010-01-01", "2019-12-01")), req=["sst_anomaly"]),
}


def emit(**kw):
    print(json.dumps(kw), flush=True)


def open_auto():
    return xr.open_zarr(ZARR, group="anomaly", chunks="auto", decode_times=True, consolidated=True)


def open_none():
    return xr.open_dataset(ZARR, engine="zarr", group="anomaly", chunks=None, cache=False,
                           decode_times=True, consolidated=True)


def pipeline(ds, variables):
    import polars as pl

    df = ds.to_dataframe().reset_index()
    mask = df["lon"] > 180
    df.loc[mask, "lon"] = df.loc[mask, "lon"] - 360
    df = df[["lon", "lat", "date"] + variables].dropna(how="all", subset=variables)
    pl_df = pl.from_pandas(df).with_columns(pl.col("date").dt.strftime("%Y-%m-%d"))
    return len(json.dumps(pl_df.to_dicts()))


for name, fn in (("open_auto", open_auto), ("open_none", open_none)):
    for i in range(3):
        t0 = time.perf_counter()
        fn().close()
        emit(kind="open", mode=name, rep=i, seconds=time.perf_counter() - t0)

stores = {"auto": open_auto(), "none": open_none()}
ALL = ["level", "sst", "sst_anomaly", "td"]
MODES = ["auto_all", "auto_req", "none_all", "none_req", "none_req_pipeline"]

for r in range(ROUNDS):
    order = MODES[r % len(MODES):] + MODES[: r % len(MODES)]  # rotate to spread warm-up bias
    for qname, q in QUERIES.items():
        for mode in order:
            store = stores["auto" if mode.startswith("auto") else "none"]
            variables = q["req"] if "_req" in mode else ALL
            t0 = time.perf_counter()
            sub = store[variables].sel(**q["sel"]).load()
            t_load = time.perf_counter() - t0
            extra = {}
            if mode == "none_req_pipeline":
                t1 = time.perf_counter()
                extra["json_bytes"] = pipeline(sub, sorted(variables))
                extra["pipeline_seconds"] = time.perf_counter() - t1
            emit(kind="read", round=r, query=qname, mode=mode, load_seconds=t_load,
                 total_seconds=time.perf_counter() - t0, nbytes=int(sub.nbytes), **extra)
