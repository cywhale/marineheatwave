"""P1-2 V-c import probes, each run in its own fresh process (``python -I -B``).

    python -I -B vc_probe.py polars      # which polars runtime loads; CPU flags; warnings
    python -I -B vc_probe.py dbdriver    # psycopg 3 + the SQLAlchemy dialect; NO connection

Standalone (stdlib + the probed packages only), so it can run with the
interpreter of the environment under test.  ``dbdriver`` blocks every socket
connect before importing anything and only builds an engine object for a
TEST-NET placeholder URL (create_engine does not connect); it never reads a
.env and never imports update_mhw.py.
"""
import importlib.util
import json
import os
import socket
import sys
import warnings


def polars_probe() -> dict:
    out: dict = {"python": sys.version, "executable": sys.executable,
                 "env": {k: os.environ.get(k) for k in ("POLARS_FORCE_PKG", "POLARS_PREFER_PKG",
                                                         "POLARS_SKIP_CPU_CHECK")},
                 "installed_runtimes": sorted(n for n in ("_polars_runtime_compat", "_polars_runtime_64",
                                                          "_polars_runtime_32") if importlib.util.find_spec(n))}
    for n in out["installed_runtimes"]:
        mod = __import__(n)
        out.setdefault("build_feature_flags", {})[n] = getattr(mod, "BUILD_FEATURE_FLAGS", None)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        try:
            import polars as pl
            out["import_ok"] = True
            out["polars"] = pl.__version__
            plr = sys.modules["polars._plr"]
            out["plr_module"] = plr.__name__
            out["plr_file"] = getattr(plr, "__file__", None)
            out["df_check"] = pl.DataFrame({"a": [1.0, float("nan"), 3.0], "b": [1, None, 3]}).select(
                pl.col("a").is_nan().sum().alias("nan"), pl.col("b").is_null().sum().alias("null"),
                pl.col("b").sum().alias("sum")).to_dicts()
        except BaseException as exc:  # recorded, not raised
            out["import_ok"] = False
            out["error"] = f"{type(exc).__name__}: {exc}"
        out["warnings"] = [f"{x.category.__name__}: {str(x.message).splitlines()[0]}"
                           + (" | " + " ".join(str(x.message).split()) if "CPU" in str(x.message) else "")
                           for x in w]
    try:
        from polars import _cpu_check as cc
        flags = cc._read_cpu_flags()
        out["cpu_flags"] = flags
        for n, ff in (out.get("build_feature_flags") or {}).items():
            want = [f.lstrip("+") for f in (ff or "").split(",") if f and not f.startswith("-")]
            want = [f for f in want if f and f != "ctr-static"]  # same rule as polars._cpu_check
            out.setdefault("missing_for_runtime", {})[n] = sorted(f for f in want if f in flags and not flags[f])
    except BaseException as exc:
        out["cpu_flags_error"] = f"{type(exc).__name__}: {exc}"
    return out


def dbdriver_probe() -> dict:
    def refuse(*a, **k):
        raise OSError("network disabled by vc_probe")
    socket.socket.connect = refuse
    socket.socket.connect_ex = refuse
    socket.create_connection = refuse
    out: dict = {}
    for name in ("psycopg", "sqlalchemy", "dask", "netCDF4", "cftime", "psycopg2"):
        spec = importlib.util.find_spec(name)
        out.setdefault("found", {})[name] = bool(spec)
    import psycopg
    out["psycopg"] = psycopg.__version__
    out["psycopg_impl"] = psycopg.pq.__impl__
    out["libpq_version"] = psycopg.pq.version()
    import sqlalchemy
    from sqlalchemy import create_engine
    out["sqlalchemy"] = sqlalchemy.__version__
    url = "postgresql://placeholder:placeholder@192.0.2.1:5432/placeholder"  # TEST-NET-1, never contacted
    eng = create_engine(url)
    out["bare_url"] = {"dialect": eng.dialect.name, "driver": eng.dialect.driver,
                       "dbapi": getattr(eng.dialect.dbapi, "__name__", None), "class": type(eng.dialect).__name__}
    eng.dispose()
    for extra in ("postgresql+psycopg", "postgresql+psycopg2"):
        try:
            e = create_engine(url.replace("postgresql", extra, 1))
            out[extra] = {"driver": e.dialect.driver}
            e.dispose()
        except Exception as exc:
            out[extra] = {"error": f"{type(exc).__name__}: {exc}"}
    import dask
    import netCDF4
    out["dask"], out["netCDF4"] = dask.__version__, netCDF4.__version__
    out["connect_attempted"] = False  # nothing above opens a connection; a connect would have raised
    return out


if __name__ == "__main__":
    what = sys.argv[1]
    res = polars_probe() if what == "polars" else dbdriver_probe()
    print(json.dumps(res, indent=1, default=str))
