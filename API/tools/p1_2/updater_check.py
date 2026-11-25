"""Isolated check of update_mhw.py against a writable COPY of the MHW store (P1-2 V-b).

No database, no NOAA, no network: the PostgreSQL query result and the downloaded
NOAA OISST daily files are fixtures built from the store's own last month.

    updater_check.py fixtures --store COPY --out FIXDIR        # db.parquet + oisst-*.nc + fixture.json
    updater_check.py snapshot --store COPY --out before.json   # per-file sha256 of the store
    updater_check.py run      --store COPY --fixtures FIXDIR --api API_DIR
    updater_check.py verify   --store COPY --fixtures FIXDIR --before before.json --before-store RO_COPY \
                              --api API_DIR --out report.json
    updater_check.py api-check --store COPY --fixtures FIXDIR --api API_DIR --out api.json  # any API env

Chained appends (P1-2, owner 2026-11-25: 538 -> 539 -> 540 months):

    updater_check.py meta-copy --store COPY --out META_BEFORE_N     # every .z* file, the pre-append metadata
    updater_check.py run ... --out run_N.json                       # timings, audit of .env / socket / writes
    updater_check.py diff --before snap_{N-1}.json --after snap_N.json --out diff_N.json
    updater_check.py verify ... --before-store META_BEFORE_N --earlier-fixtures FIX1 --no-api
    updater_check.py expect-points --fixtures FIX1 --fixtures FIX2 --out points.json  # for append_server
    updater_check.py check-result --file verify_N.json --expect ok            # 0 only when ok is true
    updater_check.py check-result --file neg.json --expect fail --expect-failed dates_last_is_fixture_month

Exit status (Codex r22: a failed check must not look like success):
    0  the command ran and its result is ok (run / verify / api-check), or
       check-result found the result as expected
    1  the command itself crashed (Python traceback; the JSON may be missing)
    2  usage error (argparse)
    3  run / verify / api-check: the JSON result was written completely and
       says ok=false (``failed_checks`` names the checks)
    4  check-result: the result is not as expected (or unreadable)

Fixture assumptions (not verified against the real DB / NOAA, which this round
may not contact):
* the DB query returns the full 0.25-degree grid (720 x 1440 rows per month),
  lat/lon as float64 on the store's grid, date as a date object, level as
  int64 without NULLs (land -1/0 as in the store), sst_anomaly and td as
  float64 with NaN where the store has NaN;
* an OISST v2.1 daily file has sst as int16, scale_factor 0.01, _FillValue
  -999, dims (time, zlev, lat, lon), time in "days since 1978-01-01 12:00:00",
  lat -89.875..89.875 and lon 0.125..359.875 (float32).  The real files also
  carry anom / err / ice, which update_mhw.py does not read.

``run`` never reads a .env: dotenv.load_dotenv is replaced before the import,
the DB settings are placeholders, create_engine / read_sql_query / requests.get
are stubs, and every socket connect raises.
"""
from __future__ import annotations

import argparse
import calendar
import contextlib
import hashlib
import io
import json
import os
import socket
import sys
import tempfile
import time
from datetime import date
from pathlib import Path

GROUP = "anomaly"
VARS = ("sst", "sst_anomaly", "level", "td")
if __package__ in (None, ""):  # run as a script (python -I tools/p1_2/updater_check.py): find tools/
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.p1_2.points import SAMPLES  # noqa: E402  (same five points everywhere)

EXIT_OK, EXIT_CRASH, EXIT_USAGE, EXIT_RESULT_FAILED, EXIT_UNEXPECTED = range(5)


def result_exit(rep: dict) -> int:
    """EXIT_OK only for a result whose ``ok`` is a real True."""
    return EXIT_OK if rep.get("ok") is True else EXIT_RESULT_FAILED


def with_checks(rep: dict, checks: dict) -> dict:
    """Record named checks; ``ok`` is true only when every check is."""
    rep["checks"] = {k: bool(v) for k, v in checks.items()}
    rep["failed_checks"] = sorted(k for k, v in rep["checks"].items() if not v)
    rep["ok"] = not rep["failed_checks"]
    return rep


def _month_after(d: date) -> date:
    return date(d.year + (d.month == 12), d.month % 12 + 1, 1)


def _open(store):
    import xarray as xr
    return xr.open_zarr(store, group=GROUP, decode_times=True, consolidated=True)


# ------------------------------------------------------------------ fixtures


def cmd_fixtures(a):
    import numpy as np
    import pandas as pd
    import xarray as xr

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    ds = _open(a.store)
    last = pd.Timestamp(ds["date"].values[-1]).date()
    month = _month_after(last)
    base = ds.isel(date=-1).load()
    lat = ds["lat"].values.astype("f8")
    lon = ds["lon"].values.astype("f8")

    # --- DB rows: last month shifted so the new month is recognisable
    LAT, LON = np.meshgrid(lat, lon, indexing="ij")
    sa = base["sst_anomaly"].values + 0.125
    td = base["td"].values + 1.0
    lv = base["level"].values.astype("int64")
    lv = np.where(np.isnan(sa), lv, np.clip(lv + 1, -1, 4))  # ocean levels moved by one (capped)
    df = pd.DataFrame({
        "lat": LAT.ravel(), "lon": LON.ravel(),
        "sst_anomaly": sa.ravel(), "level": lv.ravel(), "td": td.ravel(),
        "date": [month] * LAT.size,
    })
    df.to_parquet(out / "db.parquet", index=False)

    # --- NOAA daily files: base sst quantised to 0.01, +/- symmetric daily offsets
    ndays = calendar.monthrange(month.year, month.month)[1]
    # +0.25 so the new month differs from the last one everywhere (a copied old
    # month must fail the value check)
    sst0 = np.round(base["sst"].values.astype("f8"), 2) + 0.25
    offs = [round((d - (ndays + 1) / 2) * 0.01, 4) for d in range(1, ndays + 1)]
    t0 = np.datetime64("1978-01-01T12:00")
    for d, off in zip(range(1, ndays + 1), offs):
        day = np.datetime64(f"{month.year:04d}-{month.month:02d}-{d:02d}T12:00")
        v = (sst0 + off).astype("f4")  # packed to int16 * 0.01 by the encoding below
        da = xr.DataArray(v[None, None], dims=("time", "zlev", "lat", "lon"),
                          coords={"time": [day], "zlev": np.array([0.0], "f4"),
                                  "lat": lat.astype("f4"), "lon": lon.astype("f4")},
                          name="sst", attrs={"long_name": "Daily sea surface temperature", "units": "Celsius"})
        dsd = da.to_dataset()
        dsd["time"].encoding.update(units="days since 1978-01-01 12:00:00", dtype="f4")
        dsd["sst"].encoding.update(dtype="i2", scale_factor=0.01, add_offset=0.0, _FillValue=np.int16(-999),
                                   zlib=True, complevel=4)
        name = f"oisst-avhrr-v02r01.{month.year:04d}{month.month:02d}{d:02d}.nc"
        dsd.to_netcdf(out / name, format="NETCDF4")
    expected_sst = np.where(np.isnan(sst0), np.nan, sst0 + float(np.mean(offs)))
    np.save(out / "expected_sst.npy", expected_sst.astype("f4"))
    np.save(out / "base_lastmonth.npy", np.stack([base[v].values.astype("f8") for v in VARS]))
    meta = {"last_month_in_store": str(last), "new_month": str(month), "ndays": ndays,
            "db_rows": len(df), "db_dtypes": {k: str(v) for k, v in df.dtypes.items()},
            "noaa_files": sorted(p.name for p in out.glob("oisst-*.nc"))}
    (out / "fixture.json").write_text(json.dumps(meta, indent=1))
    print(json.dumps({k: v for k, v in meta.items() if k != "noaa_files"}, indent=1))


# ------------------------------------------------------------------ snapshot


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def snapshot(store) -> dict:
    """key -> [sha256, size] for every file of the store."""
    root = Path(store)
    return {str(p.relative_to(root)): [_sha(p), p.stat().st_size] for p in sorted(root.rglob("*")) if p.is_file()}


def _digest(v):
    return v if isinstance(v, str) else v[0]  # V-b snapshots hold the sha only


def _size(v):
    return None if isinstance(v, str) else v[1]


def manifest_diff(before: dict, after: dict) -> dict:
    """Added / changed / removed keys with sizes and hashes, grouped by array."""
    def group(k):
        parts = k.split("/")
        return "/".join(parts[:-1]) or "."
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    changed = sorted(k for k in set(after) & set(before) if _digest(after[k]) != _digest(before[k]))
    out: dict = {"n_before": len(before), "n_after": len(after), "added": {}, "changed": {}, "removed": {},
                 "unchanged": len(set(after) & set(before)) - len(changed)}
    for name, keys, src in (("added", added, after), ("changed", changed, after), ("removed", removed, before)):
        for k in keys:
            e = {"key": k, "sha256": _digest(src[k]), "size": _size(src[k])}
            if name == "changed":
                e.update(sha256_before=_digest(before[k]), size_before=_size(before[k]))
            out[name].setdefault(group(k), []).append(e)
    out["bytes_written"] = sum((e["size"] or 0) for name in ("added", "changed") for g in out[name].values() for e in g)
    out["summary"] = {name: {g: len(v) for g, v in out[name].items()} for name in ("added", "changed", "removed")}
    return out


def cmd_diff(a):
    d = manifest_diff(json.loads(Path(a.before).read_text()), json.loads(Path(a.after).read_text()))
    Path(a.out).write_text(json.dumps(d, indent=1))
    print(json.dumps({"summary": d["summary"], "bytes_written": d["bytes_written"], "unchanged": d["unchanged"]},
                     indent=1))


def cmd_meta_copy(a):
    """Copy every metadata file (.zarray / .zattrs / .zgroup / .zmetadata) of the
    store, keeping the layout: the pre-append "before store" for metadata_problems."""
    import shutil
    root, out = Path(a.store), Path(a.out)
    out.mkdir(parents=True, exist_ok=False)
    n = 0
    for p in sorted(root.rglob(".z*")):
        if p.is_file():
            q = out / p.relative_to(root)
            q.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, q)
            n += 1
    print(f"{n} metadata files")


def cmd_snapshot(a):
    snap = snapshot(a.store)
    Path(a.out).write_text(json.dumps(snap))
    print(f"{len(snap)} files")


# ----------------------------------------------------------------------- run


class _NoNetwork(OSError):
    pass


def _block_network():
    def refuse(*a, **k):
        _NET_ATTEMPTS.append("blocked")
        raise _NoNetwork("network disabled in updater_check")
    socket.socket.connect = refuse
    socket.socket.connect_ex = refuse
    socket.create_connection = refuse
    socket.getaddrinfo = refuse


_NET_ATTEMPTS: list[str] = []


def _audit_recorder(store: str, allowed_write_roots: list[str]):
    """sys.addaudithook: every socket / subprocess / urllib event, every open of a
    file named .env, and every write-mode open or rename outside the allowed roots."""
    events: dict = {"socket": [], "subprocess": [], "urllib": [], "dotenv_file_open": [], "write_outside": [],
                    "write_opens_in_store": 0, "hook_errors": []}
    roots = [os.path.realpath(r) for r in allowed_write_roots]

    def inside(path) -> bool:
        rp = os.path.realpath(os.fsdecode(path))
        return any(rp == r or rp.startswith(r + os.sep) for r in roots)

    def hook(event, args):
        try:
            _hook(event, args)
        except Exception as exc:  # an audit hook must never break the audited call
            if len(events["hook_errors"]) < 20:
                events["hook_errors"].append(f"{event}: {type(exc).__name__}: {exc}")

    def _hook(event, args):
        if event.startswith("socket.") and event not in ("socket.__new__",):
            events["socket"].append(event)
        elif event.startswith("subprocess.") or event == "os.system" or event.startswith("os.exec") \
                or event == "os.posix_spawn":
            events["subprocess"].append(event)
        elif event.startswith("urllib."):
            events["urllib"].append(event)
        elif event == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
            path, mode, flags = args[0], args[1], args[2] if len(args) > 2 else 0
            name = os.path.basename(os.fsdecode(path))
            if name == ".env":
                events["dotenv_file_open"].append(os.fsdecode(path))
            writing = (isinstance(mode, str) and any(c in mode for c in "wax+")) or \
                (isinstance(flags, int) and flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT))
            if writing:
                if inside(path):
                    if os.path.realpath(os.fsdecode(path)).startswith(os.path.realpath(store)):
                        events["write_opens_in_store"] += 1
                elif len(events["write_outside"]) < 200:
                    events["write_outside"].append(os.fsdecode(path))
        elif event in ("os.rename", "os.replace", "os.remove", "os.unlink", "shutil.rmtree", "os.rmdir"):
            for x in args[:2]:
                if isinstance(x, (str, bytes, os.PathLike)) and not inside(x) and len(events["write_outside"]) < 200:
                    events["write_outside"].append(f"{event}:{os.fsdecode(x)}")
    return events, hook


def cmd_run(a):
    import pandas as pd

    fix = Path(a.fixtures).resolve()
    store = str(Path(a.store).resolve())
    out_path = Path(a.out).resolve() if a.out else None  # the run changes cwd
    _block_network()

    calls = {"load_dotenv": 0, "create_engine": [], "read_sql_query": 0, "requests_get": 0}
    import dotenv

    def no_dotenv(*x, **k):  # never read a .env
        calls["load_dotenv"] += 1
        return False
    dotenv.load_dotenv = no_dotenv
    for k in ("DBUSER", "DBPASS", "DBHOST", "DBPORT", "DBNAME", "MHWTABLE"):
        os.environ[k] = "fixture"

    queries = []
    db = pd.read_parquet(fix / "db.parquet")

    class Engine:
        def dispose(self):
            queries.append("dispose")

    def fake_read_sql_query(q, engine, *x, **k):
        calls["read_sql_query"] += 1
        queries.append(q)
        assert isinstance(engine, Engine)
        return db.copy()

    import requests
    import sqlalchemy

    def no_requests(*x, **k):
        calls["requests_get"] += 1
        raise _NoNetwork("requests.get called: a fixture file is missing")

    def fake_engine(url, *x, **k):
        calls["create_engine"].append(str(url).split("://", 1)[0])  # the scheme only; never the URL
        return Engine()

    requests.get = no_requests
    sqlalchemy.create_engine = fake_engine
    pd.read_sql_query = fake_read_sql_query

    sys.path.insert(0, str(Path(a.api).resolve()))
    import update_mhw as upd
    upd.create_engine = sqlalchemy.create_engine  # bound at import ("from sqlalchemy import create_engine")
    upd.ZARR_PATH = store
    upd.DATA_PATH = str(fix) + "/"

    # separate timings of the append write and of the explicit consolidation
    import xarray as xr
    # to_zarr includes xarray's own consolidation of the store (consolidated=None
    # on a consolidated v2 store); update_mhw.py then consolidates once more
    timings = {"to_zarr": [], "consolidate_metadata": []}
    orig_to_zarr, orig_consolidate = xr.Dataset.to_zarr, upd.zarr.consolidate_metadata
    active = {"to_zarr": 0}

    def timed_to_zarr(self, *x, **k):
        t0 = time.perf_counter()
        active["to_zarr"] += 1
        try:
            return orig_to_zarr(self, *x, **k)
        finally:
            active["to_zarr"] -= 1
            timings["to_zarr"].append({"seconds": time.perf_counter() - t0, "mode": k.get("mode"),
                                       "append_dim": k.get("append_dim")})

    def timed_consolidate(*x, **k):
        t0 = time.perf_counter()
        inside = active["to_zarr"] > 0
        try:
            return orig_consolidate(*x, **k)
        finally:
            timings["consolidate_metadata"].append({"seconds": time.perf_counter() - t0,
                                                    "caller": "xarray to_zarr (inside)" if inside
                                                    else "update_mhw.py (explicit)"})
    xr.Dataset.to_zarr = timed_to_zarr
    upd.zarr.consolidate_metadata = timed_consolidate

    error = None
    with tempfile.TemporaryDirectory() as cwd:  # relative paths in the module must not matter
        # the allowed write roots, as resolved paths (doc 28 r6: TMPDIR is the
        # isolated work-dir tmp/ in append_chain.sh, so gettempdir() is not /tmp)
        roots = [store, cwd, tempfile.gettempdir()]
        events, hook = _audit_recorder(store, roots)
        sys.addaudithook(hook)  # cannot be removed; this process ends after the run
        start_dir = os.getcwd()
        os.chdir(cwd)
        t = time.perf_counter()
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                upd.append_to_zarr()
        except Exception as exc:  # recorded; the JSON is still written and the exit status is 3
            import traceback
            error = {"type": type(exc).__name__, "message": str(exc)[:2000],
                     "traceback_tail": traceback.format_exc()[-4000:]}
        finally:
            dt = time.perf_counter() - t
            os.chdir(start_dir)  # never stay in a directory that is about to be removed
    out = buf.getvalue()
    print(out)
    rep = {"seconds_total": round(dt, 3), "timings": timings, "queries": queries, "calls": calls,
           "network_attempts_blocked": len(_NET_ATTEMPTS), "audit": events,
           "allowed_write_roots": [os.path.realpath(r) for r in roots],
           "tempfile_gettempdir": os.path.realpath(tempfile.gettempdir()),
           "error": error, "completed": error is None and "All work append to zarr done" in out,
           "stdout_tail": out[-2000:]}
    rep["no_external_access"] = (not events["socket"] and not events["urllib"] and not events["subprocess"]
                                 and not events["dotenv_file_open"] and not _NET_ATTEMPTS
                                 and calls["requests_get"] == 0 and not events["write_outside"]
                                 and not events["hook_errors"])
    with_checks(rep, run_checks(rep))
    if out_path:
        out_path.write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps({k: v for k, v in rep.items() if k not in ("stdout_tail",)}, indent=1, default=str))
    return result_exit(rep)


def run_checks(rep: dict) -> dict:
    """Named checks of a ``run`` report (pure, for the tests)."""
    ev = rep["audit"]
    return {"no_error": rep.get("error") is None,
            "completed": rep["completed"],
            "no_external_access": rep["no_external_access"],
            "no_dotenv_open": not ev["dotenv_file_open"],
            "no_socket_urllib_subprocess": not (ev["socket"] or ev["urllib"] or ev["subprocess"]),
            "no_network_attempt": rep["network_attempts_blocked"] == 0 and rep["calls"]["requests_get"] == 0,
            "no_write_outside": not ev["write_outside"],
            "no_audit_hook_errors": not ev["hook_errors"]}


# -------------------------------------------------------------------- verify


def cmd_verify(a):
    rep = _verify_report(a)
    f = rep["files"]
    Path(a.out).write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps({"ok": rep["ok"], "failed_checks": rep["failed_checks"],
                      "files": {k: (x if not isinstance(x, list) or len(x) < 20 else f"{len(x)} items")
                                for k, x in f.items()},
                      "dates": rep["dates"], "values": rep["values"], "api": rep["api"]}, indent=1, default=str))
    return result_exit(rep)


def _verify_report(a) -> dict:
    import numpy as np
    import pandas as pd
    import zarr

    fix = Path(a.fixtures)
    meta = json.loads((fix / "fixture.json").read_text())
    before = {k: _digest(v) for k, v in json.loads(Path(a.before).read_text()).items()}
    after = {k: _digest(v) for k, v in snapshot(a.store).items()}
    rep: dict = {"fixture": meta}

    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    changed = sorted(k for k in set(after) & set(before) if after[k] != before[k])
    n_old = None

    ds = _open(a.store)
    dates = pd.DatetimeIndex(ds["date"].values)
    n_old = len(dates) - 1
    date_chunk = ds["date"].encoding.get("chunks", (501,))[0]
    data_chunk = ds["sst"].encoding["chunks"][0]
    new_ci = n_old // data_chunk
    new_chunks = {f"{GROUP}/{v}/{new_ci}.{i}.{j}" for v in VARS for i in range(4) for j in range(4)}
    date_file = f"{GROUP}/date/{n_old // date_chunk}"  # added when the new date starts a chunk, else changed
    allowed_added = new_chunks | {date_file}
    allowed_changed = new_chunks | {date_file}
    data_added = sorted(f for f in set(added) - allowed_added if not is_metadata_file(f))
    data_changed = sorted(f for f in set(changed) - allowed_changed if not is_metadata_file(f))
    meta_files = sorted(f for f in set(added) | set(changed) if is_metadata_file(f))
    meta_diff = _metadata_diff(a.before_store, a.store, meta_files)
    # the month count before this append, from the pre-append metadata (not from
    # the store under test, which would make "n == n_old + 1" a tautology).
    # --before-store is required (Codex r22 follow-up): without it there is no
    # independent baseline, so verify cannot run at all.
    zm0 = json.loads((Path(a.before_store) / ".zmetadata").read_text())["metadata"]
    n_before = zm0.get(f"{GROUP}/date/.zarray", {}).get("shape", [None])[0]
    rep["n_before"] = n_before
    rep["not_verified"] = (["api: --no-api (the API is checked separately, e.g. append_server)"]
                           if a.no_api else [])
    rep["files"] = {
        "before": len(before), "after": len(after), "added": len(added), "removed": removed,
        "changed": changed, "unexpected_data_added": data_added, "unexpected_data_changed": data_changed,
        "metadata_files": meta_files,
        "unchanged_preexisting": len(before) - len(changed) - len(removed),
        "metadata_diff": meta_diff,
    }
    # every metadata document is checked against the explicit rules in
    # metadata_problems(): only the appended arrays may change shape, and only
    # by one along date; other changes must be formatting or in KNOWN_REWRITES
    unexplained = data_added + data_changed + metadata_problems(a.before_store, a.store, n_old)
    rep["files"]["unexplained"] = unexplained

    # structure
    rep["dates"] = {"n": len(dates), "last": str(dates[-1].date()), "prev": str(dates[-2].date()),
                    "monotonic_unique": bool(dates.is_monotonic_increasing and dates.is_unique),
                    "consecutive_months": all(_month_after(dates[i].date()) == dates[i + 1].date()
                                              and dates[i].day == 1 for i in range(len(dates) - 1)),
                    "first": str(dates[0].date())}
    coord_files = sorted(k for k in set(after) | set(before) if k.split("/")[:2] in (["anomaly", "lat"], ["anomaly", "lon"]))
    rep["lat_lon_files_unchanged"] = all(after.get(k) == before.get(k) for k in coord_files) and bool(coord_files)
    rep["dims"] = {v: dict(ds[v].sizes) for v in VARS}
    rep["dtypes"] = {v: str(ds[v].dtype) for v in VARS}
    rep["chunks"] = {v: list(ds[v].encoding.get("chunks", ())) for v in VARS}
    g = zarr.open_group(a.store, mode="r")
    rep["zarr_shapes"] = {v: list(g[f"{GROUP}/{v}"].shape) for v in (*VARS, "date")}
    rep["consolidated_shapes"] = {}
    zm = json.loads((Path(a.store) / ".zmetadata").read_text())["metadata"]
    for v in (*VARS, "date"):
        rep["consolidated_shapes"][v] = zm[f"{GROUP}/{v}/.zarray"]["shape"]

    # values: new month vs fixture, previous month vs the pre-update copy
    new = ds.isel(date=-1).load()
    prev = ds.isel(date=-2).load()
    base = np.load(fix / "base_lastmonth.npy")
    exp = _expected(fix)
    vals = {}
    for k, v in enumerate(VARS):
        got = new[v].values.astype("f8")
        e = exp[v].astype("f8")
        tol = 0.006 if v == "sst" else 0.0
        same_nan = bool(np.array_equal(np.isnan(got), np.isnan(e)))
        diff = float(np.nanmax(np.abs(got - e))) if np.isfinite(got).any() else None
        prev_same = bool(np.array_equal(prev[v].values.astype("f8"), base[k], equal_nan=True))
        vals[v] = {"nan_pattern_equal": same_nan, "max_abs_diff": diff, "within_tol": diff is not None and diff <= tol,
                   "prev_month_unchanged": prev_same, "dtype": str(new[v].dtype)}
    rep["values"] = vals

    # earlier appended months (chained appends): still equal to their fixtures
    rep["earlier_months"] = {}
    for ef in a.earlier_fixtures or []:
        emeta = json.loads((Path(ef) / "fixture.json").read_text())
        eexp = _expected(Path(ef))
        mon = ds.sel(date=emeta["new_month"]).load()
        rep["earlier_months"][emeta["new_month"]] = {
            v: {"nan_pattern_equal": bool(np.array_equal(np.isnan(mon[v].values.astype("f8")), np.isnan(eexp[v]))),
                "max_abs_diff": float(np.nanmax(np.abs(mon[v].values.astype("f8") - eexp[v]))),
                "within_tol": float(np.nanmax(np.abs(mon[v].values.astype("f8") - eexp[v])))
                <= (0.006 if v == "sst" else 0.0)} for v in VARS}

    # API on the updated copy (its own startup code, relative data/mhw.zarr)
    rep["api"] = {"ok": True, "skipped": True} if a.no_api else \
        _api_check(a.api, a.store, [*(a.earlier_fixtures or []), a.fixtures])

    return with_checks(rep, verify_checks(rep, meta, n_before))


DTYPES = {"sst": "float32", "sst_anomaly": "float64", "level": "int64", "td": "float64"}


def verify_checks(rep: dict, meta: dict, n_before: int | None) -> dict:
    """Named checks of a ``verify`` report (pure, for the tests).  n_before is the
    month count of the pre-append metadata (--before-store, required); None when
    that metadata has no date shape, which fails date_count_is_before_plus_one."""
    f, v, d = rep["files"], rep["values"], rep["dates"]
    n = d["n"]
    c = {"files_none_removed": not f["removed"],
         "files_all_explained": not f["unexplained"],
         "dates_last_is_fixture_month": d["last"] == meta["new_month"],
         "dates_prev_is_fixture_base": d["prev"] == meta["last_month_in_store"],
         "dates_monotonic_unique": d["monotonic_unique"],
         "dates_consecutive_months": d["consecutive_months"],
         "lat_lon_files_unchanged": rep["lat_lon_files_unchanged"],
         "earlier_months_unchanged": all(x["nan_pattern_equal"] and x["within_tol"]
                                         for m in rep["earlier_months"].values() for x in m.values()),
         "new_month_values": all(x["nan_pattern_equal"] and x["within_tol"] for x in v.values()),
         "prev_month_unchanged": all(x["prev_month_unchanged"] for x in v.values()),
         "dtypes": rep["dtypes"] == DTYPES,
         "zarr_shapes": all(s[0] == n for s in rep["zarr_shapes"].values()),
         "consolidated_shapes": all(s[0] == n for s in rep["consolidated_shapes"].values()),
         "api": rep["api"].get("ok") is True,
         # always present: the baseline is required, its date shape must exist and be one less
         "date_count_is_before_plus_one": isinstance(n_before, int) and n == n_before + 1}
    return c


def _metadata_diff(before_store, store, files):
    """Compare metadata documents of the pre-update store and the updated copy.

    Per file: "format_only" (same JSON, different whitespace / key order),
    "added_empty_attrs" (a new ``{}`` .zattrs), or the differing top-level keys.
    For .zmetadata: per entry, the keys whose values differ.
    """
    out = {"files": {}, "zmetadata_semantic": {}}
    b0, b1 = Path(before_store), Path(store)
    for f in files:
        p0, p1 = b0 / f, b1 / f
        if not p1.name.startswith(".z"):
            continue
        j1 = json.loads(p1.read_text())
        if not p0.exists():
            out["files"][f] = "added_empty_attrs" if (p1.name == ".zattrs" and j1 == {}) else "added"
            continue
        j0 = json.loads(p0.read_text())
        if f == ".zmetadata":
            m0, m1 = j0["metadata"], j1["metadata"]
            for k in sorted(set(m0) | set(m1)):
                if k not in m0:
                    out["zmetadata_semantic"].setdefault("+file", []).append(k)
                elif k not in m1:
                    out["zmetadata_semantic"].setdefault("-file", []).append(k)
                elif m0[k] != m1[k]:
                    out["zmetadata_semantic"][k + "@.zmetadata"] = sorted(
                        x for x in set(m0[k]) | set(m1[k]) if m0[k].get(x) != m1[k].get(x))
            out["files"][f] = "semantic" if out["zmetadata_semantic"] else "format_only"
        else:
            out["files"][f] = "format_only" if j0 == j1 else sorted(
                x for x in set(j0) | set(j1) if j0.get(x) != j1.get(x))
    return out


METADATA_FILES = (".zarray", ".zattrs", ".zgroup", ".zmetadata")
APPENDED = tuple(f"{GROUP}/{v}" for v in (*VARS, "date"))  # arrays update_mhw.py extends along date
# metadata zarr 3 adds when it rewrites this v2 store (doc 29 r3 §7, V-b); each
# rule: (where, document, key) -> predicate(value before, value after)
def _is_v2_inline_envelope(a) -> bool:
    """The exact value zarr 3.4.2 wrote in V-b: {"kind": "inline", "must_understand":
    false, "metadata": {}} (no other keys, a real JSON false, an empty dict)."""
    return (type(a) is dict and set(a) == {"kind", "must_understand", "metadata"}
            and a["kind"] == "inline" and a["must_understand"] is False
            and type(a["metadata"]) is dict and a["metadata"] == {})


KNOWN_REWRITES = {
    ("zmetadata", f"{GROUP}/.zgroup", "consolidated_metadata"): lambda b, a: b is None and _is_v2_inline_envelope(a),
    ("zmetadata", f"{GROUP}/lat/.zarray", "dimension_separator"): lambda b, a: b is None and a == ".",
    ("zmetadata", f"{GROUP}/lon/.zarray", "dimension_separator"): lambda b, a: b is None and a == ".",
}
_MISSING = object()


def is_metadata_file(path: str) -> bool:
    return path.rsplit("/", 1)[-1] in METADATA_FILES


def _doc_problems(where: str, name: str, j0: dict, j1: dict, n_old: int) -> list[str]:
    """Problems of one metadata document (an on-disk file or a .zmetadata entry)."""
    probs = []
    if name.endswith("/.zarray"):
        arr = name[: -len("/.zarray")]
        s0, s1 = j0.get("shape"), j1.get("shape")
        if arr in APPENDED:
            want = [n_old + 1, *s0[1:]] if isinstance(s0, list) and s0 and s0[0] == n_old else None
            if want is None:
                probs.append(f"{where}:{name}: shape before {s0}, expected first dimension {n_old}")
            elif s1 != want:
                probs.append(f"{where}:{name}: shape {s0} -> {s1}, expected {want}")
        elif s0 != s1:
            probs.append(f"{where}:{name}: shape changed {s0} -> {s1} (not an appended array)")
    for key in sorted(set(j0) | set(j1)):
        if key == "shape" and name.endswith("/.zarray"):
            continue
        b, a = j0.get(key, _MISSING), j1.get(key, _MISSING)
        if b == a:
            continue
        rule = KNOWN_REWRITES.get((where, name, key))
        if rule is None or not rule(None if b is _MISSING else b, None if a is _MISSING else a):
            probs.append(f"{where}:{name}: key {key!r} changed {'<absent>' if b is _MISSING else b!r} -> "
                         f"{'<absent>' if a is _MISSING else a!r}")
    return probs


def metadata_problems(before_store, store, n_old: int) -> list[str]:
    """Every metadata difference between the pre-update store and the updated copy
    that is not explicitly expected.  Checks all .zmetadata entries and every
    on-disk metadata file of either store, whether or not its bytes changed."""
    b0, b1 = Path(before_store), Path(store)
    probs: list[str] = []
    z0 = json.loads((b0 / ".zmetadata").read_text())
    z1 = json.loads((b1 / ".zmetadata").read_text())
    if {k: v for k, v in z0.items() if k != "metadata"} != {k: v for k, v in z1.items() if k != "metadata"}:
        probs.append(".zmetadata: top-level keys other than 'metadata' changed")
    m0, m1 = z0["metadata"], z1["metadata"]
    for k in sorted(set(m0) - set(m1)):
        probs.append(f"zmetadata:{k}: entry removed")
    for k in sorted(set(m1) - set(m0)):
        if not (k == ".zattrs" and m1[k] == {}):
            probs.append(f"zmetadata:{k}: entry added")
    for k in sorted(set(m0) & set(m1)):
        probs += _doc_problems("zmetadata", k, m0[k], m1[k], n_old)
    for arr in APPENDED:  # also an appended array whose entry did not change at all
        if f"{arr}/.zarray" not in m0:
            probs.append(f"zmetadata:{arr}/.zarray: missing")
    files0 = {str(p.relative_to(b0)) for p in b0.rglob(".z*") if p.is_file() and p.name != ".zmetadata"}
    files1 = {str(p.relative_to(b1)) for p in b1.rglob(".z*") if p.is_file() and p.name != ".zmetadata"}
    for f in sorted(files0 - files1):
        probs.append(f"file:{f}: removed")
    for f in sorted(files1 - files0):
        if not (f == ".zattrs" and json.loads((b1 / f).read_text()) == {}):
            probs.append(f"file:{f}: added")
    for f in sorted(files0 & files1):
        probs += _doc_problems("file", f, json.loads((b0 / f).read_text()), json.loads((b1 / f).read_text()), n_old)
    for arr in APPENDED:
        if f"{arr}/.zarray" not in files0:
            probs.append(f"file:{arr}/.zarray: missing")
    return probs


def _expected(fix: Path):
    import numpy as np
    import pandas as pd
    db = pd.read_parquet(fix / "db.parquet")
    return {
        "sst": np.load(fix / "expected_sst.npy").astype("f8"),
        "sst_anomaly": db["sst_anomaly"].to_numpy().reshape(720, 1440),
        "level": db["level"].to_numpy().reshape(720, 1440),
        "td": db["td"].to_numpy().reshape(720, 1440),
    }


def expected_points(fixture_dirs) -> dict:
    """"lat,lon" -> {month: {field: value or None}} for every fixture month."""
    import numpy as np
    from tools.p1_2.points import FIELDS, grid_index
    out: dict = {}
    for fd in fixture_dirs:
        fd = Path(fd)
        month = json.loads((fd / "fixture.json").read_text())["new_month"]
        exp = _expected(fd)
        for lat, lon in SAMPLES:
            i, j = grid_index(lat, lon)
            w = {k: (None if np.isnan(exp[k][i, j]) else float(exp[k][i, j])) for k in FIELDS}
            out.setdefault(f"{lat},{lon}", {})[month] = w
    return out


def cmd_expect_points(a):
    pts = expected_points(a.fixtures)
    Path(a.out).write_text(json.dumps(pts, indent=1))
    print(f"{len(pts)} points x {len(next(iter(pts.values())))} months")


def cmd_api_check(a):
    """The API on the updated copy only; runs in any API environment (e.g. the
    zarr 2 / P1-1 one, to check that a store written by zarr 3 stays readable)."""
    import zarr
    res = _api_check(a.api, a.store, [*(a.earlier_fixtures or []), a.fixtures])
    res["zarr"] = zarr.__version__
    Path(a.out).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({"ok": res["ok"], "failed_checks": res["failed_checks"], "zarr": res["zarr"],
                      "plot": res["plot"], "cases": [(c["status"], c.get("match")) for c in res["cases"]]}))
    return result_exit(res)


def result_consistency(doc) -> list[str]:
    """Problems with the status fields of a result written by with_checks():
    ``ok`` a JSON boolean; ``checks`` a non-empty dict of name -> JSON boolean;
    ``failed_checks`` a list of names equal to the sorted false checks; ok true
    exactly when failed_checks is empty."""
    if not isinstance(doc, dict):
        return ["result is not a JSON object"]
    probs = []
    ok, checks, failed = doc.get("ok", _MISSING), doc.get("checks", _MISSING), doc.get("failed_checks", _MISSING)
    if type(ok) is not bool:
        probs.append("'ok' missing or not a boolean")
    if checks is _MISSING or type(checks) is not dict:
        probs.append("'checks' missing or not an object")
        checks = None
    elif not checks:
        probs.append("'checks' is empty")
    elif not all(isinstance(k, str) and type(v) is bool for k, v in checks.items()):
        probs.append("'checks' has a non-boolean value")
        checks = None
    if failed is _MISSING or type(failed) is not list or not all(isinstance(x, str) for x in failed):
        probs.append("'failed_checks' missing or not a list of names")
        failed = None
    if checks is not None and failed is not None and failed != sorted(k for k, v in checks.items() if not v):
        probs.append("'failed_checks' does not match the false entries of 'checks'")
    if type(ok) is bool and failed is not None and ok != (not failed):
        probs.append("'ok' contradicts 'failed_checks'")
    if type(ok) is bool and checks and ok != all(checks.values()):
        probs.append("'ok' contradicts 'checks'")
    return probs


def cmd_check_result(a):
    """Exit 0 only when the JSON result is as expected and its status fields agree
    (result_consistency).  --expect ok: ok true, every named check true,
    failed_checks empty.  --expect fail: ok false, failed_checks non-empty, and
    every --expect-failed name among them (a negative control must fail, and for
    the reason it was built for)."""
    verdict = {"file": a.file, "expect": a.expect, "expect_failed": a.expect_failed or []}
    try:
        doc = json.loads(Path(a.file).read_text())
    except (OSError, ValueError) as exc:
        verdict.update(as_expected=False, reason=f"unreadable: {type(exc).__name__}")
        print(json.dumps(verdict))
        return EXIT_UNEXPECTED
    probs = result_consistency(doc)
    ok = doc.get("ok") if isinstance(doc, dict) else None
    failed = doc.get("failed_checks") if isinstance(doc, dict) else None
    verdict.update(ok=ok, failed_checks=failed, inconsistent=probs)
    if probs:
        verdict.update(as_expected=False, reason="inconsistent status fields")
    elif a.expect == "ok":
        verdict.update(as_expected=ok is True and failed == [],
                       reason="ok" if ok else f"failed: {failed}")
    else:
        missing = [n for n in (a.expect_failed or []) if n not in failed]
        verdict.update(as_expected=ok is False and bool(failed) and not missing,
                       reason="ok is true" if ok else (f"not among failed_checks: {missing}" if missing
                                                        else "failed as expected"))
    print(json.dumps(verdict))
    return EXIT_OK if verdict["as_expected"] else EXIT_UNEXPECTED


def _api_check(api_dir, store, fixture_dirs):
    """In-process (TestClient) API on the store: every sample point returns exactly
    one row per fixture month with the fixture values; plot is a PNG."""
    from tools.p1_2.points import compare, url
    pts = expected_points(fixture_dirs)
    months = sorted(next(iter(pts.values())))
    api = Path(api_dir).resolve()
    res = {"cases": [], "months": months}
    with tempfile.TemporaryDirectory() as root:
        (Path(root) / "data").mkdir()
        os.symlink(Path(store).resolve(), Path(root) / "data" / "mhw.zarr")
        old = os.getcwd()
        os.chdir(root)
        sys.path.insert(0, str(api))
        try:
            import mhw_app
            from fastapi.testclient import TestClient
            with TestClient(mhw_app.app) as c:
                for lat, lon in SAMPLES:
                    u = url("json", lat, lon, months[0], months[-1])
                    r = c.get(u)
                    case = {"url": u, "status": r.status_code}
                    if r.status_code == 200:
                        case.update(compare(r.json(), pts[f"{lat},{lon}"]))
                    res["cases"].append(case)
                pr = c.get(f"/api/mhw/plot?bbox=(150,-15,152,-13)&period={months[-1][:4]}&sstvar=sst&mode=month")
                res["plot"] = {"status": pr.status_code, "content_type": pr.headers.get("content-type"),
                               "png": pr.content[:8] == b"\x89PNG\r\n\x1a\n"}
        finally:
            os.chdir(old)
    return with_checks(res, {"points": all(c["status"] == 200 and c.get("match") is True for c in res["cases"]),
                             "plot": res["plot"]["status"] == 200 and res["plot"]["png"]})


def main(argv=None) -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("fixtures"); s.add_argument("--store", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("snapshot"); s.add_argument("--store", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("meta-copy"); s.add_argument("--store", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("diff")
    for k in ("--before", "--after", "--out"):
        s.add_argument(k, required=True)
    s = sub.add_parser("run"); s.add_argument("--store", required=True); s.add_argument("--fixtures", required=True)
    s.add_argument("--api", required=True); s.add_argument("--out")
    s = sub.add_parser("verify")
    for k in ("--store", "--fixtures", "--before", "--api", "--out"):
        s.add_argument(k, required=True)
    s.add_argument("--before-store", required=True,
                   help="pre-update store or meta-copy (required: the metadata rules and the month-count baseline)")
    s.add_argument("--earlier-fixtures", action="append", help="fixtures of months appended before (chained)")
    s.add_argument("--no-api", action="store_true", help="skip the in-process API check (no httpx here)")
    s = sub.add_parser("api-check")
    for k in ("--store", "--fixtures", "--api", "--out"):
        s.add_argument(k, required=True)
    s.add_argument("--earlier-fixtures", action="append")
    s = sub.add_parser("expect-points"); s.add_argument("--fixtures", action="append", required=True)
    s.add_argument("--out", required=True)
    s = sub.add_parser("check-result"); s.add_argument("--file", required=True)
    s.add_argument("--expect", choices=("ok", "fail"), required=True)
    s.add_argument("--expect-failed", action="append", help="with --expect fail: a check that must have failed")
    a = p.parse_args(argv)
    rc = {"fixtures": cmd_fixtures, "snapshot": cmd_snapshot, "meta-copy": cmd_meta_copy, "diff": cmd_diff,
          "run": cmd_run, "verify": cmd_verify, "api-check": cmd_api_check, "expect-points": cmd_expect_points,
          "check-result": cmd_check_result}[a.cmd](a)
    return EXIT_OK if rc is None else rc


if __name__ == "__main__":
    sys.exit(main())
