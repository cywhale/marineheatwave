"""Isolated check of update_mhw.py against a writable COPY of the MHW store (P1-2 V-b).

No database, no NOAA, no network: the PostgreSQL query result and the downloaded
NOAA OISST daily files are fixtures built from the store's own last month.

    updater_check.py fixtures --store COPY --out FIXDIR        # db.parquet + oisst-*.nc + fixture.json
    updater_check.py snapshot --store COPY --out before.json   # per-file sha256 of the store
    updater_check.py run      --store COPY --fixtures FIXDIR --api API_DIR
    updater_check.py verify   --store COPY --fixtures FIXDIR --before before.json --before-store RO_COPY \
                              --api API_DIR --out report.json
    updater_check.py api-check --store COPY --fixtures FIXDIR --api API_DIR --out api.json  # any API env

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
SAMPLES = [(-14.875, 150.125), (0.125, 200.125), (35.125, 140.125), (-60.125, 300.125), (10.125, 125.125)]


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
    root = Path(store)
    return {str(p.relative_to(root)): _sha(p) for p in sorted(root.rglob("*")) if p.is_file()}


def cmd_snapshot(a):
    snap = snapshot(a.store)
    Path(a.out).write_text(json.dumps(snap))
    print(f"{len(snap)} files")


# ----------------------------------------------------------------------- run


class _NoNetwork(OSError):
    pass


def _block_network():
    def refuse(*a, **k):
        raise _NoNetwork("network disabled in updater_check")
    socket.socket.connect = refuse
    socket.socket.connect_ex = refuse
    socket.create_connection = refuse
    socket.getaddrinfo = refuse


def cmd_run(a):
    import pandas as pd

    fix = Path(a.fixtures).resolve()
    store = str(Path(a.store).resolve())
    _block_network()

    import dotenv
    dotenv.load_dotenv = lambda *x, **k: False  # never read a .env
    for k in ("DBUSER", "DBPASS", "DBHOST", "DBPORT", "DBNAME", "MHWTABLE"):
        os.environ[k] = "fixture"

    queries = []
    db = pd.read_parquet(fix / "db.parquet")

    class Engine:
        def dispose(self):
            queries.append("dispose")

    def fake_read_sql_query(q, engine, *x, **k):
        queries.append(q)
        assert isinstance(engine, Engine)
        return db.copy()

    import requests
    import sqlalchemy

    def no_requests(*x, **k):
        raise _NoNetwork("requests.get called: a fixture file is missing")

    requests.get = no_requests
    sqlalchemy.create_engine = lambda *x, **k: Engine()
    pd.read_sql_query = fake_read_sql_query

    sys.path.insert(0, str(Path(a.api).resolve()))
    import update_mhw as upd
    upd.create_engine = sqlalchemy.create_engine  # bound at import ("from sqlalchemy import create_engine")
    upd.ZARR_PATH = store
    upd.DATA_PATH = str(fix) + "/"

    with tempfile.TemporaryDirectory() as cwd:  # relative paths in the module must not matter
        os.chdir(cwd)
        t = time.perf_counter()
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            upd.append_to_zarr()
        dt = time.perf_counter() - t
    out = buf.getvalue()
    print(out)
    print(json.dumps({"seconds": round(dt, 1), "queries": queries,
                      "completed": "All work append to zarr done" in out}, indent=1))


# -------------------------------------------------------------------- verify


def cmd_verify(a):
    import numpy as np
    import pandas as pd
    import zarr

    fix = Path(a.fixtures)
    meta = json.loads((fix / "fixture.json").read_text())
    before = json.loads(Path(a.before).read_text())
    after = snapshot(a.store)
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
    meta_diff = _metadata_diff(a.before_store, a.store, meta_files) if a.before_store else None
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
    unexplained = data_added + data_changed
    if a.before_store:
        unexplained += metadata_problems(a.before_store, a.store, n_old)
    else:
        unexplained += [f"{f}: no pre-update store to compare metadata with" for f in meta_files]
    rep["files"]["unexplained"] = unexplained

    # structure
    rep["dates"] = {"n": len(dates), "last": str(dates[-1].date()), "prev": str(dates[-2].date()),
                    "monotonic_unique": bool(dates.is_monotonic_increasing and dates.is_unique)}
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

    # API on the updated copy (its own startup code, relative data/mhw.zarr)
    rep["api"] = _api_check(a.api, a.store, meta, exp)

    f = rep["files"]
    v = rep["values"]
    ok = (not f["removed"] and not f["unexplained"]
          and rep["dates"]["n"] == n_old + 1 and rep["dates"]["last"] == meta["new_month"]
          and rep["dates"]["prev"] == meta["last_month_in_store"] and rep["dates"]["monotonic_unique"]
          and all(x["nan_pattern_equal"] and x["within_tol"] and x["prev_month_unchanged"] for x in v.values())
          and rep["dtypes"] == {"sst": "float32", "sst_anomaly": "float64", "level": "int64", "td": "float64"}
          and all(s[0] == n_old + 1 for s in rep["zarr_shapes"].values())
          and all(s[0] == n_old + 1 for s in rep["consolidated_shapes"].values())
          and rep["api"]["ok"])
    rep["ok"] = bool(ok)
    Path(a.out).write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps({"ok": rep["ok"], "files": {k: (x if not isinstance(x, list) or len(x) < 20 else f"{len(x)} items")
                                                  for k, x in f.items()},
                      "dates": rep["dates"], "values": v, "api": rep["api"]}, indent=1, default=str))


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
KNOWN_REWRITES = {
    ("zmetadata", f"{GROUP}/.zgroup", "consolidated_metadata"): lambda b, a: b is None and a is not None,
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


def cmd_api_check(a):
    """The API on the updated copy only; runs in any API environment (e.g. the
    zarr 2 / P1-1 one, to check that a store written by zarr 3 stays readable)."""
    import zarr
    fix = Path(a.fixtures)
    meta = json.loads((fix / "fixture.json").read_text())
    res = _api_check(a.api, a.store, meta, _expected(fix))
    res["zarr"] = zarr.__version__
    Path(a.out).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({"ok": res["ok"], "zarr": res["zarr"], "plot": res["plot"],
                      "cases": [(c["status"], c.get("n_new_month"), c.get("match")) for c in res["cases"]]}))


def _api_check(api_dir, store, meta, exp):
    import numpy as np

    api = Path(api_dir).resolve()
    res = {"cases": []}
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
                    i, j = int(round((lat + 89.875) / 0.25)), int(round((lon - 0.125) / 0.25))
                    qlon = lon if lon <= 180 else lon - 360
                    u = (f"/api/mhw?lon0={qlon}&lat0={lat}&start={meta['last_month_in_store']}"
                         f"&end={meta['new_month']}&append=sst,sst_anomaly,level,td")
                    r = c.get(u)
                    case = {"url": u, "status": r.status_code}
                    if r.status_code == 200:
                        rows = r.json()
                        last = [x for x in rows if str(x.get("date", ""))[:10] == meta["new_month"]]
                        case["n_rows"], case["n_new_month"] = len(rows), len(last)
                        if last:
                            x = last[0]
                            want = {"sst": exp["sst"][i, j], "sst_anomaly": exp["sst_anomaly"][i, j],
                                    "td": exp["td"][i, j], "level": exp["level"][i, j]}
                            got = {k: x.get(k) for k in want}
                            case["got"], case["want"] = got, {k: (None if np.isnan(w) else float(w)) for k, w in want.items()}
                            case["match"] = all(
                                (got[k] is None and case["want"][k] is None)
                                or (got[k] is not None and case["want"][k] is not None
                                    and abs(got[k] - case["want"][k]) <= (0.006 if k == "sst" else 1e-9))
                                for k in want)
                    res["cases"].append(case)
                pr = c.get(f"/api/mhw/plot?bbox=(150,-15,152,-13)&period={meta['new_month'][:4]}&sstvar=sst&mode=month")
                res["plot"] = {"status": pr.status_code, "content_type": pr.headers.get("content-type"),
                               "png": pr.content[:8] == b"\x89PNG\r\n\x1a\n"}
        finally:
            os.chdir(old)
    res["ok"] = (all(c["status"] == 200 and c.get("n_new_month") == 1 and c.get("match") for c in res["cases"])
                 and res["plot"]["status"] == 200 and res["plot"]["png"])
    return res


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("fixtures"); s.add_argument("--store", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("snapshot"); s.add_argument("--store", required=True); s.add_argument("--out", required=True)
    s = sub.add_parser("run"); s.add_argument("--store", required=True); s.add_argument("--fixtures", required=True)
    s.add_argument("--api", required=True)
    s = sub.add_parser("verify")
    for k in ("--store", "--fixtures", "--before", "--api", "--out"):
        s.add_argument(k, required=True)
    s.add_argument("--before-store", help="read-only pre-update store, to classify metadata rewrites")
    s = sub.add_parser("api-check")
    for k in ("--store", "--fixtures", "--api", "--out"):
        s.add_argument(k, required=True)
    a = p.parse_args()
    {"fixtures": cmd_fixtures, "snapshot": cmd_snapshot, "run": cmd_run, "verify": cmd_verify,
     "api-check": cmd_api_check}[a.cmd](a)


if __name__ == "__main__":
    main()
