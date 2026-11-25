"""P1-1a targeted checks (spec 17 P1-1a): direct zarr reads without Dask.

The golden suite already proves byte-identical responses.  These tests cover
what golden cannot see:

* the production source no longer references Dask clients or ``.compute()``
  (a static check, so the harness patches cannot hide a leftover call);
* an *unpatched* app process (no tests.harness) serves requests without
  constructing a ``distributed.Client`` or calling the shim;
* raw queries read only the requested variables' chunks from disk;
* the new open path returns the same frames (values, dtypes, row order,
  cross-longitude concat order) as the 0.2.8 ``open_zarr(chunks='auto')`` path.
"""
import ast
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

API = Path(__file__).resolve().parents[1]
RUNTIME_SOURCES = ["mhw_app.py", "src/mhw_utils.py", "src/mhw_plot.py", "src/config.py"]
COORDS = {"lon", "lat", "date"}
DATA_VARS = {"level", "sst", "sst_anomaly", "td"}


# --------------------------------------------------------------- static source


def _imports_and_calls(path: Path):
    tree = ast.parse(path.read_text())
    imports, names, attr_calls = set(), set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add(node.module or "")
            names.update(a.name for a in node.names)
        elif isinstance(node, ast.Name):
            names.add(node.id)
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            attr_calls.add(node.func.attr)
    return imports, names, attr_calls


@pytest.mark.parametrize("rel", RUNTIME_SOURCES)
def test_runtime_source_has_no_dask_client_or_compute(rel):
    imports, names, attr_calls = _imports_and_calls(API / rel)
    assert not {m for m in imports if m.split(".")[0] in ("dask", "distributed")}, imports
    assert "src.dask_client_manager" not in imports
    assert not names & {"get_dask_client", "close_dask_client", "Client"}
    assert "compute" not in attr_calls


def test_shim_is_noop_and_does_not_import_distributed():
    imports, _, _ = _imports_and_calls(API / "src/dask_client_manager.py")
    assert not imports, imports
    from src import dask_client_manager as shim

    assert shim.get_dask_client("mhwapi") is None
    assert shim.close_dask_client("mhwapi") is None


def test_launcher_unsets_dask_env_and_no_longer_exports_it():
    sh = (API / "conf/start_app.sh").read_text()
    assert "unset DASK_SCHEDULER_ADDRESS DASK_TLS_CA_FILE DASK_TLS_CLIENT_CERT DASK_TLS_CLIENT_KEY" in sh
    assert "export DASK" not in sh and "tls_file" not in sh
    # unchanged parts (P1-2 territory)
    assert "-w 2 -k uvicorn.workers.UvicornWorker -b 127.0.0.1:8030" in sh and "--reload" in sh


# ----------------------------------------------------- unpatched app process

_UNPATCHED = r"""
import json, os, resource, sys
os.chdir(sys.argv[1])
sys.path.insert(0, sys.argv[1])
calls = {"client_init": 0, "shim_get": 0, "shim_close": 0}

# Count, do not hide: any Client construction or shim call is recorded.  The
# app modules are imported afterwards, so a leftover
# "from src.dask_client_manager import get_dask_client" would bind to these.
import distributed
_orig_init = distributed.Client.__init__
def _init(self, *a, **k):
    calls["client_init"] += 1
    raise RuntimeError("Client constructed")
distributed.Client.__init__ = _init
import src.dask_client_manager as shim
def _get(*a, **k):
    calls["shim_get"] += 1
def _close(*a, **k):
    calls["shim_close"] += 1
shim.get_dask_client, shim.close_dask_client = _get, _close

import mhw_app, src.config as config
from fastapi.testclient import TestClient
from dask.base import is_dask_collection

urls = json.loads(sys.argv[2])
statuses = {}
from tests.zarr_spy import chunk_reads  # zarr 2 DirectoryStore or zarr 3 LocalStore

with TestClient(mhw_app.app) as c:
    # Type checks only: never touch .data / .values, which would materialise a
    # whole global variable (sst 2.2 GB, f8 vars 4.5 GB) on this lazy path.
    with chunk_reads() as inspect_reads:
        lazy = {k: type(config.dz[k].variable._data).__name__ for k in config.dz.data_vars}
        dask_backed = [k for k in config.dz.data_vars if is_dask_collection(config.dz[k].variable._data)
                       or config.dz[k].chunks is not None]
    for u in urls:
        statuses[u] = c.get(u).status_code
print(json.dumps({"calls": calls, "statuses": statuses, "dask_backed": dask_backed, "lazy": lazy,
                  "inspect_reads": inspect_reads,
                  "maxrss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                  / (1024 * 1024 if sys.platform == "darwin" else 1024),
                  "default_client": distributed.client._get_global_client() is not None}))
"""

UNPATCHED_URLS = [
    "/api/mhw?lon0=135&lat0=15&start=2020-01-01&end=2020-06-01",
    "/api/mhw?lon0=-5&lon1=5&lat0=0&lat1=5&start=2020-01-01&end=2020-02-01&append=sst,td",
    "/api/mhw/csv?lon0=170&lon1=-170&lat0=-5&lat1=5&start=2020-01-01&end=2020-02-01&append=sst",
    "/api/mhw/plot?bbox=(135,15,137,17)&period=2019-202006&sstvar=sst&mode=series",
]


def test_unpatched_app_never_creates_dask_client(tmp_path):
    """Run the app without tests.harness, with a bogus scheduler address in the
    environment, and count (not refuse) every Client / shim use."""
    import json

    env = {k: v for k, v in os.environ.items() if not k.startswith("DASK_")}
    env.update(MPLBACKEND="Agg", DASK_SCHEDULER_ADDRESS="tcp://127.0.0.1:9")
    out = subprocess.run([sys.executable, "-c", _UNPATCHED, str(API), json.dumps(UNPATCHED_URLS)],
                         env=env, capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, out.stderr[-3000:]
    res = json.loads(out.stdout.strip().splitlines()[-1])
    assert res["calls"] == {"client_init": 0, "shim_get": 0, "shim_close": 0}
    assert res["default_client"] is False
    assert res["dask_backed"] == []
    assert res["inspect_reads"] == []  # the type checks read nothing from disk
    # whole-variable materialisation would be >= 2.2 GB (sst); the four small requests need far less
    assert res["maxrss_mb"] < 1500, res["maxrss_mb"]
    assert set(res["lazy"].values()) == {"CopyOnWriteArray"}
    assert all(s == 200 for s in res["statuses"].values()), res["statuses"]
    assert "Dask scheduler unavailable" not in out.stderr


def test_opened_dataset_has_no_dask_or_memory_cache(client):
    import src.config as config

    with chunk_reads() as keys:  # inspection must not materialise any variable
        for k in config.dz.data_vars:
            x, chain = config.dz[k].variable._data, []
            while x is not None:
                chain.append(type(x).__name__)
                x = getattr(x, "array", None)
            assert chain == ["CopyOnWriteArray", "LazilyIndexedArray", "ZarrArrayWrapper"], (k, chain)
            assert config.dz[k].chunks is None
    assert keys == []


# ------------------------------------------------------------ disk-level reads
# The spy lives in tests/zarr_spy.py (zarr 2 DirectoryStore and zarr 3 LocalStore);
# metadata documents (.zarray, .zmetadata, zarr.json, ...) never count as a variable read.

from tests.zarr_spy import ZARR_MAJOR, chunk_reads, data_keys  # noqa: E402
from tests.zarr_spy import vars_read as _vars_read  # noqa: E402


def test_spy_records_a_small_slice_and_nothing_for_lazy_inspection(client):
    """Self-check of the spy itself, so a silent spy cannot make the read
    assertions below pass vacuously ."""
    import src.config as config

    with chunk_reads() as keys:
        kinds = {k: type(config.dz[k].variable._data).__name__ for k in config.dz.data_vars}
        _ = config.dz["sst"].chunks, config.dz["sst"].dtype, config.dz["sst"].shape
    assert kinds and data_keys(keys) == [], data_keys(keys)[:5]  # lazy inspection reads no chunk

    with chunk_reads() as keys:
        v = config.dz["sst"].isel(date=0, lat=slice(400, 402), lon=slice(500, 502)).values
    assert v.shape == (2, 2)
    assert _vars_read(keys) - COORDS == {"sst"}, sorted(keys)[:10]
    assert all(k.startswith("anomaly/sst/") for k in data_keys(keys) if "/sst/" in k)


@pytest.mark.parametrize("url,expected", [
    ("/api/mhw?lon0=135&lat0=15&start=2020-01-01&end=2020-06-01", {"level"}),
    ("/api/mhw?lon0=135&lat0=15&start=2020-01-01&end=2020-06-01&append=sst", {"sst"}),
    ("/api/mhw?lon0=120&lon1=125&lat0=10&lat1=15&start=2020-01-01&end=2020-03-01&append=td,sst_anomaly",
     {"td", "sst_anomaly"}),
    ("/api/mhw?lon0=-5&lon1=5&lat0=0&lat1=5&start=2020-01-01&end=2020-02-01&append=sst", {"sst"}),
    ("/api/mhw/csv?lon0=170&lon1=-170&lat0=-5&lat1=5&start=2020-01-01&end=2020-02-01&append=level,sst",
     {"level", "sst"}),
])
def test_raw_reads_only_requested_variables(client, url, expected):
    with chunk_reads() as keys:
        r = client.get(url)
    assert r.status_code == 200
    read = _vars_read(keys) - COORDS
    assert read == expected, (read, sorted(keys)[:20])


def test_plot_reads_only_plotted_variable(client):
    """Was a P1-1a scope marker (plots read all 4 variables); since P1-1d a plot
    reads only its sstvar.  More plot cases are in tests/test_p1_1d.py."""
    with chunk_reads() as keys:
        r = client.get("/api/mhw/plot?bbox=(135,15,137,17)&period=2019-202006&sstvar=sst&mode=series")
    assert r.status_code == 200
    assert _vars_read(keys) - COORDS == {"sst"}


# ------------------------------------------------- equivalence with 0.2.8 path

EQUIV = [
    # (lon0, lat0, lon1, lat1, start, end, append, mode)
    (135, 15, None, None, "1982-01-01", "2026-10-01", None, None),                # point, full series
    (-179.9, 0, None, None, "2020-01-01", "2020-12-01", "sst,td", None),          # point at -180
    (120, 10, 125, 15, "2020-01-01", "2020-12-01", "level,sst,sst_anomaly,td", None),
    (-5, 0, 5, 5, "2020-01-01", "2020-03-01", "sst,level", None),                 # crosses 0
    (5, 0, -5, 5, "2020-01-01", "2020-03-01", "sst_anomaly", None),               # crosses 0, reversed
    (170, -5, -170, 5, "2020-01-01", "2021-12-01", "sst", None),                  # antimeridian (J08 behaviour)
    (-170, -5, 170, 5, "2020-01-01", "2020-02-01", "level", None),
    (150, 0, 150.1, 0.1, "2020-01-01", "2020-02-01", "td", None),                 # sub-grid bbox
    (100, 0, 100.1, 0.1, "2020-01-01", "2020-02-01", "td", None),                 # land: same 400 both ways
    (135, 15, 137, 17, "2030-01-01", "2030-02-01", "sst", None),                  # beyond data: same 400
    (135, 15, 137, 17, "2012-01-01", "2021-12-01", "sst", "area_mean_sst"),       # plot path (full ds)
    (135, 15, 137, 17, "1982-01-01", "2011-12-01", "sst_anomaly", "month_mean"),
]


@pytest.fixture(scope="module")
def old_dz():
    import xarray as xr

    ds = xr.open_zarr("data/mhw.zarr", chunks="auto", group="anomaly", decode_times=True, consolidated=True)
    yield ds
    ds.close()


@pytest.mark.parametrize("args", EQUIV, ids=[f"{a[0]}_{a[1]}_{a[2]}_{a[3]}_{a[6]}_{a[7]}" for a in EQUIV])
def test_same_frames_as_0_2_8_open_path(client, old_dz, args):
    from datetime import date

    from polars.testing import assert_frame_equal

    import src.config as config
    from src.mhw_utils import process_mhw_data

    lon0, lat0, lon1, lat1, start, end, append, mode = args
    kw = dict(lon0=lon0, lat0=lat0, lon1=lon1, lat1=lat1, start=date.fromisoformat(start),
              end=date.fromisoformat(end), append=append, mode=mode)
    from fastapi import HTTPException

    def run():
        try:
            return process_mhw_data(**kw)  # synchronous since P1-1c
        except HTTPException as e:
            return ("HTTPException", e.status_code, e.detail)

    new_res = run()
    saved = config.dz
    config.dz = old_dz
    try:
        old_res = run()
    finally:
        config.dz = saved
    if new_res[0] == "HTTPException" or old_res[0] == "HTTPException":
        assert new_res == old_res  # same error path, status and message
        return
    (new_name, new), (old_name, old) = new_res, old_res
    assert new_name == old_name
    assert new.schema == old.schema
    assert_frame_equal(new, old, check_exact=True)  # values, dtypes and row order
    if "lon" in new.columns:
        assert new["lon"].to_list() == old["lon"].to_list()  # concat order across 0 / 180
