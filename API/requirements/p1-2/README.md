# MHW API dependency locks (current)

These files are the **official dependency set** of the MHW API. Production has used them since 2026-11-25; they ship with v0.3.0.
- One uv venv serves both the API (`mhw_app.py`) and the monthly updater (`update_mhw.py`).
- The older `../../requirements.txt` and `../../Pipfile` belong to the environment used before 2026-11-25. They are kept unchanged for reference and the legacy `conf/simu.sh`; they are **not** current.

## Files

| file | used for | sha256 (`SHA256SUMS`) |
|---|---|---|
| `runtime.in` → `runtime.lock` | **Production**: the API plus the updater's packages (`sqlalchemy`, `psycopg[binary]`, `dask[array]`, `netcdf4`). 56 packages | `623aaab8…` → `73d273ce…` |
| `dev.in` → `dev.lock` | **Tests only**, installed together with `runtime.lock`: pytest, psutil, httpx, `dask[complete]`, `distributed`. Constrained by the runtime pins (`-c runtime.lock`); all packages the two locks share have the same versions | `7b8dcae5…` → `451f225e…` |
| `tools.in` → `tools.lock` | **Benchmark and check runners only** (`tools/bench`, `tools/p1_2`): psutil, httpx. A separate small venv, never the API's | `677eb367…` → `60d5704f…` |

**Why the runtime lock contains Dask.** The updater opens the store and the NOAA daily files through xarray with Dask chunks, so `dask[array]` is a runtime dependency. The **API does not use Dask** for queries: it reads the store with `chunks=None` and does not connect to the shared Dask scheduler. `distributed` is **not** in the runtime lock; it is only in `dev.lock` for tests. The launchers also clear the old scheduler settings (`DASK_SCHEDULER_ADDRESS` and the three `DASK_TLS_*` variables).

**polars.** The extra is `polars[rtcompat]`, which installs the `compat` runtime. The launchers set `POLARS_FORCE_PKG=compat`: the production CPU lacks AVX2, and this setting makes a missing compat runtime fail at import instead of falling back to a build that needs AVX2.

## Python and uv

- **Python:** CPython 3.14.10, uv-managed (`uv python install 3.14.10`).
- **Locks:** generated with uv 0.13.10, as universal locks with hashes, and resolved against the package index as of 2026-11-25:

  ```bash
  uv pip compile --universal --python-version 3.14 --generate-hashes --exclude-newer 2026-11-25T00:00:00Z runtime.in -o runtime.lock
  ```

  The header of `runtime.lock` names `prod-B.in` / `prod-B.lock`. The file was validated under that name before it was copied here unchanged.
- **Changing a lock.** Edit the `.in` file, recompile with the command above (new `--exclude-newer` date), and update `SHA256SUMS`. Then rerun the full test suite and the production acceptance checks before deploying. `tests/test_p1_2_launch_config.py` checks the sums and that the shared pins agree.

## Install

Always into a **new** venv; never into a shared or system Python.

Verify the files first:

```bash
cd API/requirements/p1-2
sha256sum -c SHA256SUMS        # macOS: shasum -a 256 -c SHA256SUMS
```

Production (API and updater); `MHW_VENV` must match `conf/p12_env.sh`:

```bash
uv venv --python 3.14.10 "$MHW_VENV"
uv pip sync --python "$MHW_VENV/bin/python" --require-hashes --no-build runtime.lock
```

Tests (from the `API/` directory: `pytest`, see `pytest.ini`):

```bash
uv venv --python 3.14.10 .venv-test
uv pip sync --python .venv-test/bin/python --require-hashes runtime.lock dev.lock
```

The Phase 1 test venv was built this way with an earlier runtime lock, which lacked only the updater's packages. This combination with the current `runtime.lock` is expected to resolve, because the shared pins agree, but it has not been installed as such.

Benchmark/check runners:

```bash
uv venv --python 3.14.10 .venv-tools
uv pip sync --python .venv-tools/bin/python --require-hashes --no-build tools.lock
```

## How the venv is used

- `../../conf/p12_env.sh`: the venv path (`MHW_VENV`, default `/home/odbadmin/venvs/mhw-p12`; an `MHW_VENV` set in the environment wins) and `POLARS_FORCE_PKG=compat`. It refuses to continue if the venv has no python.
- `../../conf/start_app.sh`: the API (gunicorn + uvicorn-worker), started by PM2 through `../../conf/ecosystem.config.js`.
- `../../conf/run_update_mhw.sh`: the monthly updater, called by the host's cron wrapper. That wrapper lives on the host; its template is `../../conf/update_mhw_cron.sh.template`.
