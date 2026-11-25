# P1-2 dependency set (official for the P1-2 environment; doc 29 r5 §3)

One uv venv for the API **and** the monthly updater (`update_mhw.py`), owner
decision 2026-11-25. The older `../../requirements.txt` and `../../Pipfile` are
kept unchanged for the pre-P1-2 environment; they are **not** the P1-2 lock.

| file | role | sha256 |
|---|---|---|
| `runtime.in` | production input: API set + updater set (`sqlalchemy`, `dask[array]`, `netcdf4`, `psycopg[binary]`); polars as `polars[rtcompat]` | `623aaab8…41802` |
| `runtime.lock` | **production lock** (hashes, universal), 56 packages | `73d273ce…6d952d27` |
| `dev.in` / `dev.lock` | tests only (`-c runtime.lock`): pytest, psutil, httpx, dask[complete], distributed | `7b8dcae5…` / `451f225e…cfd63a77` |
| `tools.in` / `tools.lock` | benchmark runner only (`tools/bench`, `tools/p1_2`): psutil, httpx | `677eb367…` / `60d5704f…bf446e` |

Provenance (no re-resolution in this commit):
- **`runtime.in` / `runtime.lock`** are byte-identical copies of V-c's `prod-B.in` / `prod-B.lock`, which is why the lock's header still names `prod-B.*`. That lock was installed on vm24 and validated in V-c (doc 28 r4 §9).
- **`dev.in` / `dev.lock`** are byte-identical to V-a, which ran full pytest on the Mac and vm214. All 20 packages it shares with `runtime.lock` have the same pinned versions.
- **Resolver:** uv 0.13.10, `uv pip compile --universal --python-version 3.14 --generate-hashes --exclude-newer 2026-11-25T00:00:00Z`.
- **Python:** CPython 3.14.10, uv-managed.

Install (into a new uv venv only; never into a shared Python):

```bash
uv venv --python 3.14.10 "$MHW_VENV"
uv pip sync --python "$MHW_VENV/bin/python" --require-hashes --no-build runtime.lock
```

Launchers: `../../conf/p12_env.sh` (venv path, `POLARS_FORCE_PKG=compat`), `../../conf/start_app.sh` (API), `../../conf/run_update_mhw.sh` (updater).
Polars: the extra is `rtcompat`; the runtime it loads is `compat`, forced by `POLARS_FORCE_PKG=compat` so a missing compat runtime fails at import instead of falling back to runtime-32.
