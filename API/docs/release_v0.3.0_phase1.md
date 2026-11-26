# v0.3.0 (draft): Phase 1 core milestone — API performance and runtime environment

Status: **draft, not tagged.** Date of the summary: 2026-11-26. The changes below have been in production since 2026-11-25; the version number v0.3.0 is planned and not yet released.

The planned v0.3.0 covers the core of Phase 1 of the MHW improvement roadmap:
- the API code fixes that keep outputs unchanged (P1-0 and P1-1);
- the new runtime environment (P1-2).

Both are deployed on the production API host and passed acceptance there. **Phase 1 as a whole is not finished**: the monthly-update hardening (P1-3) is still open, the output-contract changes (P1-4) are deferred, and the zarr rewrite (P1-5) is optional.

## 1. Scope

| Roadmap item | Status |
|---|---|
| P1-0 test net and baseline | Done: golden response tests and baseline measurements (`tests/`, `tools/bench/`) |
| P1-1 API fixes, output unchanged | Done and deployed (2026-11-25) |
| P1-2 runtime environment | Deployed (2026-11-25). Reboot recovery is not verified (postponed). Moving TLS termination to nginx (needs root) is deferred |
| P1-3 monthly update flow | Open. Only the cron wrapper part is in place (§5) |
| P1-4 output-contract changes | Deferred, each item needs an owner decision |
| P1-5 zarr store rewrite | Optional; no decision |
| P2-0 GeoServer measurement tools | Can start (§7) |

## 2. What changed

API (P1-1), all with byte-identical data responses:
- **Data path.** The zarr store is opened without Dask chunks (`chunks=None`), and only the requested variables are read. The API no longer uses the shared Dask scheduler; `src/dask_client_manager.py` is a no-op shim kept for old imports.
- **CSV.** The temporary CSV file is deleted after the response is sent.
- **Plots.** Plots read only the plotted variable and draw on `matplotlib.figure.Figure` objects, not on global pyplot state.
- **Heavy requests.** Plots, and JSON/CSV selections larger than `MHW_HEAVY_BYTES` (default 2 MiB), run in their own thread lane with `MHW_HEAVY_SLOTS` tokens per worker (default 1). They wait in the event loop, not in the default threadpool, so light queries keep running. A lock serialises matplotlib drawing. See `src/guard.py`.

Runtime (P1-2):
- **One uv venv** is shared by the API and the monthly updater (`update_mhw.py`): CPython 3.14.10 with hash-locked dependencies. The installation is described in [`requirements/p1-2/README.md`](../requirements/p1-2/README.md).
- **Launchers.**
  - `conf/p12_env.sh` is the single place that defines the venv path and sets `POLARS_FORCE_PKG=compat`.
  - `conf/start_app.sh` execs gunicorn: two `uvicorn_worker.UvicornWorker` workers, no `--reload`.
  - `conf/run_update_mhw.sh` runs the updater with the same venv.
- **PM2.** `conf/ecosystem.config.js` manages only the `mhwapi` app, with an explicit `cwd` and `kill_timeout`, and no `pre_stop`.
- **polars on the production CPU.** The production CPU lacks AVX2. `polars[rtcompat]` installs the compat runtime, and `POLARS_FORCE_PKG=compat` makes a missing compat runtime fail at import instead of falling back to the AVX2 build.

## 3. Acceptance evidence (summary)

| Check | Result |
|---|---|
| Golden data cases on production (17 JSON/CSV cases) | Status, headers and body sha256 identical to the golden set, after both the P1-1 and the P1-2 switch. The golden set records the 0.2.8 behaviour, so "identical" means *unchanged*, not necessarily *correct* (J08, §4) |
| CSV for the 180-degree case (J08), Range requests, PNG | The J08 CSV equals the output of the P1-1 candidate; it was **not** compared with a 0.2.8 CSV golden (the J08 golden is JSON). Range 206/400 identical, 416 see §4; plot is a 1200×800 PNG |
| Concurrency 1/2/4 (production loopback) | All requests 200 |
| Single-point latency (production loopback, warm, median) | About 3.0 s before Phase 1 → 0.093 s with P1-1 on the old environment → **0.121 s** with P1-2 |
| P1-2 vs old environment (paired, 30 rounds, interleaved) | New environment **about +32 ms slower** (95% CI upper bound about +33 ms). The old environment was measured as an isolated single-worker plain-HTTP candidate, so TLS and worker count differ from production. Within the **provisional** limits (median ≤ 0.3 s; CI upper bound ≤ +50 ms); the owner has **not** formally adopted these limits |
| Production DB read check (read-only, one session) | psycopg 3 connects; column types and NULL handling of the last month match the updater's expectations; that month's DB data equals the store |
| Rollback | Rehearsed three times in an isolated PM2 instance before the switch; not needed |
| Full test suite (last run at `a2a2f64`) | 459 passed / 7 skipped in the CPython 3.14.10 test venv; 458 passed / 8 skipped in the P1-1 environment |

## 4. Accepted differences and known issues

- **CSV Range 416.** The response's `Content-Range` is now `bytes */<size>`, the Starlette 1.x format. The owner accepted this.
- **PNG plots** are verified by format and size, not by bytes (owner decision for the new library versions).
- **Dateline: known bug, not fixed.**
  - A single bounding-box query crossing 180 degrees (golden case J08: `lon0=170&lon1=-170`) currently returns the **wrong range**: the 340-degree complement instead of the 20-degree box.
  - Because that area is larger, the requested period is also cut (24 months requested, 12 returned).
  - Workshop users work around this with two queries, one on each side of 180 degrees.
  - Phase 1 deliberately kept this behaviour unchanged: the J08 JSON golden records this wrong output, and the J08 CSV was only shown to equal the P1-1 candidate's output.
- **`lon0=0`.** A query with `lon0=0` is treated as crossing 0 degrees and reads the data twice: the output is correct but costs more. This is a P1-4 candidate.
- **Shared Dask services.** They stay running for other users; the MHW API and updater no longer need them. Whether to keep them is not an MHW decision.

## 5. Monthly update: what is verified and what is not

**Fixture tests on a store copy** (isolated test host, two consecutive appends with the production lock):
- The store's date chunk size is 2.
- The first append (538 → 539 months) added 64 data chunks: 4 variables × 16 spatial chunks.
- The second append (539 → 540) rewrote the same 64 chunks and added none. So every second month rewrites the chunks that hold the previous, already published month.
- Files added or rewritten totalled about 8.5 MB and 14.6 MB (sum of file sizes, not measured I/O).
- The explicit `consolidate_metadata` call took about 0.007 s each time.
- In both runs: the previous month's values were unchanged, the new month equalled its fixture, and the API (new environment, and the zarr 2 rollback environment) read both months.

These are two isolated fixture runs, **not** a real NOAA/DB monthly update and **not** a long-term performance result. They neither show that a rechunk is needed nor that it is not.

**Production cron wrapper.** The production cron job calls a host-side wrapper; revision r2 was installed on 2026-11-26. It:
- runs `conf/run_update_mhw.sh` (the new venv);
- restarts `mhwapi` only after a successful append, through PM2's existing daemon;
- waits for HTTP 200;
- exits non-zero on every failure;
- when some NOAA daily files fail to download, still restarts but prints a WARNING. When all of them fail, the updater aborts and the wrapper does not restart.

The wrapper is tested with fakes and one controlled restart, but **the first real update on the new venv is 2026-12-18**. Its log, cron mail and the published month still have to be checked.

**Version-control gap.** The wrapper itself is not in this repository. `API/update_mhw.sh` is a symlink to its host path, so a clone only gets a dangling link (§6).

**Still open for P1-3:**
- a daily retry window and an idempotent run;
- day-count completeness checks, `.part` downloads, request timeouts and retries;
- verify-then-publish ordering;
- backup and restore of the 64 chunks a second-month append rewrites;
- API read consistency during an append;
- an interrupted append;
- clean-up of downloaded NOAA files.

## 6. Follow-up work and owner decisions

1. Check the 2026-12-18 monthly update (log `# result:` line, cron mail, new month visible in the API).
2. Formally adopt or change the latency limits (0.3 s; +50 ms).
3. Production checkout: it still points at the pre-Phase-1 commit, with the deployed files copied in. Decide on push or tag so that files and HEAD agree. Since this documentation commit, the repo also differs from production in comments of `conf/p12_env.sh`, `conf/run_update_mhw.sh` and `conf/ecosystem.config.js` (comments only).
4. Cron wrapper in version control. Proposal: replace the `API/update_mhw.sh` symlink with a parameterised template (paths and health URL as variables at the top), copied from the verified r2 wrapper, which contains host paths but no secrets. This needs owner approval.
5. Reboot recovery check (PM2 resurrect with the new app), postponed.
6. When to retire the pre-Phase-1 Python environment, which the rollback script still uses.
7. `conf/simu.sh` must not be run as it is (it pip-installs into the current Python, kills processes by name and calls `pm2 start`); see [`README.md`](../README.md#development-and-deployment-maintainers). Moving it to the uv locks would be a small separate change.

## 7. Next: P2-0

P2-0 is read-only measurement and can start: repeatable GeoServer timing at three levels (public URL, GeoServer loopback, database), cache hit state, and PNG byte comparison. GeoServer layer performance is the main remaining target. Any change to GeoServer, the database or their schedules needs separate authorization.
