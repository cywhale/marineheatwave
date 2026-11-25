#!/bin/bash
# Monthly updater with the SAME P1-2 uv venv as the API (doc 29 r5 D1).  For the
# vm24 cron wrapper to call in V-d (the wrapper itself lives outside the repo and
# is not changed here).  cwd = API dir: update_mhw.py uses data/mhw.zarr and its
# .env relative to it.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source conf/p12_env.sh
exec "$MHW_VENV/bin/python" update_mhw.py "$@"
