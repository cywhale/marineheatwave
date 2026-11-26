#!/bin/bash
# MHW API launcher, P1-2 (doc 29 r5 D2): the shared P1-2 uv venv, exec (PM2's
# process becomes the gunicorn master), uvicorn-worker, no auto-reload.
# P1-1a: the API no longer connects to the shared Dask scheduler; p12_env.sh
# unsets DASK_SCHEDULER_ADDRESS DASK_TLS_CA_FILE DASK_TLS_CLIENT_CERT DASK_TLS_CLIENT_KEY.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source conf/p12_env.sh
exec "$MHW_VENV/bin/gunicorn" mhw_app:app -w 2 -k uvicorn_worker.UvicornWorker -b 127.0.0.1:8030 --keyfile conf/privkey.pem --certfile conf/fullchain.pem --timeout 180
