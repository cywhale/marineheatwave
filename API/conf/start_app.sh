#!/bin/bash

# P1-1a: the API no longer connects to the shared Dask scheduler, so the Dask
# TLS check and exports are gone.  Clear any inherited values so no library
# picks up a scheduler address or client credentials from the environment.
unset DASK_SCHEDULER_ADDRESS DASK_TLS_CA_FILE DASK_TLS_CLIENT_CERT DASK_TLS_CLIENT_KEY

gunicorn mhw_app:app -w 2 -k uvicorn.workers.UvicornWorker -b 127.0.0.1:8030 --keyfile conf/privkey.pem --certfile conf/fullchain.pem --timeout 180 --reload
