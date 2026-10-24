#!/bin/bash

DASK_TLS_DIR=/home/odbadmin/.config/dask/tls
export DASK_SCHEDULER_ADDRESS=tls://127.0.0.1:8786
export DASK_TLS_CA_FILE="$DASK_TLS_DIR/certs/ca.crt"
export DASK_TLS_CLIENT_CERT="$DASK_TLS_DIR/certs/mhw-client.crt"
export DASK_TLS_CLIENT_KEY="$DASK_TLS_DIR/private/mhw-client.key"
for tls_file in "$DASK_TLS_CA_FILE" "$DASK_TLS_CLIENT_CERT" "$DASK_TLS_CLIENT_KEY"; do
    if [[ ! -r "$tls_file" ]]; then
        echo "Required Dask TLS file is not readable: $tls_file" >&2
        exit 1
    fi
done

gunicorn mhw_app:app -w 2 -k uvicorn.workers.UvicornWorker -b 127.0.0.1:8030 --keyfile conf/privkey.pem --certfile conf/fullchain.pem --timeout 180 --reload
