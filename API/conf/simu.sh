#!/usr/bin/env bash
# Local development only. Do not run on shared/production hosts.

# Refuse to start unless the scheduler port can be checked reliably.
if ! command -v ss >/dev/null 2>&1; then
    echo 'Refusing to start local Dask: ss is required to check port 8786.' >&2
    exit 1
fi
if ! scheduler_listeners=$(ss -H -ltn 'sport = :8786'); then
    echo 'Refusing to start local Dask: failed to inspect port 8786.' >&2
    exit 1
fi
if [[ -n "$scheduler_listeners" ]]; then
    echo 'Refusing to start local Dask: port 8786 is already in use.' >&2
    exit 1
fi

# install packages
pip3 install -r ./requirements.txt

# Run local-development Dask only on loopback. The worker dashboard is disabled.
dask scheduler --host 127.0.0.1 --port 8786 --dashboard-address 127.0.0.1:0 &
dask worker tcp://127.0.0.1:8786 --host 127.0.0.1 --memory-limit 8GB --no-dashboard --dashboard-address 127.0.0.1:0 &

# run API server
## localhost: gunicorn mhw_app:app -k uvicorn.workers.UvicornWorker -b 127.0.0.1:8030 --timeout 120
gunicorn mhw_app:app -w 2 -k uvicorn.workers.UvicornWorker -b 127.0.0.1:8030 --keyfile conf/privkey.pem --certfile conf/fullchain.pem --reload --timeout 180

# kill process
ps -ef | grep 'mhw_app' | grep -v grep | awk '{print $2}' | xargs -r kill -9

# kill dask
ps -ef | grep -w 'dask scheduler' | grep -v grep | awk '{print $2}' | xargs -r kill -9 && ps -ef | grep -w 'dask worker' | grep -v grep | awk '{print $2}' | xargs -r kill -9 && ps -ef | grep -w 'mhw_app' | grep -v grep | awk '{print $2}' | xargs -r kill -9

# pm2 start
pm2 start ./conf/ecosystem.config.js

