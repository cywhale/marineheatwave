module.exports = {
  apps : [
  {
    name: 'dask-scheduler',
    script: 'dask-scheduler &',
    autorestart: true,
    restart_delay: 5000,
  },
  {
    name: 'dask-worker',
    script: 'dask-worker tcp://localhost:8786 &',
    autorestart: true,
    restart_delay: 10000,
  },
  {
    name: 'mhwapi',
    script: 'gunicorn mhw_app:app -w 4 -k uvicorn.workers.UvicornWorker -b 127.0.0.1:8020 --keyfile conf/privkey.pem --certfile conf/fullchain.pem --reload',
    args: '',
    merge_logs: true,
    autorestart: true,
    log_file: "tmp/combined.outerr.log",
    out_file: "tmp/out.log",
    error_file: "tmp/err.log",
    log_date_format : "YYYY-MM-DD HH:mm Z",
    append_env_to_name: true,
    watch: false,
    max_memory_restart: '4G',
    pre_stop: 'kill $(ps aux | grep \'mhw_app\' | awk \'{print $2}\') && kill $(ps aux | grep \'dask-scheduler\' | awk \'{print $2}\') && kill $(ps aux | grep \'dask-worker\' | awk \'{print $2}\')',
  }],
};
