module.exports = {
  apps : [
/*{
    name: 'dask-scheduler',
    script: 'dask-scheduler',
    autorestart: true,
    restart_delay: 5000,
  },
  {
    name: 'dask-worker',
    script: 'dask-worker tcp://localhost:8786',
    autorestart: true,
    restart_delay: 10000,
  },*/
  {
    name: 'mhwapi',
    script: './conf/start_app.sh',
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
    pre_stop: "ps -ef | grep -w 'mhw_app' | grep -v grep | awk '{print $2}' | xargs -r kill -9"
  }],
};
