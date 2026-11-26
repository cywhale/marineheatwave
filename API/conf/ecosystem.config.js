// PM2 app for the MHW API, P1-2 (doc 29 r5 D2).  In production since 2026-11-25.
// Only mhwapi: the shared Dask PM2 services are not managed here.
module.exports = {
  apps: [
    {
      name: 'mhwapi',
      cwd: '/home/odbadmin/proj/marineheatwave/API',
      script: './conf/start_app.sh',
      interpreter: 'none',       // the script's own #!/bin/bash; it execs gunicorn
      args: '',
      merge_logs: true,
      autorestart: true,
      log_file: 'tmp/combined.outerr.log',
      out_file: 'tmp/out.log',
      error_file: 'tmp/err.log',
      log_date_format: 'YYYY-MM-DD HH:mm Z',
      append_env_to_name: false,
      watch: false,
      max_memory_restart: '4G',
      kill_timeout: 20000,       // ms: PM2's stop signal (default SIGINT) reaches the gunicorn master directly
                                 // because start_app.sh execs; SIGKILL only after this timeout
      // pre_stop removed (doc 29 D2: not an effective PM2 setting; it was a kill -9 of every mhw_app process)
    },
  ],
};
