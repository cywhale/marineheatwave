# Sourced by conf/start_app.sh (API) and conf/run_update_mhw.sh (monthly updater):
# both use ONE uv venv built from requirements/p1-2/runtime.lock (owner 2026-11-25;
# doc 29 r5 D1).  In production since 2026-11-25.
#
# MHW_VENV: the default below is the production venv path; an explicit MHW_VENV
# in the environment wins.
: "${MHW_VENV:=/home/odbadmin/venvs/mhw-p12}"
export MHW_VENV
# polars[rtcompat] installs the "compat" runtime; force it so a missing compat
# runtime fails at import instead of falling back to the AVX2 build (doc 28 r4 §9.3)
export POLARS_FORCE_PKG=compat
# the API no longer uses the shared Dask scheduler (P1-1a): unset the four Dask
# settings known from the old launcher environment.  Only these four; any other
# DASK_* variable (e.g. DASK_CONFIG) is NOT cleared (doc 28 r6 §11.4).
unset DASK_SCHEDULER_ADDRESS DASK_TLS_CA_FILE DASK_TLS_CLIENT_CERT DASK_TLS_CLIENT_KEY
if [ ! -x "$MHW_VENV/bin/python" ]; then
  echo "p12_env.sh: no python in MHW_VENV=$MHW_VENV" >&2
  exit 1
fi
