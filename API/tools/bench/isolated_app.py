"""Benchmark-only gunicorn entry point (spec 16 r3 §5.1, Codex round 3 #3).

    cd API && .venv/bin/gunicorn tools.bench.isolated_app:app -w 2 \
        -k uvicorn.workers.UvicornWorker -b 127.0.0.1:18030 --timeout 180

pytest's conftest does not reach separately started gunicorn workers, so this
module applies the same isolation (tests/harness.py) at import time.  Without
--preload every worker imports it, so every worker is isolated.  Each worker
prints its environment as one JSON line on stderr.  Production code is not
changed.
"""
import json
import sys

from tests.harness import isolate_env

isolate_env()  # before dask / distributed / matplotlib are imported

from tests.harness import describe_env, patch_app  # noqa: E402

app = patch_app()
print("MHW_BENCH_WORKER_ENV " + json.dumps(describe_env()), file=sys.stderr, flush=True)
