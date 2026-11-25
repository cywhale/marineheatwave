"""P1-2 launch configuration (doc 29 r5 D1/D2): one uv venv for the API and the
updater, exec, uvicorn-worker, no --reload, POLARS_FORCE_PKG=compat, a PM2 app
with explicit cwd / kill_timeout and no pre_stop.  The scripts run against a
fake venv (fake gunicorn / python that only report what they received)."""
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

API = Path(__file__).resolve().parents[1]
REQ = API / "requirements" / "p1-2"

FAKE = r"""#!/bin/bash
python3 - "$$" "$PWD" "$@" <<'PY'
import json, os, sys
print(json.dumps({"pid": int(sys.argv[1]), "cwd": sys.argv[2], "argv": sys.argv[3:],
                  "env": {k: os.environ.get(k) for k in ("POLARS_FORCE_PKG", "MHW_VENV", "DASK_SCHEDULER_ADDRESS",
                                                         "DASK_TLS_CA_FILE", "DASK_TLS_CLIENT_CERT",
                                                         "DASK_TLS_CLIENT_KEY")}}))
PY
"""


@pytest.fixture()
def venv(tmp_path):
    v = tmp_path / "venv"
    (v / "bin").mkdir(parents=True)
    for name in ("gunicorn", "python"):
        f = v / "bin" / name
        f.write_text(FAKE)
        f.chmod(0o755)
    return v


def _run(script, venv, extra_env=None):
    env = {"PATH": "/usr/bin:/bin", "MHW_VENV": str(venv), "DASK_SCHEDULER_ADDRESS": "tls://x:8786",
           "DASK_TLS_CA_FILE": "/x/ca", "DASK_TLS_CLIENT_CERT": "/x/c", "DASK_TLS_CLIENT_KEY": "/x/k",
           **(extra_env or {})}
    p = subprocess.Popen([str(API / "conf" / script)], cwd="/", env=env, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True)
    out, err = p.communicate(timeout=30)
    return p, out, err


def test_api_launcher_execs_gunicorn_from_the_shared_venv(venv):
    p, out, err = _run("start_app.sh", venv)
    assert p.returncode == 0, err
    got = json.loads(out)
    assert got["pid"] == p.pid  # exec: PM2's process IS the gunicorn master
    assert got["cwd"] == str(API)
    a = got["argv"]
    assert a[:1] == ["mhw_app:app"] and a[a.index("-w") + 1] == "2"
    assert a[a.index("-k") + 1] == "uvicorn_worker.UvicornWorker" and "--reload" not in a
    assert a[a.index("-b") + 1] == "127.0.0.1:8030"
    assert got["env"]["POLARS_FORCE_PKG"] == "compat" and got["env"]["MHW_VENV"] == str(venv)
    assert all(got["env"][k] is None for k in got["env"] if k.startswith("DASK"))


def test_updater_launcher_uses_the_same_venv_and_settings(venv):
    p, out, err = _run("run_update_mhw.sh", venv)
    assert p.returncode == 0, err
    got = json.loads(out)
    assert got["pid"] == p.pid and got["cwd"] == str(API) and got["argv"] == ["update_mhw.py"]
    assert got["env"]["POLARS_FORCE_PKG"] == "compat" and got["env"]["MHW_VENV"] == str(venv)
    assert all(got["env"][k] is None for k in got["env"] if k.startswith("DASK"))


@pytest.mark.parametrize("script", ["start_app.sh", "run_update_mhw.sh"])
def test_missing_venv_fails_before_starting_anything(tmp_path, script):
    p, out, err = _run(script, tmp_path / "nowhere")
    assert p.returncode == 1 and out == "" and "no python in MHW_VENV" in err


def test_one_venv_definition_for_both_launchers():
    for script in ("start_app.sh", "run_update_mhw.sh"):
        sh = (API / "conf" / script).read_text()
        assert "source conf/p12_env.sh" in sh and "set -euo pipefail" in sh
        assert "pyenv" not in sh and "py311" not in sh and "MHW_VENV=" not in sh
    assert (API / "conf/p12_env.sh").read_text().count("MHW_VENV:=") == 1


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_ecosystem_manages_only_mhwapi_with_explicit_settings():
    js = subprocess.run(["node", "-e", "process.stdout.write(JSON.stringify(require(process.argv[1])))",
                         str(API / "conf/ecosystem.config.js")], capture_output=True, text=True, check=True)
    apps = json.loads(js.stdout)["apps"]
    assert [a["name"] for a in apps] == ["mhwapi"]
    a = apps[0]
    assert a["cwd"] == "/home/odbadmin/proj/marineheatwave/API" and a["script"] == "./conf/start_app.sh"
    assert a["interpreter"] == "none" and a["kill_timeout"] >= 10000
    assert a["append_env_to_name"] is False and "pre_stop" not in a and a["watch"] is False
    assert not any("dask" in str(v).lower() for app in apps for v in app.values())


def test_official_requirements_are_the_validated_set():
    sums = dict(reversed(l.split()) for l in (REQ / "SHA256SUMS").read_text().splitlines())
    for name, digest in sums.items():
        assert hashlib.sha256((REQ / name).read_bytes()).hexdigest() == digest, name
    # the locks validated in V-c (vm24, B) and V-a (dev)
    assert sums["runtime.lock"].startswith("73d273ce") and sums["dev.lock"].startswith("451f225e")
    rin = [l.split("#")[0].strip() for l in (REQ / "runtime.in").read_text().splitlines()]
    for need in ("polars[rtcompat]", "psycopg[binary]", "sqlalchemy", "dask[array]", "netcdf4", "uvicorn-worker"):
        assert need in rin, need
    assert "polars" not in rin  # only the rtcompat extra


def _pins(path):
    return {l.split("==")[0]: l.split("==")[1].split()[0] for l in path.read_text().splitlines()
            if l and l[0].isalnum() and "==" in l}


def test_locks_agree_and_carry_the_compat_runtime():
    rt, dev = _pins(REQ / "runtime.lock"), _pins(REQ / "dev.lock")
    assert rt["polars"] == rt["polars-runtime-compat"] and "psycopg-binary" in rt and "psycopg2-binary" not in rt
    shared = set(rt) & set(dev)
    assert shared and all(rt[k] == dev[k] for k in shared)
    assert {k: v for k, v in _pins(REQ / "tools.lock").items() if k in dev} == \
        {k: dev[k] for k in _pins(REQ / "tools.lock") if k in dev}


def test_old_dependency_files_are_kept():
    assert (API / "requirements.txt").exists() and (API / "Pipfile").exists()
