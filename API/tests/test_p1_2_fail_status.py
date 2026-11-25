"""Codex r22 (Major): a failed check must not look like success.

updater_check run / verify / api-check write the whole JSON first and then exit
3 when ok=false; check-result turns "ok as expected" into an exit status; and
append_chain.sh stops (exit 1, "CHAIN FAILED", diagnostics) on any failed step,
treating the negative control as passed only when it fails as built to.
"""
import importlib.util
import json
import os
import stat
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from tools.p1_2 import updater_check as uc

API = Path(__file__).resolve().parents[1]
UC = API / "tools" / "p1_2" / "updater_check.py"
CHAIN = API / "tools" / "p1_2" / "append_chain.sh"


# --------------------------------------------------------------- check-result


def _cr(tmp_path, doc, *args):
    f = tmp_path / "r.json"
    if doc is not None:
        f.write_text(doc if isinstance(doc, str) else json.dumps(doc))
    return uc.main(["check-result", "--file", str(f), *args])


@pytest.mark.parametrize("doc,expect,extra,rc", [
    ({"ok": True}, "ok", [], 0),
    ({"ok": False, "failed_checks": ["x"]}, "ok", [], uc.EXIT_UNEXPECTED),
    ({"ok": "true"}, "ok", [], uc.EXIT_UNEXPECTED),                      # not a JSON boolean
    ({"ok": 1}, "ok", [], uc.EXIT_UNEXPECTED),
    ({}, "ok", [], uc.EXIT_UNEXPECTED),
    (None, "ok", [], uc.EXIT_UNEXPECTED),                                # missing file
    ("{truncated", "ok", [], uc.EXIT_UNEXPECTED),                        # partial JSON
    # negative control: must fail, and for the reason it was built for
    ({"ok": False, "failed_checks": ["dates_last_is_fixture_month", "new_month_values"]}, "fail",
     ["--expect-failed", "dates_last_is_fixture_month"], 0),
    ({"ok": True, "failed_checks": []}, "fail", ["--expect-failed", "dates_last_is_fixture_month"],
     uc.EXIT_UNEXPECTED),                                                # it passed: not as expected
    ({"ok": False, "failed_checks": ["api"]}, "fail", ["--expect-failed", "dates_last_is_fixture_month"],
     uc.EXIT_UNEXPECTED),                                                # failed for another reason
    (None, "fail", [], uc.EXIT_UNEXPECTED),                              # crashed, no JSON
])
def test_check_result(tmp_path, capsys, doc, expect, extra, rc):
    assert _cr(tmp_path, doc, "--expect", expect, *extra) == rc
    verdict = json.loads(capsys.readouterr().out.strip().splitlines()[-1])
    assert verdict["as_expected"] is (rc == 0)


# --------------------------------------------------------------------- verify


def _good_verify_rep():
    vals = {v: {"nan_pattern_equal": True, "within_tol": True, "prev_month_unchanged": True} for v in uc.VARS}
    return {"files": {"removed": [], "unexplained": [], "metadata_diff": {}},
            "dates": {"n": 540, "last": "2026-12-01", "prev": "2026-11-01", "monotonic_unique": True,
                      "consecutive_months": True},
            "lat_lon_files_unchanged": True,
            "earlier_months": {"2026-11-01": {v: {"nan_pattern_equal": True, "within_tol": True} for v in uc.VARS}},
            "values": vals, "dtypes": dict(uc.DTYPES),
            "zarr_shapes": {v: [540, 720, 1440] for v in uc.VARS} | {"date": [540]},
            "consolidated_shapes": {v: [540, 720, 1440] for v in uc.VARS} | {"date": [540]},
            "api": {"ok": True, "skipped": True}, "n_before": 539, "before_store_given": True}


META2 = {"new_month": "2026-12-01", "last_month_in_store": "2026-11-01"}


def test_verify_checks_all_pass_on_a_good_report():
    rep = _good_verify_rep()
    uc.with_checks(rep, uc.verify_checks(rep, META2, 539))
    assert rep["ok"] is True and rep["failed_checks"] == [] and "date_count_is_before_plus_one" in rep["checks"]


@pytest.mark.parametrize("mutate,failed", [
    (lambda r: r["files"]["removed"].append("x"), "files_none_removed"),
    (lambda r: r["files"]["unexplained"].append("x"), "files_all_explained"),
    (lambda r: r["dates"].update(last="2026-11-01"), "dates_last_is_fixture_month"),
    (lambda r: r["dates"].update(prev="2026-10-01"), "dates_prev_is_fixture_base"),
    (lambda r: r["dates"].update(consecutive_months=False), "dates_consecutive_months"),
    (lambda r: r.update(lat_lon_files_unchanged=False), "lat_lon_files_unchanged"),
    (lambda r: r["earlier_months"]["2026-11-01"]["sst"].update(within_tol=False), "earlier_months_unchanged"),
    (lambda r: r["values"]["td"].update(nan_pattern_equal=False), "new_month_values"),
    (lambda r: r["values"]["level"].update(prev_month_unchanged=False), "prev_month_unchanged"),
    (lambda r: r["dtypes"].update(level="float64"), "dtypes"),
    (lambda r: r["zarr_shapes"].update(date=[539]), "zarr_shapes"),
    (lambda r: r["consolidated_shapes"].update(sst=[539, 720, 1440]), "consolidated_shapes"),
    (lambda r: r["api"].update(ok=False), "api"),
    (lambda r: r.update(n_before=538), "date_count_is_before_plus_one"),
])
def test_verify_checks_each_failure_is_named(mutate, failed):
    rep = _good_verify_rep()
    mutate(rep)
    uc.with_checks(rep, uc.verify_checks(rep, META2, rep["n_before"]))
    assert rep["ok"] is False and failed in rep["failed_checks"]


def test_verify_date_count_needs_the_before_shape():
    rep = _good_verify_rep()
    uc.with_checks(rep, uc.verify_checks(rep, META2, None))  # --before-store given, no date shape in it
    assert "date_count_is_before_plus_one" in rep["failed_checks"]


def _verify_argv(tmp_path):
    return ["verify", "--store", "s", "--fixtures", "f", "--before", "b", "--api", "a",
            "--out", str(tmp_path / "v.json")]


@pytest.mark.parametrize("ok,rc", [(True, 0), (False, uc.EXIT_RESULT_FAILED)])
def test_verify_command_writes_the_json_then_exits_by_ok(tmp_path, monkeypatch, ok, rc):
    def fake_report(a):
        rep = _good_verify_rep()
        if not ok:
            rep["dates"]["last"] = "2026-11-01"
        return uc.with_checks(rep, uc.verify_checks(rep, META2, 539))
    monkeypatch.setattr(uc, "_verify_report", fake_report)
    assert uc.main(_verify_argv(tmp_path)) == rc
    doc = json.loads((tmp_path / "v.json").read_text())       # complete JSON in both cases
    assert doc["ok"] is ok and ("dates_last_is_fixture_month" in doc["failed_checks"]) is (not ok)


# ------------------------------------------------------------------ api-check


@pytest.mark.parametrize("case,plot,rc", [
    ({"status": 200, "match": True}, {"status": 200, "png": True}, 0),
    ({"status": 200, "match": False}, {"status": 200, "png": True}, uc.EXIT_RESULT_FAILED),
    ({"status": 500}, {"status": 200, "png": True}, uc.EXIT_RESULT_FAILED),
    ({"status": 200, "match": True}, {"status": 200, "png": False}, uc.EXIT_RESULT_FAILED),
])
def test_api_check_command_exit_status(tmp_path, monkeypatch, case, plot, rc):
    pytest.importorskip("zarr")

    def fake_api_check(api_dir, store, fixture_dirs):
        res = {"cases": [dict(case)], "plot": dict(plot, content_type="image/png"), "months": ["2026-11-01"]}
        return uc.with_checks(res, {"points": all(c["status"] == 200 and c.get("match") is True
                                                  for c in res["cases"]),
                                    "plot": res["plot"]["status"] == 200 and res["plot"]["png"]})
    monkeypatch.setattr(uc, "_api_check", fake_api_check)
    out = tmp_path / "api.json"
    assert uc.main(["api-check", "--store", "s", "--fixtures", "f", "--api", "a", "--out", str(out)]) == rc
    assert json.loads(out.read_text())["ok"] is (rc == 0)


# ------------------------------------------------------------------------ run


def _good_run_rep():
    return {"error": None, "completed": True, "no_external_access": True, "network_attempts_blocked": 0,
            "calls": {"requests_get": 0},
            "audit": {"socket": [], "subprocess": [], "urllib": [], "dotenv_file_open": [], "write_outside": [],
                      "hook_errors": []}}


@pytest.mark.parametrize("mutate,failed", [
    (lambda r: r.update(error={"type": "ValueError"}, completed=False), "no_error"),
    (lambda r: r.update(completed=False), "completed"),
    (lambda r: r["audit"]["dotenv_file_open"].append("/x/.env"), "no_dotenv_open"),
    (lambda r: r["audit"]["socket"].append("socket.connect"), "no_socket_urllib_subprocess"),
    (lambda r: r["audit"]["subprocess"].append("subprocess.Popen"), "no_socket_urllib_subprocess"),
    (lambda r: r.update(network_attempts_blocked=1), "no_network_attempt"),
    (lambda r: r["calls"].update(requests_get=1), "no_network_attempt"),
    (lambda r: r["audit"]["write_outside"].append("/etc/x"), "no_write_outside"),
    (lambda r: r["audit"]["hook_errors"].append("open: TypeError"), "no_audit_hook_errors"),
])
def test_run_checks(mutate, failed):
    rep = _good_run_rep()
    assert uc.with_checks(dict(rep), uc.run_checks(rep))["ok"] is True
    mutate(rep)
    out = uc.with_checks(rep, uc.run_checks(rep))
    assert out["ok"] is False and failed in out["failed_checks"] and uc.result_exit(out) == uc.EXIT_RESULT_FAILED


_FAKE_UPDATER = '''
import os, socket
import zarr
from sqlalchemy import create_engine
ZARR_PATH = DATA_PATH = None
MODE, OUTSIDE = {mode!r}, {outside!r}


def append_to_zarr():
    open(os.path.join(ZARR_PATH, "chunk"), "w").close()      # a write inside the store
    if MODE == "raise":
        raise ValueError("append failed")
    if MODE == "socket":
        try:
            socket.create_connection(("192.0.2.1", 5432), timeout=0.01)
        except OSError:
            pass
    if MODE == "outside":
        open(os.path.join(OUTSIDE, "x"), "w").close()
    if MODE == "dotenv":
        open(os.path.join(ZARR_PATH, ".env"), "w").close()
    if MODE != "silent":
        print("All work append to zarr done")
'''

# the run's own imports, with stubs only for modules this venv lacks (V-a has no
# sqlalchemy); never the updater's network or DB code
_RUN_CHILD = r"""
import sys, types
sys.path.insert(0, sys.argv[1])
try:
    import sqlalchemy  # noqa: F401
except ImportError:
    m = types.ModuleType("sqlalchemy"); m.create_engine = lambda *a, **k: None; sys.modules["sqlalchemy"] = m
from tools.p1_2 import updater_check as uc
sys.exit(uc.main(sys.argv[2:]))
"""


@pytest.mark.parametrize("mode,rc,failed", [
    ("ok", 0, []),
    ("raise", 3, ["completed", "no_error"]),
    ("silent", 3, ["completed"]),
    ("socket", 3, ["no_external_access", "no_network_attempt"]),  # the patched connect refuses first
    ("outside", 3, ["no_external_access", "no_write_outside"]),
    ("dotenv", 3, ["no_dotenv_open", "no_external_access"]),
])
def test_run_command_end_to_end_exit_status(tmp_path, mode, rc, failed):
    for m in ("pandas", "pyarrow", "dotenv", "requests", "xarray", "zarr"):
        pytest.importorskip(m)
    import pandas as pd
    store, fix, api, outside, tmpdir = (tmp_path / n for n in ("store", "fix", "api", "outside", "tmpdir"))
    for d in (store, fix, api, outside, tmpdir):
        d.mkdir()
    pd.DataFrame({"lat": [0.125], "lon": [0.125]}).to_parquet(fix / "db.parquet")
    (api / "update_mhw.py").write_text(_FAKE_UPDATER.format(mode=mode, outside=str(outside)))
    out = tmp_path / "run.json"
    env = {"PATH": "/usr/bin:/bin", "HOME": str(tmp_path), "TMPDIR": str(tmpdir), "LANG": "C.UTF-8"}
    p = subprocess.run([sys.executable, "-B", "-c", _RUN_CHILD, str(API), "run", "--store", str(store),
                        "--fixtures", str(fix), "--api", str(api), "--out", str(out)],
                       capture_output=True, text=True, env=env, timeout=300)
    assert p.returncode == rc, p.stderr[-3000:]
    doc = json.loads(out.read_text())                         # written in every case
    assert doc["ok"] is (rc == 0) and doc["failed_checks"] == sorted(failed), doc["failed_checks"]
    # the isolated TMPDIR is what the audit allows, not the system temp directory
    assert doc["tempfile_gettempdir"] == os.path.realpath(tmpdir)
    assert os.path.realpath(outside) not in doc["allowed_write_roots"]
    if mode == "raise":
        assert doc["error"]["type"] == "ValueError"


# ------------------------------------------------------------ append_chain.sh

# Fake interpreters for append_chain.sh: the heavy updater_check commands, the
# loopback server and recover are simulated per label (no store, no network);
# check-result and the snapshot comparison run for real.
_FAKE = r'''#!{py}
import json, os, sys
REAL_PY, REAL_UC, CFG, CALLS = {py!r}, {uc!r}, {cfg!r}, {calls!r}
cfg = json.load(open(CFG))
argv = sys.argv[1:]
while argv and argv[0] in ("-I", "-B"):
    argv = argv[1:]
if argv[0] == "-c":
    os.execv(REAL_PY, [REAL_PY, "-I", "-B", *argv])


def opt(name):
    return argv[argv.index(name) + 1] if name in argv else None


def log(label):
    with open(CALLS, "a") as f:
        f.write(label + "\n")


def finish(label, doc=None, default_rc=0):
    c = cfg.get(label, {{}})
    doc = c.get("doc", doc)
    out = opt("--out")
    if doc is not None and out and not out.endswith(("meta0", "meta1", "fix1", "fix2")):
        with open(out, "w") as f:
            json.dump(doc, f)
    if doc is not None and not out:
        print(json.dumps(doc))
    sys.exit(c.get("rc", default_rc))


if argv[0] == "-m":                                        # tools interpreter
    mod = argv[1]
    if mod == "tools.p1_2.append_server":
        label = "server" + opt("--step")[-1]
        log(label); finish(label)
    if mod == "tools.bench.candidate":
        label = f"recover{{opt('--dir')[-1]}}-{{opt('--mode')}}"
        log(label); finish(label, {{"mode": opt("--mode"), "exit_status": 0}})
    sys.exit(99)
sub = argv[1]
if sub == "check-result":
    os.execv(REAL_PY, [REAL_PY, "-I", "-B", REAL_UC, *argv[1:]])
out = opt("--out")
label = os.path.basename(out).removesuffix(".json")
log(label)
if sub in ("meta-copy", "fixtures"):
    os.makedirs(out)
    finish(label)
if sub == "snapshot":
    finish(label, {{"a": ["sha", 1]}})
if sub == "run":
    finish(label, {{"ok": True, "failed_checks": [], "completed": True, "no_external_access": True}})
if sub == "verify":
    if label == "verify_negative_wrong_n":
        finish(label, {{"ok": False, "failed_checks": ["dates_last_is_fixture_month", "new_month_values"]}}, 3)
    finish(label, {{"ok": True, "failed_checks": []}})
finish(label, {{}})
'''


@pytest.fixture()
def chain(tmp_path):
    w, api, src, fk = (tmp_path / n for n in ("w", "api", "src", "fake"))
    for d in (w, api / "src", src, fk):
        d.mkdir(parents=True)
    (api / "mhw_app.py").write_text("")
    (src / ".zmetadata").write_text("{}")
    os.chmod(src / ".zmetadata", 0o444)
    os.chmod(src, 0o555)
    calls, cfg = tmp_path / "calls.txt", tmp_path / "cfg.json"
    for name in ("py", "tpy"):
        f = fk / name
        f.write_text(_FAKE.format(py=sys.executable, uc=str(UC), cfg=str(cfg), calls=str(calls)))
        f.chmod(f.stat().st_mode | stat.S_IXUSR)

    def run(config=None):
        cfg.write_text(json.dumps(config or {}))
        p = subprocess.run(["bash", str(CHAIN), str(w), str(fk / "py"), str(fk / "tpy"), str(api), str(src), "18099"],
                           capture_output=True, text=True, timeout=120)
        steps = (w / "logs" / "steps.txt").read_text()
        status = (w / "out" / "chain_status.txt").read_text().strip()
        return p.returncode, steps, status, calls.read_text().split()

    yield run
    os.chmod(src, 0o755)


def test_chain_all_steps_proven(chain):
    rc, steps, status, calls = chain()
    assert rc == 0 and status == "OK" and "CHAIN OK" in steps and "FAIL" not in steps
    assert "negative control: failed as expected" in steps
    assert calls[-2:] == ["verify_negative_wrong_n", "src_snap_after"]


def _failed(rc, steps, status, reason):
    assert rc == 1 and status.startswith("FAILED") and reason in status, (status, steps)
    assert "CHAIN FAILED" in steps and "CHAIN OK" not in steps


@pytest.mark.parametrize("config,reason,last,diag", [
    ({"verify1": {"rc": 3, "doc": {"ok": False, "failed_checks": ["new_month_values"]}}}, "verify1 rc=3", "verify1",
     []),
    ({"verify1": {"rc": 0, "doc": {"ok": False}}}, "verify1: out/verify1.json is not ok=true", "verify1", []),
    ({"run1": {"rc": 3, "doc": {"ok": False, "failed_checks": ["no_dotenv_open"]}}}, "run1 rc=3", "run1", []),
    # after the first server step, diagnose lists the API dir already used
    ({"verify2": {"rc": 1, "doc": None}}, "verify2 rc=1", "verify2", ["recover1-list"]),   # crashed, no JSON
    ({"run2": {"rc": 0, "doc": {"ok": False}}}, "run2: out/run2.json is not ok=true", "run2", ["recover1-list"]),
])
def test_chain_stops_at_a_failed_check(chain, config, reason, last, diag):
    rc, steps, status, calls = chain(config)
    _failed(rc, steps, status, reason)
    assert calls[calls.index(last):] == [last, *diag]  # nothing else after the failed step


@pytest.mark.parametrize("config,reason,anchor,tail", [
    ({"server1": {"rc": 1}}, "server1 rc=1", "server1", ["recover1-list"]),
    ({"server2": {"rc": 5}}, "server2 rc=5", "server2", ["recover1-list", "recover2-list"]),
    # the step's recover finds the candidate alive (exit 1): diagnose lists it
    # again and only then stops it (recover stop signals proven identities only)
    ({"recover2-list": {"rc": 1}}, "recover2 rc=1", "server2",
     ["recover2-list", "recover1-list", "recover2-list", "recover2-stop"]),
    # identity undetermined (exit 6): listed, never stopped
    ({"recover1-list": {"rc": 6}}, "recover1 rc=6", "server1", ["recover1-list", "recover1-list"]),
])
def test_chain_server_or_recover_failure_runs_diagnostics(chain, config, reason, anchor, tail):
    rc, steps, status, calls = chain(config)
    _failed(rc, steps, status, reason)
    assert calls[calls.index(anchor) + 1:] == tail
    assert "src_snap_after" not in calls and "verify_negative_wrong_n" not in calls


def test_chain_source_changed_fails(chain):
    rc, steps, status, calls = chain({"src_snap_after": {"doc": {"a": ["other", 1]}}})
    _failed(rc, steps, status, "source before vs after: DIFFERENT")


def test_chain_copy_differs_fails(chain):
    rc, steps, status, calls = chain({"snap0": {"doc": {"a": ["other", 1]}}})
    _failed(rc, steps, status, "copy vs source: DIFFERENT")
    assert calls[-1] == "snap0"


@pytest.mark.parametrize("neg,reason", [
    ({"rc": 0, "doc": {"ok": True, "failed_checks": []}}, "negative control: verify rc=0"),          # passed
    ({"rc": 1, "doc": None}, "negative control: verify rc=1"),                                       # crashed
    ({"rc": 3, "doc": {"ok": False, "failed_checks": ["api"]}}, "negative control: not the expected failure"),
    ({"rc": 3, "doc": {"ok": True, "failed_checks": []}}, "negative control: not the expected failure"),
])
def test_chain_negative_control_must_fail_as_built(chain, neg, reason):
    rc, steps, status, calls = chain({"verify_negative_wrong_n": neg})
    _failed(rc, steps, status, reason)
    assert "src_snap_after" not in calls


def test_chain_r22_record_would_now_fail():
    """The r22 record: negative verify wrote ok=false but the chain logged rc=0;
    with exit status 3 for ok=false that rc=0 can no longer pass the chain."""
    assert uc.result_exit({"ok": False}) == uc.EXIT_RESULT_FAILED and uc.result_exit({"ok": True}) == 0
    assert uc.result_exit({}) == uc.EXIT_RESULT_FAILED and uc.result_exit({"ok": "true"}) == uc.EXIT_RESULT_FAILED
