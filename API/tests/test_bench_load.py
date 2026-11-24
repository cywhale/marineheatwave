"""Failure handling of tools/bench/load.py against a local mock server (Codex round 4 #1).

No real MHW server or data is needed.  Every failure mode must make the run
fail (LoadError / non-zero exit / FAILED file) and nothing may hang.
"""
import json
import os
import socket
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from tools.bench import load


class _Handler(BaseHTTPRequestHandler):
    def do_GET(self):  # noqa: N802
        if self.path.startswith("/slow"):
            time.sleep(1.5)
        if self.path.startswith("/fail"):
            code, body = 500, b'{"detail":"boom"}'
        else:
            code, body = 200, b'{"ok":true}'
        if self.path.startswith("/heavy"):
            time.sleep(0.3)
        self.send_response(code)
        self.send_header("content-type", "application/json")
        self.send_header("content-length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass


@pytest.fixture(scope="module")
def server():
    srv = ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
    srv.daemon_threads = True
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    yield f"http://127.0.0.1:{srv.server_address[1]}"
    srv.shutdown()


@pytest.fixture()
def dead_base():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()  # nothing listens on this port any more
    return f"http://127.0.0.1:{port}"


def test_success_counts_and_real_start_end(server):
    res = load.run_threads(server, [["/ok"] * 3, ["/ok"] * 2], timeout=5)
    assert len(res) == 5 and all(r["status"] == 200 for r in res)
    assert all(r["start"] <= r["end"] for r in res)
    # start is taken before the request: the slow path must show a real gap
    slow = load.run_threads(server, [["/slow"]], timeout=5)[0]
    assert slow["end"] - slow["start"] >= 1.4


def test_http_500_fails(server):
    with pytest.raises(load.LoadError) as e:
        load.run_threads(server, [["/ok", "/fail"], ["/ok"]], timeout=5)
    assert any(d.get("status") == 500 for d in e.value.details)


def test_timeout_fails_without_hanging(server):
    t0 = time.perf_counter()
    with pytest.raises(load.LoadError) as e:
        load.run_threads(server, [["/slow"], ["/ok"]], timeout=0.2)
    assert time.perf_counter() - t0 < 3
    assert any(d["kind"] == "exception" and "Timeout" in d["type"] for d in e.value.details)
    assert any(d["kind"] == "count" for d in e.value.details)


def test_connection_failure_fails(dead_base):
    with pytest.raises(load.LoadError) as e:
        load.run_threads(dead_base, [["/ok"], ["/ok"]], timeout=1)
    assert sum(d["kind"] == "exception" for d in e.value.details) == 2


def test_mixed_success(server):
    res = load.scenario_mixed(server, None, heavy_path="/heavy", light_path="/ok", timeout=5)
    assert res["heavy"]["status"] == 200
    assert len(res["light"]) >= 9


@pytest.mark.parametrize("heavy_path", ["/fail", "/slow"])
def test_mixed_heavy_failure_stops_light_threads(server, heavy_path):
    t0 = time.perf_counter()
    with pytest.raises(load.LoadError):
        load.scenario_mixed(server, None, heavy_path=heavy_path, light_path="/ok", timeout=0.5)
    assert time.perf_counter() - t0 < 5


def test_mixed_heavy_connection_failure(dead_base):
    t0 = time.perf_counter()
    with pytest.raises(load.LoadError):
        load.scenario_mixed(dead_base, None, heavy_path="/heavy", light_path="/ok", timeout=1)
    assert time.perf_counter() - t0 < 5


def test_light_500_fails_mixed(server):
    with pytest.raises(load.LoadError):
        load.scenario_mixed(server, None, heavy_path="/heavy", light_path="/fail", timeout=5)


# A dedicated "master" process stands in for gunicorn: the pytest process itself
# may own unrelated children (e.g. multiprocessing.resource_tracker after the
# golden tests ran), which the sampler correctly reports as unclassified.
_FAKE_MASTER = r"""
import os, subprocess, sys, time
n, odd, replace_after = int(sys.argv[1]), sys.argv[2] == "1", float(sys.argv[3])

def fork_worker():
    pid = os.fork()
    if pid == 0:
        time.sleep(30)
        os._exit(0)
    return pid

workers = [fork_worker() for _ in range(n)]
if odd:  # a direct child with a different title, which itself has a child
    subprocess.Popen([sys.executable, "-c",
                      "import subprocess, sys, time; "
                      "subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(30)']); time.sleep(30)"])
print(" ".join(map(str, workers)), flush=True)
if replace_after > 0:
    time.sleep(replace_after)
    os.kill(workers[0], 9)
    os.waitpid(workers[0], 0)
    print(fork_worker(), flush=True)
time.sleep(30)
"""


@pytest.fixture()
def fake_master():
    import signal
    import subprocess
    import sys

    started = []

    def start(workers=0, odd_child=False, replace_after=0.0):
        proc = subprocess.Popen([sys.executable, "-c", _FAKE_MASTER, str(workers), "1" if odd_child else "0",
                                 str(replace_after)], stdout=subprocess.PIPE, text=True, start_new_session=True)
        started.append(proc)
        line = proc.stdout.readline().split()
        time.sleep(0.5 if odd_child else 0.1)  # let the odd child start its grandchild
        return proc, [int(x) for x in line]

    yield start
    for proc in started:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


def _ok(base, sampler, timeout):
    return {"r": load.run_threads(base, [["/ok"]], timeout=timeout)}


def _read_failed(tmp_path, name="light"):
    assert not (tmp_path / f"local_{name}.json").exists()
    return json.loads((tmp_path / f"local_{name}.FAILED.json").read_text())


def test_run_scenario_failure_writes_failed_file_and_stops_sampler(server, tmp_path, fake_master):
    def broken(base, sampler, timeout):
        return {"r": load.run_threads(base, [["/fail"]], timeout=timeout)}

    master, _ = fake_master()
    rc = load.run_scenario("light", server, master.pid, tmp_path, timeout=2, scenario_fn=broken, expected_workers=0)
    assert rc == 1
    doc = _read_failed(tmp_path)
    assert doc["valid"] is False and doc["http_valid"] is False and doc["memory_valid"] is True
    assert doc["sampler_alive_after_stop"] is False
    assert doc["error"]["type"] == "LoadError"


def test_run_scenario_success_records_process_roles(server, tmp_path, fake_master):
    master, workers = fake_master(workers=2)
    rc = load.run_scenario("light", server, master.pid, tmp_path, timeout=2, scenario_fn=_ok, expected_workers=2)
    assert rc == 0
    doc = json.loads((tmp_path / "local_light.json").read_text())
    procs = {p["pid"]: p for p in doc["processes"]}
    assert procs[master.pid]["role"] == "master" and procs[master.pid]["ppid"] == os.getpid()
    for w in workers:
        assert procs[w]["role"] == "worker" and procs[w]["ppid"] == master.pid
        assert procs[w]["cmdline"] == procs[master.pid]["cmdline"] and procs[w]["name"]
    assert doc["sampler_alive_after_stop"] is False
    info = doc["sampler"]
    assert info["n_samples"] >= 1 and info["error"] is None
    assert info["started_at"] <= info["first_sample_t"] <= info["last_sample_t"] <= info["ended_at"]
    assert info["observed_workers"] == sorted(workers)
    assert "not a total scenario deadline" in doc["note"]


def test_run_scenario_without_expected_workers_is_invalid(server, tmp_path, fake_master):
    master, _ = fake_master(workers=1)
    rc = load.run_scenario("light", server, master.pid, tmp_path, timeout=2, scenario_fn=_ok)
    assert rc == 1
    doc = _read_failed(tmp_path)
    assert doc["http_valid"] is True and doc["memory_valid"] is False
    assert [p["kind"] for p in doc["memory_problems"]] == ["expected_workers_missing"]


def test_cli_requires_expected_workers(tmp_path):
    with pytest.raises(SystemExit):
        load.main(["--master-pid", str(os.getpid()), "--scenario", "light", "--out", str(tmp_path)])


# ------------------------------------------------------------ sampler failures


def test_master_vanishing_mid_run_fails_with_reason(server, tmp_path, fake_master):
    """HTTP succeeds but the sampled process exits: the memory result must be invalid."""
    import signal

    master, _ = fake_master(workers=1)

    def ok_then_master_dies(base, sampler, timeout):
        r = load.run_threads(base, [["/ok"]], timeout=timeout)
        time.sleep(0.3)
        os.killpg(master.pid, signal.SIGKILL)
        master.wait()
        time.sleep(0.5)  # let the sampler hit the vanished master
        return {"r": r}

    rc = load.run_scenario("light", server, master.pid, tmp_path, timeout=2, scenario_fn=ok_then_master_dies,
                           expected_workers=1)
    assert rc == 1
    doc = _read_failed(tmp_path)
    assert doc["http_valid"] is True and doc["memory_valid"] is False
    err = [p for p in doc["memory_problems"] if p["kind"] == "sampler_error"]
    assert err and err[0]["type"] in ("NoSuchProcess", "ZombieProcess")
    assert doc["sampler"]["error"]["samples_before_error"] >= 1


def test_sampler_exception_mid_run_fails_with_reason(server, tmp_path, monkeypatch, fake_master):
    real = load.RssSampler.snapshot
    calls = {"n": 0}

    def flaky(self):
        calls["n"] += 1
        if calls["n"] >= 3:
            raise RuntimeError("simulated sampler fault")
        return real(self)

    monkeypatch.setattr(load.RssSampler, "snapshot", flaky)

    def wait_ok(base, sampler, timeout):
        r = load.run_threads(base, [["/ok"]], timeout=timeout)
        time.sleep(0.8)
        return {"r": r}

    master, _ = fake_master()
    rc = load.run_scenario("light", server, master.pid, tmp_path, timeout=2, scenario_fn=wait_ok,
                           expected_workers=0)
    assert rc == 1
    doc = _read_failed(tmp_path)
    assert doc["http_valid"] is True and doc["memory_valid"] is False
    err = [p for p in doc["memory_problems"] if p["kind"] == "sampler_error"][0]
    assert err["type"] == "RuntimeError" and "simulated sampler fault" in err["message"]


def test_sampler_stop_timeout_fails(server, tmp_path, monkeypatch, fake_master):
    real = load.RssSampler.snapshot

    def stuck(self):
        if threading.current_thread() is self:
            time.sleep(1.5)
        return real(self)

    monkeypatch.setattr(load.RssSampler, "snapshot", stuck)
    master, _ = fake_master()
    t0 = time.perf_counter()
    rc = load.run_scenario("light", server, master.pid, tmp_path, timeout=2, scenario_fn=_ok,
                           expected_workers=0, sampler_kwargs={"stop_timeout": 0.1})
    assert time.perf_counter() - t0 < 3
    assert rc == 1
    doc = _read_failed(tmp_path)
    assert doc["http_valid"] is True and doc["memory_valid"] is False
    assert any(p["kind"] == "sampler_stop_timeout" for p in doc["memory_problems"])
    assert doc["sampler"]["stop_timed_out"] is True


# ------------------------------------------------------- worker classification


def _sample_for(sampler, seconds=0.5):
    sampler.start()
    time.sleep(seconds)
    sampler.stop()
    return sampler


def test_role_classification(fake_master):
    """Same-cmdline direct child -> worker; different-title direct child -> unclassified;
    grandchild -> helper.  PID/PPID/name/cmdline are kept for manual checks."""
    master, workers = fake_master(workers=1, odd_child=True)
    s = load.RssSampler(master.pid)
    s.snapshot()
    procs = s.processes
    assert procs[master.pid]["role"] == "master"
    assert procs[workers[0]]["role"] == "worker"
    odd = [p for p in procs.values() if p["ppid"] == master.pid and p["pid"] not in workers]
    assert len(odd) == 1 and odd[0]["role"] == "unclassified"
    grand = [p for p in procs.values() if p["ppid"] == odd[0]["pid"]]
    assert len(grand) == 1 and grand[0]["role"] == "helper"
    for p in procs.values():
        assert set(p) >= {"pid", "ppid", "name", "cmdline", "role"}


def test_expected_worker_count_matches(fake_master):
    master, workers = fake_master(workers=2)
    s = _sample_for(load.RssSampler(master.pid, interval=0.05))
    assert s.observed_workers() == sorted(workers)
    assert s.memory_problems(2) == []


def test_worker_count_mismatch_is_invalid(fake_master):
    master, _ = fake_master(workers=1)
    s = _sample_for(load.RssSampler(master.pid, interval=0.05))
    kinds = [p["kind"] for p in s.memory_problems(2)]
    assert "worker_count" in kinds and "worker_count_in_sample" in kinds


def test_different_process_title_is_invalid(fake_master):
    """A direct child whose title differs from the master (a renamed worker or an
    unknown helper) is not silently dropped: it is unclassified and invalidates
    the memory result, with its PID/cmdline kept for manual checking."""
    master, workers = fake_master(workers=1, odd_child=True)
    s = _sample_for(load.RssSampler(master.pid, interval=0.05))
    probs = s.memory_problems(1)
    unc = [p for p in probs if p["kind"] == "unclassified_child"]
    assert unc and len(unc[0]["processes"]) == 1
    odd = unc[0]["processes"][0]
    assert odd["ppid"] == master.pid and odd["pid"] not in workers and odd["cmdline"]


def test_setproctitle_style_worker_is_recognised():
    s = load.RssSampler(os.getpid())
    assert s._role(123, os.getpid(), ["gunicorn: worker [tools.bench.isolated_app:app]"]) == "worker"
    assert s._role(123, os.getpid(), ["gunicorn: master [x]"]) == "unclassified"
    assert s._role(123, 1, s.master_cmdline) == "helper"


def test_worker_replaced_mid_run_is_invalid(fake_master):
    """Same count at the end but a different PID part-way (e.g. a restarted worker)."""
    master, workers = fake_master(workers=1, replace_after=0.4)
    s = load.RssSampler(master.pid, interval=0.05)
    s.start()
    new_pid = int(master.stdout.readline())
    time.sleep(0.3)
    s.stop()
    assert s.observed_workers() == sorted([workers[0], new_pid])
    assert "worker_count" in [p["kind"] for p in s.memory_problems(1)]
