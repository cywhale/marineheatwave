"""Pure parts of tools/p1_2/vc_server.py (P1-2 V-c).  The servers themselves are
exercised by the candidate-runner tests and by the V-c run."""
import hashlib
import http.server
import socket
import threading

import pytest

from tools.bench import load
from tools.p1_2 import vc_server as vc


def test_runtime_of_maps():
    sp = "/v/lib/python3.14/site-packages"
    assert vc.runtime_of([f"{sp}/_polars_runtime_compat/_polars_runtime.abi3.so"]) == "compat"
    assert vc.runtime_of([f"{sp}/_polars_runtime_32/_polars_runtime.abi3.so"]) == "32"
    assert vc.runtime_of([]) is None  # nothing mapped
    assert vc.runtime_of([f"{sp}/_polars_runtime_32/x.so", f"{sp}/_polars_runtime_compat/x.so"]) is None  # both
    assert vc.runtime_of("unreadable: PermissionError") is None


def test_polars_maps_of_a_missing_process():
    out = vc.polars_maps([2 ** 22 + 12345])
    assert list(out.values())[0].startswith("unreadable")


def test_bootstrap_ci_brackets_the_median():
    d = [0.001 * i for i in range(-10, 21)]
    lo, hi = vc.bootstrap_ci(d, n=2000)
    assert lo <= sorted(d)[len(d) // 2] <= hi
    assert vc.bootstrap_ci(d, n=2000) == (lo, hi)  # seeded


def _serve(body: bytes):
    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *a):
            pass

    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, f"http://127.0.0.1:{srv.server_address[1]}"


def test_measure_paired_rotates_and_checks_bodies(monkeypatch):
    monkeypatch.setattr(load, "GOLDEN_J01_SHA256", hashlib.sha256(b"same").hexdigest())
    s1, b1 = _serve(b"same")
    s2, b2 = _serve(b"same")
    s3, b3 = _serve(b"same")
    try:
        m = vc.measure_paired({"prod": b1, "A": b2, "B": b3}, rounds=6)
    finally:
        for s in (s1, s2, s3):
            s.shutdown()
    assert m["problems"] == []
    meas = [r for r in m["rows"] if r["phase"] == "measure"]
    assert len(meas) == 6 * 2 * 3
    firsts = [next(r["label"] for r in meas if r["round"] == k and r["query"] == "Q1") for k in range(6)]
    assert sorted(firsts) == ["A", "A", "B", "B", "prod", "prod"]  # each label first equally often
    assert set(m["stats"]["pairs"]["Q1"]) == {"A-prod", "B-prod", "B-A"}
    assert m["stats"]["per_label"]["Q3"]["B"]["n"] == 6


def test_measure_paired_flags_a_differing_body(monkeypatch):
    monkeypatch.setattr(load, "GOLDEN_J01_SHA256", hashlib.sha256(b"same").hexdigest())
    s1, b1 = _serve(b"same")
    s2, b2 = _serve(b"different")
    try:
        m = vc.measure_paired({"A": b1, "B": b2}, rounds=2)
    finally:
        s1.shutdown()
        s2.shutdown()
    assert len([p for p in m["problems"] if "rows" in p]) == 4


def test_gate_reports_a_bound_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    s.listen(1)
    try:
        assert any("in use" in p for p in vc.gate(s.getsockname()[1]))
    finally:
        s.close()


@pytest.mark.parametrize("step", ["echo", "golden", "smoke", "combo"])
def test_scenario_output_name_follows_the_step_suffix(monkeypatch, tmp_path, step):
    """V-c deviation: golden and golden_forcecompat wrote the same out/local_*.json."""
    seen = []
    monkeypatch.setattr(vc, "gate", lambda port: [])
    monkeypatch.setattr(vc.cand, "run_step", lambda D, name, port, **k: seen.append((name, k["scenario_name"])) or 0)
    import argparse
    for suffix in ("", "_forcecompat"):
        a = argparse.Namespace(port=1, python="/p", gunicorn="/g", worker_class="w.W", step_suffix=suffix)
        vc.run_one(tmp_path, step, a, {}, None)
    (n1, s1), (n2, s2) = seen
    assert n2 == n1 + "_forcecompat" and s2 == s1 + "_forcecompat" and s1 != s2
