"""tools/p2/probe.py against a local fake HTTP server only (MHW Phase 2, Stage 0).

Every socket connect is checked: anything other than loopback fails the test, so
these tests cannot reach a real host.  The fake server records each request line,
which shows what the tool actually sent (no retries, no added parameters).
"""
import json
import os
import socket
import statistics
import struct
import subprocess
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

from tools.p2 import probe

API_DIR = Path(__file__).resolve().parents[1]
SECRET = "session-secret-value-123"


def fake_png(w=256, h=128, pad=1000):
    return probe.PNG_SIG + struct.pack(">I", 13) + b"IHDR" + struct.pack(">II", w, h) + b"\x08\x06\x00\x00\x00" + b"\0" * 4 + b"x" * pad


class _Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    seen: list = []

    def do_GET(self):  # noqa: N802
        type(self).seen.append(self.path)
        q = {k: v[0] for k, v in parse_qs(urlsplit(self.path).query, keep_blank_values=True).items()}
        mode = q.get("mode", "png")
        time.sleep(int(q.get("delay_ms", "0")) / 1000)
        if mode == "drop":
            self.close_connection = True
            return
        status = int(q.get("status", "200"))
        if mode == "xml":
            body, ctype = b"<ServiceExceptionReport/>", "application/vnd.ogc.se_xml"
        else:
            body, ctype = fake_png(int(q.get("w", 256)), int(q.get("h", 128)), int(q.get("pad", 1000))), "image/png"
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(1000 if mode == "partial" else len(body)))
        self.send_header("geowebcache-cache-result", q.get("gwc", "MISS"))
        self.send_header("X-Test-Cache", "BYPASS")
        self.send_header("Set-Cookie", f"sid={SECRET}; Path=/")
        self.send_header("Connection", "close")
        self.end_headers()
        if mode == "partial":
            self.wfile.write(b"0123456789")
            self.close_connection = True
            return
        time.sleep(int(q.get("body_delay_ms", "0")) / 1000)
        self.wfile.write(body)

    def log_message(self, *args):
        pass


@pytest.fixture(autouse=True)
def loopback_only(monkeypatch):
    real = socket.socket.connect

    def guarded(self, address):
        host = address[0] if isinstance(address, tuple) else address
        if host not in ("127.0.0.1", "::1", "localhost"):
            raise AssertionError(f"test tried to connect to non-loopback address {address!r}")
        return real(self, address)

    monkeypatch.setattr(socket.socket, "connect", guarded)


@pytest.fixture()
def server():
    _Handler.seen = []
    srv = ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
    srv.daemon_threads = True
    srv.handle_error = lambda *a: None  # the client may close first (timeout, drop): not a test failure
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}", _Handler.seen
    srv.shutdown()
    srv.server_close()


def manifest(tmp_path, base, *, cases, arms=None, design=None, cap=100, timeout=10, **extra):
    d = {
        "schema": probe.SCHEMA, "batch": "local-test", "max_requests": cap, "timeout_s": timeout,
        "targets": {"fake": {"base_url": base + "/geoserver/marineheatwave/wms"}},
        "arms": arms or {"A": {"target": "fake"}},
        "cases": cases,
        "design": design or {"warmup": 1, "repeats": 3},
        **extra,
    }
    if cap is None:
        del d["max_requests"]
    p = tmp_path / "manifest.json"
    p.write_text(json.dumps(d))
    return p


def read_run(out):
    recs = [json.loads(line) for line in (out / "requests.jsonl").read_text().splitlines()]
    return recs, json.loads((out / "summary.json").read_text())


# ------------------------------------------------------------------ timing accuracy


def test_injected_header_delay_is_measured_within_5pct(server, tmp_path):
    base, _ = server
    p = manifest(tmp_path, base, cases=[{"id": "d300", "query": "SERVICE=WMS&REQUEST=GetMap&delay_ms=300"}],
                 design={"warmup": 1, "repeats": 7})
    rc = probe.run(probe.load_manifest(p), tmp_path / "out", log=lambda *_: None)
    assert rc == probe.EXIT_OK
    recs, summ = read_run(tmp_path / "out")
    timed = [r for r in recs if r["phase"] == "repeat"]
    assert len(timed) == 7
    for r in timed:
        for k in ("first_byte", "headers", "total"):
            assert abs(r["t"][k] - 0.300) / 0.300 <= 0.05, (k, r["t"])
        assert r["t"]["connect"] <= r["t"]["first_byte"] <= r["t"]["headers"] <= r["t"]["total"]
    totals = [r["t"]["total"] for r in timed]
    assert (max(totals) - min(totals)) / statistics.median(totals) <= 0.05
    assert summ["cases"]["d300"]["arms"]["A"]["total"]["n"] == 7


def test_body_delay_separates_headers_from_total(server, tmp_path):
    base, _ = server
    # Delays of a few hundred ms: the fake server's own sleep overshoot (a few ms on
    # macOS) stays well inside the 5% criterion of the Stage 0 plan.
    p = manifest(tmp_path, base, cases=[{"id": "b", "query": "delay_ms=300&body_delay_ms=400"}],
                 design={"warmup": 1, "repeats": 3})
    assert probe.run(probe.load_manifest(p), tmp_path / "out", log=lambda *_: None) == probe.EXIT_OK
    recs, _ = read_run(tmp_path / "out")
    for r in recs:
        assert abs(r["t"]["headers"] - 0.300) / 0.300 <= 0.05
        assert abs((r["t"]["total"] - r["t"]["headers"]) - 0.400) / 0.400 <= 0.05


# ------------------------------------------------------------------ records


def test_status_headers_size_sha_and_redaction(server, tmp_path):
    base, _ = server
    p = manifest(tmp_path, base, cases=[{"id": "c", "query": "w=1024&h=512&pad=5000&gwc=HIT"}],
                 design={"warmup": 0, "repeats": 1})
    assert probe.run(probe.load_manifest(p), tmp_path / "out", save_bodies=True, log=lambda *_: None) == 0
    recs, summ = read_run(tmp_path / "out")
    (r,) = recs
    body = fake_png(1024, 512, 5000)
    import hashlib

    assert r["status"] == 200 and r["ok"] is True and r["failure"] is None
    assert r["body_bytes"] == len(body)
    assert r["sha256"] == hashlib.sha256(body).hexdigest()
    assert r["png_size"] == [1024, 512]
    hdr = {k.lower(): v for k, v in r["headers"]}
    assert hdr["geowebcache-cache-result"] == "HIT"
    assert hdr["x-test-cache"] == "BYPASS"
    assert hdr["content-length"] == str(len(body))
    assert hdr["set-cookie"] == probe.REDACTED
    raw = (tmp_path / "out" / "requests.jsonl").read_text() + (tmp_path / "out" / "summary.json").read_text()
    assert SECRET not in raw
    assert "geowebcache-cache-result: HIT" in summ["cases"]["c"]["arms"]["A"]["cache_headers"]
    saved = list((tmp_path / "out" / "bodies").iterdir())
    assert len(saved) == 1 and saved[0].read_bytes() == body


def test_url_is_verbatim_and_never_retried(server, tmp_path):
    base, seen = server
    q = "SERVICE=WMS&BBOX=-180,-90,180,90&TIME=2026-10-01T15:59:59Z&LAYERS=marineheatwave:mhw"
    arms = {"A": {"target": "fake"}, "B": {"target": "fake", "extra_query": "STYLES=alt"}}
    p = manifest(tmp_path, base, cases=[{"id": "v", "query": q}], arms=arms,
                 design={"warmup": 1, "pairs": 2})
    assert probe.run(probe.load_manifest(p), tmp_path / "out", log=lambda *_: None) == 0
    recs, _ = read_run(tmp_path / "out")
    path = "/geoserver/marineheatwave/wms?"
    want = [path + q + ("&STYLES=alt" if r["arm"] == "B" else "") for r in recs]
    assert seen == want  # exact bytes, in order, one request per record, nothing added
    assert len(seen) == len(recs) == 2 + 4


# ------------------------------------------------------------------ pairing and drift


def test_plan_order_is_balanced_and_warmups_first():
    m = probe.Manifest("b", 100, 5, True, {}, {"A": None, "B": None},
                       [probe.Case("c1", "x=1", (200,), None, False), probe.Case("c2", "x=2", (200,), None, False)],
                       warmup=1, pairs=4, repeats=0, a_reference=1)
    plan = probe.build_plan(m)
    c1 = [(s.arm, s.phase, s.pair) for s in plan if s.case == "c1"]
    assert c1 == [("A", "warmup", None), ("B", "warmup", None), ("A", "a_ref_before", None),
                  ("A", "pair", 0), ("B", "pair", 0), ("B", "pair", 1), ("A", "pair", 1),
                  ("B", "pair", 2), ("A", "pair", 2), ("A", "pair", 3), ("B", "pair", 3),
                  ("A", "a_ref_after", None)]
    assert [s.seq for s in plan] == list(range(1, len(plan) + 1))
    assert plan[len(c1)].case == "c2"  # cases run one after another


def test_drift_calculation():
    d = probe.drift([1.0, 1.0, 1.2, 1.2])
    assert d["basis"] == "halves" and d["pct"] == pytest.approx(20.0)
    d = probe.drift([9, 9], before=[2.0, 2.0], after=[1.9, 1.9])
    assert d["basis"] == "a_reference" and d["pct"] == pytest.approx(-5.0)
    assert probe.drift([1.0]) is None


def test_paired_run_differences_and_drift(server, tmp_path):
    base, _ = server
    arms = {"A": {"target": "fake", "extra_query": "delay_ms=100"},
            "B": {"target": "fake", "extra_query": "delay_ms=250"}}
    p = manifest(tmp_path, base, cases=[{"id": "ab", "query": "SERVICE=WMS"}], arms=arms,
                 design={"warmup": 1, "pairs": 4, "a_reference": 1})
    assert probe.run(probe.load_manifest(p), tmp_path / "out", log=lambda *_: None) == 0
    recs, summ = read_run(tmp_path / "out")
    assert [(r["arm"], r["phase"]) for r in recs][:3] == [("A", "warmup"), ("B", "warmup"), ("A", "a_ref_before")]
    s = summ["cases"]["ab"]
    assert [x["order"] for x in s["pairs"]] == ["AB", "BA", "BA", "AB"]
    for x in s["pairs"]:
        assert abs(x["b_minus_a"] - 0.150) <= 0.015
    assert s["paired_b_minus_a"]["n"] == 4
    assert s["a_drift"]["basis"] == "a_reference" and abs(s["a_drift"]["pct"]) < 10
    assert s["a_b_bytes_identical"] is True
    assert summ["status"] == "complete" and summ["sent"] == summ["planned"] == 2 + 1 + 8 + 1


# ------------------------------------------------------------------ cap


def test_cap_is_required(server, tmp_path):
    base, seen = server
    p = manifest(tmp_path, base, cases=[{"id": "c", "query": "x=1"}], cap=None)
    with pytest.raises(probe.ManifestError, match="cap is required"):
        probe.load_manifest(p)
    assert probe.main(["run", str(p), "--out", str(tmp_path / "o")]) == probe.EXIT_MANIFEST
    assert seen == []


def test_over_cap_refuses_to_start(server, tmp_path):
    base, seen = server
    p = manifest(tmp_path, base, cases=[{"id": "c", "query": "x=1"}], design={"warmup": 1, "repeats": 3}, cap=3)
    assert probe.main(["run", str(p), "--out", str(tmp_path / "o")]) == probe.EXIT_MANIFEST
    with pytest.raises(probe.ManifestError, match="exceed the batch cap"):
        probe.run(probe.load_manifest(p), tmp_path / "o2", log=lambda *_: None)
    assert seen == [] and not (tmp_path / "o" / "requests.jsonl").exists()


def test_cap_from_command_line_and_no_builtin_limit(server, tmp_path):
    base, seen = server
    p = manifest(tmp_path, base, cases=[{"id": "c", "query": "pad=10"}], design={"warmup": 0, "repeats": 60}, cap=3)
    assert probe.main(["run", str(p), "--out", str(tmp_path / "o1"), "--max-requests", "59"]) == probe.EXIT_MANIFEST
    assert seen == []
    # a larger batch is allowed when the cap says so: no hidden low limit in the tool
    assert probe.main(["run", str(p), "--out", str(tmp_path / "o2"), "--max-requests", "60"]) == probe.EXIT_OK
    _, summ = read_run(tmp_path / "o2")
    assert summ["sent"] == 60 == len(seen) and summ["max_requests"] == 60


# ------------------------------------------------------------------ failures stop the batch


@pytest.mark.parametrize(
    "query,kind",
    [
        ("delay_ms=2000", "timeout"),
        ("status=500", "http_status"),
        ("mode=drop", "connection"),
        ("mode=partial", "incomplete_body"),
        ("mode=xml", "content_type"),
    ],
)
def test_failure_stops_batch_and_is_recorded(server, tmp_path, query, kind):
    base, seen = server
    cases = [{"id": "good", "query": "pad=10"}, {"id": "bad", "query": query}, {"id": "after", "query": "pad=11&after=1"}]
    p = manifest(tmp_path, base, cases=cases, design={"warmup": 0, "repeats": 2}, timeout=0.5)
    t0 = time.monotonic()
    rc = probe.run(probe.load_manifest(p), tmp_path / "out", log=lambda *_: None)
    assert time.monotonic() - t0 < 5
    assert rc == probe.EXIT_STOPPED
    recs, summ = read_run(tmp_path / "out")
    assert [r["case"] for r in recs] == ["good", "good", "bad"]
    assert [r["ok"] for r in recs] == [True, True, False]
    assert recs[-1]["failure"]["kind"] == kind
    assert summ["status"] == "stopped" and summ["stop"]["kind"] == kind and summ["stop"]["seq"] == 3
    assert summ["sent"] == 3 and summ["failed"] == 1 and summ["not_sent"] == 3
    assert not any("after" in s for s in seen)
    time.sleep(0.05)
    assert len(seen) == 3  # one attempt for the failed request, no retry


def test_connection_refused_is_a_failure(tmp_path):
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    p = manifest(tmp_path, f"http://127.0.0.1:{port}", cases=[{"id": "c", "query": "x=1"}],
                 design={"warmup": 0, "repeats": 2})
    assert probe.run(probe.load_manifest(p), tmp_path / "out", log=lambda *_: None) == probe.EXIT_STOPPED
    recs, summ = read_run(tmp_path / "out")
    assert len(recs) == 1 and recs[0]["failure"]["kind"] == "connection" and recs[0]["sha256"] is None
    assert summ["ok"] == 0 and summ["status"] == "stopped"


# ------------------------------------------------------------------ manifest safety


@pytest.mark.parametrize(
    "patch,msg",
    [
        ({"targets": {"fake": {"base_url": "http://user:pw@127.0.0.1:1/wms"}}}, "credentials in the URL"),
        ({"request_headers": {"Authorization": "Basic eA=="}}, "credentials are not allowed"),
        ({"request_headers": {"Cookie": "a=b"}}, "credentials are not allowed"),
        ({"cases": [{"id": "c", "query": "?x=1"}]}, "raw query string"),
        ({"cases": [{"id": "c", "query": "x=1 2"}]}, "raw query string"),
        ({"arms": {"B": {"target": "fake"}}}, "'arms' must be"),
        ({"max_requests": 0}, "positive integer"),
        ({"schema": "other"}, "schema"),
    ],
)
def test_manifest_rejections(tmp_path, patch, msg):
    d = json.loads(manifest(tmp_path, "http://127.0.0.1:1", cases=[{"id": "c", "query": "x=1"}]).read_text())
    d.update(patch)
    p = tmp_path / "m2.json"
    p.write_text(json.dumps(d))
    with pytest.raises(probe.ManifestError, match=msg):
        probe.load_manifest(p)


def test_nonempty_output_dir_is_refused(server, tmp_path):
    base, seen = server
    out = tmp_path / "out"
    out.mkdir()
    (out / "old.txt").write_text("keep")
    p = manifest(tmp_path, base, cases=[{"id": "c", "query": "x=1"}])
    assert probe.main(["run", str(p), "--out", str(out)]) == probe.EXIT_OUTDIR
    assert seen == [] and (out / "old.txt").read_text() == "keep"


def test_example_manifest_validates_without_network():
    m = probe.load_manifest(API_DIR / "tools" / "p2" / "example_manifest.json")
    assert len(m.plan) <= m.max_requests
    for arm in m.arms.values():
        assert urlsplit(arm.base_url).hostname in ("localhost", "127.0.0.1")


def test_cli_plan_sends_nothing_and_run_works_in_subprocess(server, tmp_path):
    base, seen = server
    p = manifest(tmp_path, base, cases=[{"id": "c", "query": "pad=10"}], design={"warmup": 1, "repeats": 2})
    env = {k: v for k, v in os.environ.items() if not k.lower().endswith("_proxy")}
    r = subprocess.run([sys.executable, "-m", "tools.p2.probe", "plan", str(p)], cwd=API_DIR,
                       capture_output=True, text=True, env=env, timeout=60)
    assert r.returncode == 0 and "3 requests planned, cap 100" in r.stdout and seen == []
    r = subprocess.run([sys.executable, "-m", "tools.p2.probe", "run", str(p), "--out", str(tmp_path / "o")],
                       cwd=API_DIR, capture_output=True, text=True, env=env, timeout=60)
    assert r.returncode == 0, r.stderr
    assert len(seen) == 3
