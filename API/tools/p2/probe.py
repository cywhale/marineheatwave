"""Fixed-request GetMap probe with paired A/B comparison (MHW Phase 2, Stage 0).

    python -m tools.p2.probe plan MANIFEST.json            # validate + show the sequence; no network
    python -m tools.p2.probe run  MANIFEST.json --out DIR  # execute one batch

Everything about the target comes from the external manifest (README.md has the
format): base URLs, raw query strings, the batch's request cap and timeout.
The tool itself has no host names, no credentials and no built-in request cap.

Rules the tool enforces:
- The planned request count must not exceed the batch cap (``max_requests``,
  required; ``--max-requests`` may override it).  Over the cap the batch does
  not start; it is never silently truncated.
- URLs are ``base_url`` + the raw query strings from the manifest, verbatim.
  No cache-busting parameter is added; no cache is purged; HTTP proxies from
  the environment are not used (``http.client`` connects directly).
- One attempt per request, no retry.  The first failed request (connection
  error, timeout, unexpected status/content type, invalid PNG, short body)
  stops the batch; the failure is recorded as such and the exit code is 2.
- Records keep what interpretation needs (timings, status, all response
  headers, body size, SHA-256) and redact cookies and authorization headers.
  Manifests with credentials (URL userinfo, auth/cookie request headers) are
  rejected.

Timings (seconds, ``time.perf_counter`` relative to the request start):
``connect`` (TCP, and TLS for https), ``first_byte`` (the socket became readable
after the request was sent), ``headers`` (status line and headers parsed) and
``total`` (whole body received).

Exit codes: 0 complete, 2 stopped on a failed request, 3 invalid manifest or over
the cap, 4 output directory not empty.
"""
from __future__ import annotations

import argparse
import hashlib
import http.client
import json
import select
import socket
import ssl
import statistics
import struct
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlsplit

SCHEMA = "mhw-p2-probe/1"
PNG_SIG = b"\x89PNG\r\n\x1a\n"
CHUNK = 64 * 1024
REDACTED = "[redacted]"
# Response headers whose values are never written to the record.
SECRET_RESPONSE_HEADERS = {"set-cookie", "set-cookie2", "authorization", "proxy-authorization", "cookie"}
# Request headers a manifest may not carry (the record must not contain secrets).
FORBIDDEN_REQUEST_HEADERS = {"authorization", "proxy-authorization", "cookie", "x-api-key", "x-auth-token"}

EXIT_OK, EXIT_STOPPED, EXIT_MANIFEST, EXIT_OUTDIR = 0, 2, 3, 4


class ManifestError(ValueError):
    pass


# --------------------------------------------------------------------------- manifest


@dataclass(frozen=True)
class Arm:
    name: str
    target: str
    base_url: str
    extra_query: str


@dataclass(frozen=True)
class Case:
    id: str
    query: str
    expect_status: tuple[int, ...]
    expect_content_type: str | None
    expect_png: bool


@dataclass(frozen=True)
class Step:
    seq: int
    case: str
    arm: str
    phase: str  # warmup | a_ref_before | pair | a_ref_after | repeat
    pair: int | None = None
    position: int | None = None  # 0 or 1 inside a pair


@dataclass
class Manifest:
    batch: str
    max_requests: int
    timeout_s: float
    tls_verify: bool
    headers: dict[str, str]
    arms: dict[str, Arm]
    cases: list[Case]
    warmup: int
    pairs: int
    repeats: int
    a_reference: int
    raw_sha256: str = ""
    plan: list[Step] = field(default_factory=list)


def _req(d: dict, key: str, typ, where: str):
    if key not in d:
        raise ManifestError(f"{where}: '{key}' is required")
    v = d[key]
    if typ is int and isinstance(v, bool):
        raise ManifestError(f"{where}: '{key}' must be an integer")
    if not isinstance(v, typ):
        raise ManifestError(f"{where}: '{key}' must be {getattr(typ, '__name__', typ)}")
    return v


def _check_query(q: str, where: str) -> str:
    if not isinstance(q, str):
        raise ManifestError(f"{where}: query must be a string")
    if any(c.isspace() for c in q) or "#" in q or q.startswith(("?", "&")):
        raise ManifestError(f"{where}: query must be a raw query string without whitespace, '#' or a leading ?/&")
    return q


def _check_base_url(u: str, where: str) -> str:
    p = urlsplit(u)
    if p.scheme not in ("http", "https") or not p.hostname:
        raise ManifestError(f"{where}: base_url must be an absolute http(s) URL")
    if p.username is not None or p.password is not None:
        raise ManifestError(f"{where}: credentials in the URL are not allowed")
    if p.fragment or any(c.isspace() for c in u):
        raise ManifestError(f"{where}: base_url must not contain a fragment or whitespace")
    return u


def load_manifest(path: Path, max_requests_override: int | None = None) -> Manifest:
    raw = Path(path).read_bytes()
    try:
        d = json.loads(raw)
    except json.JSONDecodeError as e:
        raise ManifestError(f"manifest is not valid JSON: {e}") from None
    if not isinstance(d, dict) or d.get("schema") != SCHEMA:
        raise ManifestError(f"manifest 'schema' must be '{SCHEMA}'")
    batch = _req(d, "batch", str, "manifest")
    cap = max_requests_override if max_requests_override is not None else d.get("max_requests")
    if cap is None:
        raise ManifestError("the batch request cap is required: set 'max_requests' or pass --max-requests")
    if not isinstance(cap, int) or isinstance(cap, bool) or cap < 1:
        raise ManifestError("the batch request cap must be a positive integer")
    timeout = _req(d, "timeout_s", (int, float), "manifest")
    if isinstance(timeout, bool) or timeout <= 0:
        raise ManifestError("'timeout_s' must be a positive number")
    tls_verify = d.get("tls_verify", True)
    if not isinstance(tls_verify, bool):
        raise ManifestError("'tls_verify' must be true or false")

    headers = d.get("request_headers", {}) or {}
    if not isinstance(headers, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in headers.items()):
        raise ManifestError("'request_headers' must map strings to strings")
    bad = sorted(k for k in headers if k.lower() in FORBIDDEN_REQUEST_HEADERS)
    if bad:
        raise ManifestError(f"request headers carrying credentials are not allowed: {bad}")

    targets = _req(d, "targets", dict, "manifest")
    if not targets:
        raise ManifestError("'targets' must not be empty")
    base = {}
    for name, t in targets.items():
        where = f"targets.{name}"
        if not isinstance(t, dict):
            raise ManifestError(f"{where} must be an object")
        base[name] = _check_base_url(_req(t, "base_url", str, where), where)

    arms_d = _req(d, "arms", dict, "manifest")
    if set(arms_d) not in ({"A"}, {"A", "B"}):
        raise ManifestError("'arms' must be {'A'} (series) or {'A', 'B'} (paired)")
    arms = {}
    for name in sorted(arms_d):
        a = arms_d[name]
        where = f"arms.{name}"
        if not isinstance(a, dict):
            raise ManifestError(f"{where} must be an object")
        tgt = _req(a, "target", str, where)
        if tgt not in base:
            raise ManifestError(f"{where}: unknown target '{tgt}'")
        extra = a.get("extra_query", "")
        if extra:
            _check_query(extra, where)
        arms[name] = Arm(name, tgt, base[tgt], extra)

    cases_l = _req(d, "cases", list, "manifest")
    if not cases_l:
        raise ManifestError("'cases' must not be empty")
    cases, seen = [], set()
    for i, c in enumerate(cases_l):
        where = f"cases[{i}]"
        if not isinstance(c, dict):
            raise ManifestError(f"{where} must be an object")
        cid = _req(c, "id", str, where)
        if cid in seen:
            raise ManifestError(f"{where}: duplicate id '{cid}'")
        seen.add(cid)
        q = _check_query(_req(c, "query", str, where), where)
        exp = c.get("expect", {}) or {}
        status = exp.get("status", [200])
        if not isinstance(status, list) or not status or not all(isinstance(s, int) for s in status):
            raise ManifestError(f"{where}: expect.status must be a non-empty list of integers")
        ctype = exp.get("content_type", "image/png")
        png = exp.get("png", ctype == "image/png")
        cases.append(Case(cid, q, tuple(status), ctype, bool(png)))

    design = _req(d, "design", dict, "manifest")
    warmup = design.get("warmup", 1)
    a_ref = design.get("a_reference", 0)
    paired = "B" in arms
    pairs = design.get("pairs", 0) if paired else 0
    repeats = design.get("repeats", 0) if not paired else 0
    for k, v in (("warmup", warmup), ("a_reference", a_ref), ("pairs", pairs), ("repeats", repeats)):
        if not isinstance(v, int) or isinstance(v, bool) or v < 0:
            raise ManifestError(f"design.{k} must be a non-negative integer")
    if paired and pairs < 1:
        raise ManifestError("design.pairs must be >= 1 when arm B is present")
    if not paired and repeats < 1:
        raise ManifestError("design.repeats must be >= 1 for a single-arm (series) batch")
    if not paired and a_ref:
        raise ManifestError("design.a_reference applies only to paired batches")

    m = Manifest(batch, cap, float(timeout), tls_verify, dict(headers), arms, cases,
                 warmup, pairs, repeats, a_ref, hashlib.sha256(raw).hexdigest())
    m.plan = build_plan(m)
    return m


def build_plan(m: Manifest) -> list[Step]:
    """Per case: warm-up (A then B, excluded from timing), A-only references before,
    pairs in balanced order (AB, BA, BA, AB, ...: pair i is AB when i % 4 in {0, 3}),
    A-only references after.  Single-arm batches: warm-up then repeats."""
    steps: list[Step] = []

    def add(case, arm, phase, pair=None, pos=None):
        steps.append(Step(len(steps) + 1, case, arm, phase, pair, pos))

    for c in m.cases:
        for _ in range(m.warmup):
            for arm in m.arms:
                add(c.id, arm, "warmup")
        if "B" in m.arms:
            for _ in range(m.a_reference):
                add(c.id, "A", "a_ref_before")
            for i in range(m.pairs):
                order = ("A", "B") if i % 4 in (0, 3) else ("B", "A")
                for pos, arm in enumerate(order):
                    add(c.id, arm, "pair", i, pos)
            for _ in range(m.a_reference):
                add(c.id, "A", "a_ref_after")
        else:
            for _ in range(m.repeats):
                add(c.id, "A", "repeat")
    return steps


def build_url(arm: Arm, case: Case) -> str:
    sep = "&" if "?" in arm.base_url else "?"
    url = arm.base_url + sep + case.query
    if arm.extra_query:
        url += "&" + arm.extra_query
    return url


# --------------------------------------------------------------------------- one request


def redact_headers(items: list[tuple[str, str]]) -> list[list[str]]:
    return [[k, REDACTED if k.lower() in SECRET_RESPONSE_HEADERS else v] for k, v in items]


def png_size(body: bytes) -> tuple[int, int] | None:
    if len(body) >= 24 and body[:8] == PNG_SIG and body[12:16] == b"IHDR":
        return struct.unpack(">II", body[16:24])
    return None


def fetch(url: str, *, timeout_s: float, headers: dict[str, str], tls_verify: bool) -> dict:
    """One GET, one attempt.  Returns a record; ``error`` is set when the request
    could not complete (status/content checks are applied by the caller)."""
    p = urlsplit(url)
    host, port = p.hostname, p.port
    path = (p.path or "/") + (f"?{p.query}" if p.query else "")
    rec: dict = {"url": url, "started_at": datetime.now(timezone.utc).isoformat(timespec="milliseconds"),
                 "t": {}, "status": None, "reason": None, "headers": [], "body_bytes": 0,
                 "sha256": None, "error": None}
    deadline = None
    conn = None
    t0 = time.perf_counter()

    def left() -> float:
        r = deadline - time.perf_counter()
        if r <= 0:
            raise TimeoutError("batch timeout_s exceeded")
        return r

    try:
        deadline = t0 + timeout_s
        if p.scheme == "https":
            ctx = ssl.create_default_context()
            if not tls_verify:
                ctx.check_hostname = False
                ctx.verify_mode = ssl.CERT_NONE
            conn = http.client.HTTPSConnection(host, port, timeout=timeout_s, context=ctx)
        else:
            conn = http.client.HTTPConnection(host, port, timeout=timeout_s)
        conn.connect()
        rec["t"]["connect"] = time.perf_counter() - t0
        conn.putrequest("GET", path, skip_accept_encoding=True)
        for k, v in headers.items():
            conn.putheader(k, v)
        conn.endheaders()
        sock = conn.sock
        pending = sock.pending() if isinstance(sock, ssl.SSLSocket) else 0
        if not pending:
            ready, _, _ = select.select([sock], [], [], left())
            if not ready:
                raise TimeoutError("no response byte before timeout_s")
        rec["t"]["first_byte"] = time.perf_counter() - t0
        sock.settimeout(left())
        resp = conn.getresponse()
        rec["t"]["headers"] = time.perf_counter() - t0
        rec["status"], rec["reason"] = resp.status, resp.reason
        rec["headers"] = redact_headers(resp.getheaders())
        h = hashlib.sha256()
        body = bytearray()
        # With "Connection: close" http.client hands the socket to the response
        # and closes its own reference, so per-read timeouts go to the socket
        # behind the response's file object (the deadline is checked after each
        # chunk either way).
        body_sock = getattr(getattr(resp.fp, "raw", None), "_sock", None)
        while not resp.isclosed():  # http.client closes the response once Content-Length is read
            if body_sock is not None:
                body_sock.settimeout(left())
            chunk = resp.read(CHUNK)
            if not chunk:
                break
            h.update(chunk)
            body.extend(chunk)
            left()
        rec["t"]["total"] = time.perf_counter() - t0
        rec["body_bytes"] = len(body)
        rec["sha256"] = h.hexdigest()
        rec["_body"] = bytes(body)
        # http.client's read(amt) returns b"" at EOF even when Content-Length was
        # not reached, so a short body has to be detected here.
        clen = resp.getheader("Content-Length")
        if clen is not None and resp.getheader("Transfer-Encoding") is None:
            try:
                expected = int(clen)
            except ValueError:
                expected = None
            if expected is not None and expected != len(body):
                rec["error"] = {"kind": "incomplete_body", "detail": f"got {len(body)} of {expected} bytes"}
    except (TimeoutError, socket.timeout) as e:
        rec["error"] = {"kind": "timeout", "detail": str(e) or "socket timeout"}
    except http.client.IncompleteRead as e:
        rec["error"] = {"kind": "incomplete_body", "detail": f"got {len(e.partial)} bytes, {e.expected} more expected"}
    except (ConnectionError, http.client.RemoteDisconnected, http.client.BadStatusLine) as e:
        rec["error"] = {"kind": "connection", "detail": f"{type(e).__name__}: {e}"}
    except ssl.SSLError as e:
        rec["error"] = {"kind": "tls", "detail": f"{type(e).__name__}: {e}"}
    except OSError as e:
        rec["error"] = {"kind": "connection", "detail": f"{type(e).__name__}: {e}"}
    finally:
        rec["t"]["elapsed"] = time.perf_counter() - t0
        if conn is not None:
            conn.close()
    return rec


def check(rec: dict, case: Case) -> dict | None:
    """Status/content checks.  Returns the failure, or None."""
    if rec["error"]:
        return rec["error"]
    if rec["status"] not in case.expect_status:
        return {"kind": "http_status", "detail": f"status {rec['status']} not in {list(case.expect_status)}"}
    if case.expect_content_type:
        ct = next((v for k, v in rec["headers"] if k.lower() == "content-type"), "")
        if not ct.lower().startswith(case.expect_content_type.lower()):
            return {"kind": "content_type", "detail": f"content-type '{ct}' is not {case.expect_content_type}"}
    if case.expect_png and png_size(rec.get("_body", b"")) is None:
        return {"kind": "not_png", "detail": "body is not a PNG with an IHDR chunk"}
    return None


# --------------------------------------------------------------------------- summary


def _stats(xs: list[float]) -> dict | None:
    if not xs:
        return None
    return {"n": len(xs), "median": statistics.median(xs), "min": min(xs), "max": max(xs)}


def drift(a_times_in_order: list[float], before: list[float] | None = None, after: list[float] | None = None) -> dict | None:
    """Relative change of A over the batch: (median(late) - median(early)) / median(early).
    Uses the A-only reference blocks when present, else the first and second half
    of A's timed observations in request order."""
    if before and after:
        early, late, basis = before, after, "a_reference"
    else:
        if len(a_times_in_order) < 2:
            return None
        h = len(a_times_in_order) // 2
        early, late, basis = a_times_in_order[:h], a_times_in_order[-h:], "halves"
    e, l_ = statistics.median(early), statistics.median(late)
    return {"basis": basis, "early_median": e, "late_median": l_, "pct": (l_ - e) / e * 100 if e else None}


def summarize(m: Manifest, records: list[dict]) -> dict:
    out = {}
    ok = [r for r in records if r["ok"]]
    for c in m.cases:
        rs = [r for r in ok if r["case"] == c.id]
        timed = [r for r in rs if r["phase"] != "warmup"]
        per_arm = {}
        for arm in m.arms:
            ar = [r for r in timed if r["arm"] == arm and r["phase"] in ("pair", "repeat")]
            per_arm[arm] = {
                "total": _stats([r["t"]["total"] for r in ar]),
                "headers": _stats([r["t"]["headers"] for r in ar]),
                "first_byte": _stats([r["t"]["first_byte"] for r in ar]),
                "sha256": sorted({r["sha256"] for r in rs if r["arm"] == arm}),
                "body_bytes": sorted({r["body_bytes"] for r in rs if r["arm"] == arm}),
                "cache_headers": sorted({f"{k}: {v}" for r in rs if r["arm"] == arm
                                         for k, v in r["headers"] if "cache" in k.lower() or k.lower() == "age"}),
            }
        entry = {"arms": per_arm}
        if "B" in m.arms:
            pairs = {}
            for r in timed:
                if r["phase"] == "pair":
                    pairs.setdefault(r["pair"], {})[r["arm"]] = r["t"]["total"]
            diffs = [{"pair": i, "order": "AB" if i % 4 in (0, 3) else "BA", "A": p["A"], "B": p["B"],
                      "b_minus_a": p["B"] - p["A"], "b_over_a": p["B"] / p["A"] if p["A"] else None}
                     for i, p in sorted(pairs.items()) if "A" in p and "B" in p]
            d = [x["b_minus_a"] for x in diffs]
            entry["pairs"] = diffs
            entry["paired_b_minus_a"] = _stats(d)
            a_order = [r["t"]["total"] for r in timed if r["arm"] == "A" and r["phase"] == "pair"]
            entry["a_drift"] = drift(
                a_order,
                [r["t"]["total"] for r in timed if r["phase"] == "a_ref_before"],
                [r["t"]["total"] for r in timed if r["phase"] == "a_ref_after"],
            )
            entry["a_b_bytes_identical"] = bool(per_arm["A"]["sha256"]) and per_arm["A"]["sha256"] == per_arm["B"]["sha256"]
        else:
            entry["a_drift"] = drift([r["t"]["total"] for r in timed if r["phase"] == "repeat"])
        out[c.id] = entry
    return out


# --------------------------------------------------------------------------- run


def run(m: Manifest, out: Path, *, save_bodies: bool = False, log=print) -> int:
    planned = len(m.plan)
    if planned > m.max_requests:
        raise ManifestError(f"planned {planned} requests exceed the batch cap {m.max_requests}; nothing was sent")
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        return EXIT_OUTDIR
    out.mkdir(parents=True, exist_ok=True)
    if save_bodies:
        (out / "bodies").mkdir()
    cases = {c.id: c for c in m.cases}
    records: list[dict] = []
    sent = 0
    stop = None
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with open(out / "requests.jsonl", "w", encoding="utf-8") as fh:
        for step in m.plan:
            if sent >= m.max_requests:  # defensive; the plan was checked above
                stop = {"seq": step.seq, "kind": "cap", "detail": "batch cap reached"}
                break
            arm, case = m.arms[step.arm], cases[step.case]
            rec = fetch(build_url(arm, case), timeout_s=m.timeout_s, headers=m.headers, tls_verify=m.tls_verify)
            sent += 1
            fail = check(rec, case)
            body = rec.pop("_body", None)
            if save_bodies and body is not None:
                (out / "bodies" / f"{step.seq:05d}_{step.case}_{step.arm}.bin").write_bytes(body)
            size = png_size(body or b"")
            rec.update({"seq": step.seq, "case": step.case, "arm": step.arm, "target": arm.target,
                        "phase": step.phase, "pair": step.pair, "position": step.position,
                        "png_size": list(size) if size else None, "ok": fail is None, "failure": fail})
            fh.write(json.dumps(rec, sort_keys=True) + "\n")
            fh.flush()
            records.append(rec)
            log(f"{step.seq:4d}/{planned} {step.case} {step.arm} {step.phase:12s} "
                f"{'ok' if fail is None else 'FAIL ' + fail['kind']} total={rec['t'].get('total')}")
            if fail is not None:
                stop = {"seq": step.seq, **fail}
                break
    summary = {
        "schema": SCHEMA, "batch": m.batch, "manifest_sha256": m.raw_sha256,
        "started_at": started, "finished_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status": "complete" if stop is None else "stopped", "stop": stop,
        "max_requests": m.max_requests, "planned": planned, "sent": sent,
        "ok": sum(r["ok"] for r in records), "failed": sum(not r["ok"] for r in records),
        "not_sent": planned - sent, "design": {"warmup": m.warmup, "pairs": m.pairs, "repeats": m.repeats,
                                               "a_reference": m.a_reference, "arms": sorted(m.arms)},
        "timer": "time.perf_counter; seconds from request start",
        "cases": summarize(m, records),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return EXIT_OK if stop is None else EXIT_STOPPED


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m tools.p2.probe", description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("plan", "run"):
        s = sub.add_parser(name)
        s.add_argument("manifest", type=Path)
        s.add_argument("--max-requests", type=int, default=None,
                       help="batch request cap; overrides the manifest's max_requests")
        if name == "run":
            s.add_argument("--out", type=Path, required=True, help="new or empty output directory")
            s.add_argument("--save-bodies", action="store_true", help="also keep every response body")
    a = ap.parse_args(argv)
    try:
        m = load_manifest(a.manifest, a.max_requests)
        if len(m.plan) > m.max_requests:
            raise ManifestError(f"planned {len(m.plan)} requests exceed the batch cap {m.max_requests}; nothing was sent")
    except ManifestError as e:
        print(f"manifest error: {e}", file=sys.stderr)
        return EXIT_MANIFEST
    if a.cmd == "plan":
        print(f"batch {m.batch}: {len(m.plan)} requests planned, cap {m.max_requests}, timeout {m.timeout_s} s")
        for st in m.plan:
            arm = m.arms[st.arm]
            case = next(c for c in m.cases if c.id == st.case)
            print(f"{st.seq:4d} {st.case} {st.arm} {st.phase} pair={st.pair} {build_url(arm, case)}")
        return EXIT_OK
    rc = run(m, a.out, save_bodies=a.save_bodies)
    if rc == EXIT_OUTDIR:
        print(f"output directory {a.out} is not empty", file=sys.stderr)
    return rc


if __name__ == "__main__":
    sys.exit(main())
