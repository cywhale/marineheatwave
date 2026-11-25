"""P1-2 V-b: the API under a real gunicorn + uvicorn-worker, isolated (doc 29 r3 §5).

Uses the stage A candidate runner (tools/bench/candidate.run_step: fresh
environment, proven-identity stop, exit / port check) with an explicit argv:

    python -B -m tools.p1_2.vb_server --dir D --port P --python PY --gunicorn GUNICORN \\
        [--worker-class uvicorn_worker.UvicornWorker] [--workers 2]

D/app/API must be (a link to) the API checkout whose data/mhw.zarr is the
read-only copy.  Steps (each a fresh server; a step directory is never reused):

  asgi_echo  tools.p1_2.asgi_echo:app, -w 1: the scope the real worker sends
  golden     mhw_app:app, -w N: every golden case (JSON / CSV / errors byte-
             compared, PNG status + headers + signature + size), CSV Range
             200 / 206 / 416 / 400 checked as a whole, client disconnects in
             the middle of a CSV body, and no temporary file left in D/tmp

No TLS (production's certificate is not used here).  The real clock is used:
golden was recorded with the harness clock fixed at 2026-11-20; cases whose
default dates depend on "today" give the same result within November 2026.
"""
from __future__ import annotations

import argparse
import json
import socket
import struct
import sys
import time
from pathlib import Path

import httpx

from tools.bench import candidate as cand
from tools.bench import load

API = Path(__file__).resolve().parents[2]


def argv_for(app: str, a, D: Path, step: str, workers: int) -> list[str]:
    run_dir = D / "run" / step
    return [a.python, "-B", a.gunicorn, app, "-w", str(workers), "-k", a.worker_class,
            "-b", f"127.0.0.1:{a.port}", "--timeout", "180", "--graceful-timeout", "30",
            "--pid", str(run_dir / "gunicorn.pid"), "--worker-tmp-dir", str(D / "tmp"),
            "--error-logfile", str(D / "logs" / f"gunicorn_{step}.err"),
            "--access-logfile", str(D / "logs" / f"access_{step}.log")]


def scenario_echo(base, sampler, timeout=load.DEFAULT_TIMEOUT):
    with httpx.Client(base_url=base, timeout=timeout) as c:
        r = c.get("/")
    r.raise_for_status()
    return {"scope": r.json(), "server_header": r.headers.get("server")}


def _raw_partial_read(port: int, path: str, nbytes: int) -> dict:
    """Read the first bytes of a response, then close the socket (client gone)."""
    s = socket.create_connection(("127.0.0.1", port), timeout=30)
    s.sendall(f"GET {path} HTTP/1.1\r\nHost: 127.0.0.1\r\nConnection: close\r\n\r\n".encode())
    got = b""
    while len(got) < nbytes:
        chunk = s.recv(nbytes - len(got))
        if not chunk:
            break
        got += chunk
    s.setsockopt(socket.SOL_SOCKET, socket.SO_LINGER, struct.pack("ii", 1, 0))  # RST: abrupt close
    s.close()
    return {"read_bytes": len(got), "status_line": got.split(b"\r\n", 1)[0].decode(errors="replace")}


def make_scenario_golden(D: Path, port: int):
    sys.path.insert(0, str(API))
    from tests.golden.cases import CASES
    from tests.golden.runner import expect_problems, fetch, header_problems, load_manifest

    manifest = load_manifest()
    tmp = D / "tmp"

    def leftovers():
        return sorted(p.name for p in tmp.iterdir())

    def scenario(base, sampler, timeout=load.DEFAULT_TIMEOUT):
        out = {"cases": {}, "range": {}, "disconnect": {}, "problems": []}
        P = out["problems"]
        with httpx.Client(base_url=base, timeout=timeout) as c:
            for case in CASES:
                g = manifest["cases"][case["id"]]
                t = time.perf_counter()
                res = fetch(c, case)
                dt = time.perf_counter() - t
                probs = expect_problems(case, res) + header_problems(case, res)
                if res["status"] != g["status"]:
                    probs.append(f"status {res['status']} != {g['status']}")
                if res["headers"] != g["headers"]:
                    probs.append(f"headers {res['headers']} != {g['headers']}")
                if case["kind"] == "png":
                    if res["body"][:8] != b"\x89PNG\r\n\x1a\n":
                        probs.append("not a PNG")
                    elif struct.unpack(">II", res["body"][16:24]) != (g["summary"]["width"], g["summary"]["height"]):
                        probs.append("PNG size differs")
                    compared = "png-format"
                else:
                    compared = "sha256"
                    if res["sha256"] != g["sha256"]:
                        probs.append(f"body sha256 differs ({res['size']} vs {g['size']} bytes)")
                out["cases"][case["id"]] = {"status": res["status"], "size": res["size"], "seconds": round(dt, 3),
                                            "compared": compared, "ok": not probs}
                P += [f"{case['id']}: {p}" for p in probs]

            c01 = manifest["cases"]["C01"]
            size = c01["size"]
            r = c.get(c01["url"], headers={"Range": "bytes=999999999-"})
            got = {"status": r.status_code, "body": r.content.decode(), "headers": dict(r.headers)}
            got["headers"].pop("date", None)
            got["headers"].pop("server", None)
            want = {"status": 416, "body": "", "headers": {"content-length": "0", "content-range": f"bytes */{size}",
                                                           "content-type": "text/plain; charset=utf-8"}}
            out["range"]["416"] = {"got": got, "ok": got == want}
            if got != want:
                P.append(f"416 differs: {got}")
            r = c.get(c01["url"], headers={"Range": "bytes=0-9"})
            full = c.get(c01["url"]).content
            ok206 = (r.status_code == 206 and r.content == full[:10]
                     and r.headers.get("content-range") == f"bytes 0-9/{size}")
            out["range"]["206"] = {"status": r.status_code, "ok": ok206}
            if not ok206:
                P.append("206 differs")
            r = c.get(c01["url"], headers={"Range": "items=0-1"})
            ok400 = r.status_code == 400 and r.text == "Only support bytes range"
            out["range"]["400"] = {"status": r.status_code, "ok": ok400}
            if not ok400:
                P.append("400 differs")

        # client goes away in the middle of the largest CSV body (1.29 MB)
        c03 = manifest["cases"]["C03"]
        out["disconnect"]["C03"] = [_raw_partial_read(port, c03["url"], 2048) for _ in range(3)]
        time.sleep(2.0)
        with httpx.Client(base_url=base, timeout=timeout) as c:  # still serving afterwards
            out["disconnect"]["after_status"] = c.get(manifest["cases"]["J01"]["url"]).status_code
        time.sleep(1.0)
        out["tmp_leftovers_while_running"] = leftovers()
        if out["tmp_leftovers_while_running"]:
            P.append(f"temporary files left: {out['tmp_leftovers_while_running']}")
        if out["disconnect"]["after_status"] != 200:
            P.append("server not serving after the disconnects")
        if P:
            raise load.LoadError("V-b golden scenario failed", P)
        return out

    return scenario


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--port", type=int, required=True)
    ap.add_argument("--python", required=True)
    ap.add_argument("--gunicorn", required=True)
    ap.add_argument("--worker-class", default="uvicorn_worker.UvicornWorker")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--steps", default="asgi_echo,golden")
    a = ap.parse_args(argv)
    D = Path(a.dir).resolve()
    lim = dict(observe_mb=2500, abort_rss_mb=4000, abort_avail_mb=4096)  # V-b: not a memory gate (V-c is)
    results = {}
    for step in a.steps.split(","):
        if step == "asgi_echo":
            rc = cand.run_step(D, step, a.port, python=a.python, gunicorn=a.gunicorn, expected_workers=1,
                               argv=argv_for("tools.p1_2.asgi_echo:app", a, D, step, 1),
                               scenario_fn=scenario_echo, **lim)
        elif step == "golden":
            rc = cand.run_step(D, step, a.port, python=a.python, gunicorn=a.gunicorn, expected_workers=a.workers,
                               argv=argv_for("mhw_app:app", a, D, step, a.workers),
                               scenario_fn=make_scenario_golden(D, a.port), step_deadline=1800, **lim)
        else:
            raise SystemExit(f"unknown step {step}")
        left = sorted(p.name for p in (D / "tmp").iterdir())
        results[step] = {"exit": rc, "tmp_after_stop": left}
        print(json.dumps({"step": step, "exit": rc, "tmp_after_stop": left}), flush=True)
        if rc != 0 and rc != cand.EXIT_OBSERVE:
            break
    (D / "vb_server_summary.json").write_text(json.dumps(results, indent=1))
    return max((r["exit"] for r in results.values()), default=0)


if __name__ == "__main__":
    sys.exit(main())
