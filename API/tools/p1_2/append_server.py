"""P1-2: the API under a real gunicorn + uvicorn-worker on an appended store COPY.

    python -B -m tools.p1_2.append_server --dir D --port P --python PY --gunicorn G \\
        --points points.json --step after_append2

D/app/API is the API checkout whose data/mhw.zarr points at the writable test
copy (never a production store).  One fresh candidate per call (stage A runner:
readiness on every worker's startup marker + HTTP probe, watchdog, proven-
identity stop, exit / port check), 127.0.0.1 only, POLARS_FORCE_PKG=compat.

Checks: every sample point (tools/p1_2/points.py) as JSON and as CSV returns
exactly one row per appended month with the fixture values; a month plot is a
PNG (status, content-type, signature, size only).  Expectations come from
``updater_check expect-points`` so this runs without numpy.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import struct
import sys
from pathlib import Path

import httpx

from tools.bench import candidate as cand
from tools.bench import load
from tools.p1_2 import points as P

LIMITS = dict(observe_mb=1500, abort_rss_mb=3500, abort_avail_mb=4096)


def png_info(body: bytes) -> dict:
    ok = body[:8] == b"\x89PNG\r\n\x1a\n" and body[12:16] == b"IHDR"
    w, h = struct.unpack(">II", body[16:24]) if ok else (None, None)
    return {"png_signature": body[:8] == b"\x89PNG\r\n\x1a\n", "width": w, "height": h}


def make_scenario(pts: dict):
    months = sorted(next(iter(pts.values())))

    def scenario(base, sampler, timeout=load.DEFAULT_TIMEOUT):
        out: dict = {"months": months, "json": {}, "csv": {}, "plot": None}
        problems = []
        with httpx.Client(base_url=base, timeout=timeout) as c:
            for lat, lon in P.SAMPLES:
                want = pts[f"{lat},{lon}"]
                r = c.get(P.url("json", lat, lon, months[0], months[-1]))
                res = {"status": r.status_code}
                if r.status_code == 200:
                    res.update(P.compare(r.json(), want))
                out["json"][f"{lat},{lon}"] = res
                if not (r.status_code == 200 and res["match"]):
                    problems.append(f"json {lat},{lon}: {res}")
                r = c.get(P.url("csv", lat, lon, months[0], months[-1]))
                res = {"status": r.status_code, "content_type": r.headers.get("content-type")}
                if r.status_code == 200:
                    res.update(P.compare(list(csv.DictReader(io.StringIO(r.text))), want, default_tol=1e-5))
                out["csv"][f"{lat},{lon}"] = res
                if not (r.status_code == 200 and res["match"] and "text/csv" in (res["content_type"] or "")):
                    problems.append(f"csv {lat},{lon}: {res}")
            r = c.get(f"/api/mhw/plot?bbox=(150,-15,152,-13)&period={months[-1][:4]}&sstvar=sst&mode=month")
            out["plot"] = {"status": r.status_code, "content_type": r.headers.get("content-type"),
                           "bytes": len(r.content), **png_info(r.content)}
            pl = out["plot"]
            if not (pl["status"] == 200 and pl["content_type"] == "image/png" and pl["png_signature"]
                    and (pl["width"] or 0) > 0 and (pl["height"] or 0) > 0):
                problems.append(f"plot: {pl}")
        if problems:
            raise load.LoadError("appended-store API check failed", problems)
        return out
    return scenario


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True)
    ap.add_argument("--port", type=int, required=True)
    ap.add_argument("--python", required=True)
    ap.add_argument("--gunicorn", required=True)
    ap.add_argument("--points", required=True)
    ap.add_argument("--step", required=True)
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args(argv)
    pts = json.loads(Path(a.points).read_text())
    D = Path(a.dir).resolve()
    rc = cand.run_step(D, a.step, a.port, python=a.python, gunicorn=a.gunicorn, expected_workers=a.workers,
                       workers=a.workers, extra_env={"POLARS_FORCE_PKG": "compat"}, scenario_fn=make_scenario(pts),
                       scenario_name=f"append_api_{a.step}", step_deadline=900, **LIMITS)
    return rc


if __name__ == "__main__":
    sys.exit(main())
