"""P1-2 V-d: data acceptance of the switched production API (L2 loopback).

    python -B -m tools.p1_2.vd_accept --base https://127.0.0.1:8030 --out accept.json

Narrow, sequential queries only (owner 2026-11-25; the case set of doc 27 §4.1):
* the 17 golden data cases of doc 27 (JSON / CSV / errors): status, compared
  headers, body sha256 equal to tests/golden/manifest.json (values, dates,
  order, NaN/null are all in those bytes);
* J08 CSV (crosses 180, dateline behaviour unchanged): body sha256 equal to the
  57071aa output recorded in doc 27 (2e605c48...);
* P01 plot: status, content-type, PNG signature and image size only (PNG bytes
  and drawing are not compared, owner decision);
* C01 Range: 416 with ``content-range: bytes */<size>`` (accepted change),
  206 for bytes=0-9, 400 for a non-bytes unit.

Exit 0 only when every check passes; the JSON (checks / failed_checks / ok, as
updater_check.with_checks) is always written first.  TLS is not verified
(loopback to production's own certificate).
"""
from __future__ import annotations

import argparse
import json
import struct
import sys
import time
from pathlib import Path

import httpx

API = Path(__file__).resolve().parents[2]
DOC27_CASES = ("J01", "J02", "J04", "J04b", "J05b", "J06d", "J07", "J12", "E01", "E02", "E02b", "E03",
               "C01", "C02", "C03", "C04", "C05")
J08_CSV = ("/api/mhw/csv?lon0=170&lon1=-170&lat0=-5&lat1=5&start=2020-01-01&end=2021-12-01&append=sst",
           "2e605c481b2d8c4642df49de0c594c5ea52c766189abb90326d35544ee1c9355")
PNG_SIG = b"\x89PNG\r\n\x1a\n"


def png_size(body: bytes):
    if body[:8] != PNG_SIG or body[12:16] != b"IHDR":
        return None
    return list(struct.unpack(">II", body[16:24]))


def range_checks(c01: dict, r416, r206, full: bytes, r400) -> dict:
    """Pure: the three C01 Range answers against the golden C01 size."""
    size = c01["size"]
    return {
        "range_416": (r416["status"] == 416 and r416["body"] == b""
                      and r416["headers"].get("content-range") == f"bytes */{size}"),
        "range_206": (r206["status"] == 206 and r206["body"] == full[:10]
                      and r206["headers"].get("content-range") == f"bytes 0-9/{size}"),
        "range_400": r400["status"] == 400 and r400["body"] == b"Only support bytes range",
    }


def run(base: str, cases: list[str], timeout: float = 120.0) -> dict:
    sys.path.insert(0, str(API))
    from tests.golden.cases import CASES
    from tests.golden.runner import expect_problems, fetch, header_problems, load_manifest
    from tools.p1_2.updater_check import with_checks

    manifest = load_manifest()
    by_id = {c["id"]: c for c in CASES}
    out: dict = {"base": base, "cases": {}, "j08_csv": None, "p01": None, "range": {}, "started": time.time()}
    checks: dict = {}

    def raw(r):
        h = {k.lower(): v for k, v in r.headers.items()}
        return {"status": r.status_code, "body": r.content, "headers": h}

    with httpx.Client(base_url=base, timeout=timeout, verify=False) as c:
        for cid in cases:
            case, g = by_id[cid], manifest["cases"][cid]
            t = time.perf_counter()
            res = fetch(c, case)
            dt = time.perf_counter() - t
            probs = expect_problems(case, res) + header_problems(case, res)
            if res["status"] != g["status"]:
                probs.append(f"status {res['status']} != {g['status']}")
            if res["headers"] != g["headers"]:
                probs.append(f"headers {res['headers']} != {g['headers']}")
            if res["sha256"] != g["sha256"]:
                probs.append(f"body sha256 differs ({res['size']} vs {g['size']} bytes)")
            out["cases"][cid] = {"status": res["status"], "size": res["size"], "sha256": res["sha256"],
                                 "seconds": round(dt, 3), "problems": probs}
            checks[f"case_{cid}"] = not probs

        t = time.perf_counter()
        r = c.get(J08_CSV[0])
        out["j08_csv"] = {"status": r.status_code, "size": len(r.content), "content_type": r.headers.get("content-type"),
                          "sha256": __import__("hashlib").sha256(r.content).hexdigest(),
                          "seconds": round(time.perf_counter() - t, 3)}
        checks["j08_csv_equals_57071aa"] = r.status_code == 200 and out["j08_csv"]["sha256"] == J08_CSV[1]

        p01 = manifest["cases"]["P01"]
        r = c.get(by_id["P01"]["url"])
        out["p01"] = {"status": r.status_code, "content_type": r.headers.get("content-type"), "bytes": len(r.content),
                      "png_size": png_size(r.content), "golden_size": [p01["summary"]["width"], p01["summary"]["height"]]}
        checks["p01_png"] = (r.status_code == 200 and r.headers.get("content-type") == "image/png"
                             and out["p01"]["png_size"] == out["p01"]["golden_size"])

        c01 = manifest["cases"]["C01"]
        full = c.get(c01["url"]).content
        r416 = raw(c.get(c01["url"], headers={"Range": "bytes=999999999-"}))
        r206 = raw(c.get(c01["url"], headers={"Range": "bytes=0-9"}))
        r400 = raw(c.get(c01["url"], headers={"Range": "items=0-1"}))
        rc = range_checks(c01, r416, r206, full, r400)
        out["range"] = {"416": {"status": r416["status"], "content_range": r416["headers"].get("content-range")},
                        "206": {"status": r206["status"], "content_range": r206["headers"].get("content-range")},
                        "400": {"status": r400["status"]}}
        checks.update(rc)
    out["ended"] = time.time()
    return with_checks(out, checks)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cases", default=",".join(DOC27_CASES))
    a = ap.parse_args(argv)
    res = run(a.base, a.cases.split(","))
    Path(a.out).write_text(json.dumps(res, indent=1, default=str))
    print(json.dumps({"ok": res["ok"], "failed_checks": res["failed_checks"],
                      "cases": {k: (v["status"], v["seconds"]) for k, v in res["cases"].items()},
                      "j08_csv": res["j08_csv"], "p01": res["p01"], "range": res["range"]}, indent=1))
    return 0 if res["ok"] is True else 3


if __name__ == "__main__":
    sys.exit(main())
