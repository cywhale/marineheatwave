"""L1/L2 of spec 16 r2 §4 and light-query concurrency on L2 (§5.3-style, read-only).

    .venv/bin/python -m tools.bench.latency_layers layers --rounds 5 --out DIR
    .venv/bin/python -m tools.bench.latency_layers loopback-concurrency --per-level 8 --out DIR

L1 = this machine -> https://eco.odb.ntu.edu.tw (nginx).  Every request carries a
unique ``_nc`` parameter so nginx cannot answer from cache; ``x-api-cache`` is
recorded to prove it.  L2 = ``ssh vm24 curl -sk https://127.0.0.1:8030`` (the
app directly; curl's own time_total excludes the ssh overhead).  Rounds
alternate L1/L2 order and rotate the query order.  Raw values are kept.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time
import uuid
from pathlib import Path

PUBLIC = "https://eco.odb.ntu.edu.tw"
LOOPBACK = "https://127.0.0.1:8030"
FULL = "start=1982-01-01&end=2026-10-01"
QUERIES = {
    "Q1_point_level": f"/api/mhw?lon0=135&lat0=15&{FULL}",
    "Q2_point_4vars": f"/api/mhw?lon0=135&lat0=15&{FULL}&append=sst,sst_anomaly,level,td",
    "Q3_point2_level": f"/api/mhw?lon0=135.3&lat0=15.3&{FULL}",
    "Q4_bbox_10x10_10y": "/api/mhw?lon0=130&lon1=140&lat0=10&lat1=20&start=2010-01-01&end=2019-12-01&append=sst_anomaly",
}
CURL_FMT = "%{http_code} %{time_starttransfer} %{time_total} %{size_download}"


def nonce() -> str:
    return "p10" + uuid.uuid4().hex[:10]


def parse(line: str, headers: str) -> dict:
    code, ttfb, total, size = line.split()
    cache = None
    for h in headers.splitlines():
        if h.lower().startswith("x-api-cache:"):
            cache = h.split(":", 1)[1].strip()
    return {"http": int(code), "ttfb": float(ttfb), "total": float(total), "size": int(size), "x_api_cache": cache}


def l1(path: str) -> dict:
    url = f"{PUBLIC}{path}&_nc={nonce()}"
    p = subprocess.run(["curl", "-sS", "-o", "/dev/null", "-D", "-", "-w", "\n" + CURL_FMT, url],
                       capture_output=True, text=True, check=True)
    headers, _, line = p.stdout.rpartition("\n")
    return {"url": url, **parse(line, headers)}


def l2(path: str) -> dict:
    url = f"{LOOPBACK}{path}&_nc={nonce()}"
    remote = f"curl -sk -o /dev/null -D - -w '\\n{CURL_FMT}' '{url}'"
    p = subprocess.run(["ssh", "-o", "BatchMode=yes", "vm24", remote], capture_output=True, text=True, check=True)
    headers, _, line = p.stdout.rpartition("\n")
    return {"url": url, **parse(line, headers)}


def cmd_layers(args) -> None:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    names = list(QUERIES)
    for r in range(args.rounds):
        order = names[r % len(names):] + names[: r % len(names)]
        for q in order:
            layers = [("L1", l1), ("L2", l2)] if r % 2 == 0 else [("L2", l2), ("L1", l1)]
            for layer, fn in layers:
                row = {"round": r, "query": q, "layer": layer, "at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), **fn(QUERIES[q])}
                rows.append(row)
                print(json.dumps(row), flush=True)
    (out / "layers.json").write_text(json.dumps(rows, indent=1))


def cmd_loopback_concurrency(args) -> None:
    """Light point queries only (J01 shape), concurrency 1/2/4, run with xargs on vm24."""
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    results = {}
    for c in (1, 2, 4):
        tag = nonce()
        url = f"{LOOPBACK}{QUERIES['Q1_point_level']}&_nc={tag}_"
        remote = (f"start=$(date +%s.%N); seq {args.per_level} | xargs -P {c} -I{{}} "
                  f"curl -sk -o /dev/null -w '%{{http_code}} %{{time_total}}\\n' '{url}{{}}'; "
                  f"end=$(date +%s.%N); awk -v s=$start -v e=$end 'BEGIN{{printf \"wall %.3f\\n\", e-s}}'")
        p = subprocess.run(["ssh", "-o", "BatchMode=yes", "vm24", remote], capture_output=True, text=True, check=True)
        lines = p.stdout.strip().splitlines()
        wall = float(lines[-1].split()[1])
        reqs = [{"http": int(l.split()[0]), "total": float(l.split()[1])} for l in lines[:-1]]
        results[str(c)] = {"wall": wall, "requests": reqs}
        print(c, json.dumps(results[str(c)]), flush=True)
    (out / "loopback_concurrency.json").write_text(json.dumps(results, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("layers")
    a.add_argument("--rounds", type=int, default=5)
    a.add_argument("--out", required=True)
    b = sub.add_parser("loopback-concurrency")
    b.add_argument("--per-level", type=int, default=8)
    b.add_argument("--out", required=True)
    args = ap.parse_args()
    {"layers": cmd_layers, "loopback-concurrency": cmd_loopback_concurrency}[args.cmd](args)


if __name__ == "__main__":
    main()
