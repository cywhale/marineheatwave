"""Record and inspect golden outputs (spec 16 r2 §3).

Run from the API directory with the project venv:

    .venv/bin/python -m tests.golden.runner record --tag 0.2.8 --runs 3
    .venv/bin/python -m tests.golden.runner run-once --out /some/dir [--only J01,P01]
    .venv/bin/python -m tests.golden.runner amend-headers --runs 3

``amend-headers`` adds newly compared headers to an existing manifest.  It
re-runs every case in fresh processes, refuses to proceed unless every body
hash and status equals the manifest, and never rewrites golden bodies.

``record`` starts one fresh Python process per run, so every run re-imports the
app and re-enters its lifespan.  Only cases whose bytes are identical across
all runs become byte baselines; the others are written to the manifest with
``stable: false`` and must be handled per spec 16 §3.3 before use.

Bodies are stored outside the repo (``$MHW_GOLDEN_DIR``, default
``~/backup/mhw/golden``); the repo keeps only ``tests/golden/manifest.json``.
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import io
import json
import os
import shutil
import struct
import subprocess
import sys
from pathlib import Path

API_DIR = Path(__file__).resolve().parents[2]
MANIFEST = Path(__file__).with_name("manifest.json")
GOLDEN_ROOT = Path(os.environ.get("MHW_GOLDEN_DIR", "~/backup/mhw/golden")).expanduser()
EXT = {"json": "json", "csv": "csv", "png": "png"}
BASE_HEADERS = ("content-type", "content-disposition")
# FileResponse headers for CSV (Codex round 4 #2).  etag, last-modified and date
# change on every response and stay excluded.
CSV_HEADERS = BASE_HEADERS + ("accept-ranges", "content-length")


def compared_headers(case: dict, status: int) -> tuple:
    return CSV_HEADERS if case["kind"] == "csv" and status == 200 else BASE_HEADERS


class ManifestError(RuntimeError):
    pass


def load_manifest(path: Path = MANIFEST) -> dict:
    if not path.exists():
        raise ManifestError(f"golden manifest missing: {path}")
    return json.loads(path.read_text())


def validate_manifest(manifest: dict, cases: list[dict]) -> list[str]:
    """Problems that must fail a regression run: case/URL drift and unapproved instability."""
    problems = []
    recorded = manifest.get("cases", {})
    approved = set(manifest.get("approved_unstable", []))
    ids = [c["id"] for c in cases]
    for cid in ids:
        if cid not in recorded:
            problems.append(f"{cid}: in CASES but missing from manifest")
    for cid in recorded:
        if cid not in ids:
            problems.append(f"{cid}: in manifest but not in CASES")
    for c in cases:
        g = recorded.get(c["id"])
        if g is None:
            continue
        if g.get("url") != c["url"]:
            problems.append(f"{c['id']}: URL differs (manifest {g.get('url')!r} vs CASES {c['url']!r})")
        if g.get("kind") != c["kind"]:
            problems.append(f"{c['id']}: kind differs ({g.get('kind')} vs {c['kind']})")
        if not g.get("stable", False) and c["id"] not in approved:
            problems.append(f"{c['id']}: recorded as unstable and not in approved_unstable")
        for h in compared_headers(c, g.get("status", 0)):
            if h not in g.get("headers", {}):
                problems.append(f"{c['id']}: manifest lacks header evidence for {h!r}")
    return problems


def stability_problems(per_run: list[dict]) -> list[str]:
    """Compare status, compared headers and body hash of one case across runs."""
    problems = []
    first = per_run[0]
    for i, r in enumerate(per_run[1:], start=2):
        if r["status"] != first["status"]:
            problems.append(f"run{i} status {r['status']} != run1 {first['status']}")
        if r["headers"] != first["headers"]:
            problems.append(f"run{i} headers {r['headers']} != run1 {first['headers']}")
        if r["sha256"] != first["sha256"]:
            problems.append(f"run{i} body sha256 differs from run1")
    return problems


def header_problems(case: dict, res: dict) -> list[str]:
    """Internal consistency of a single response (content-length matches the body)."""
    problems = []
    cl = res["headers"].get("content-length")
    if case["kind"] == "csv" and res["status"] == 200:
        if cl is None:
            problems.append("CSV response has no content-length")
        elif int(cl) != res["size"]:
            problems.append(f"content-length {cl} != body length {res['size']}")
        if res["headers"].get("accept-ranges") is None:
            problems.append("CSV response has no accept-ranges")
    return problems


def body_path(directory: Path, case: dict, status: int) -> Path:
    ext = EXT[case["kind"]] if status == 200 else "json"
    return directory / f"{case['id']}.{ext}"


def summarize(case: dict, status: int, headers: dict, body: bytes) -> dict:
    """Small, human-readable digest used to locate differences quickly."""
    ctype = headers.get("content-type", "")
    if ctype.startswith("application/json"):
        data = json.loads(body)
        if isinstance(data, list):
            columns = list(data[0].keys()) if data else []
            nulls = {c: sum(1 for r in data if r.get(c) is None) for c in columns}
            dates = sorted({r["date"] for r in data if "date" in r})
            return {"rows": len(data), "columns": columns, "nulls": nulls,
                    "months": len(dates), "date_range": [dates[0], dates[-1]] if dates else None,
                    "first": data[:3], "last": data[-3:]}
        return {"object": data}
    if ctype.startswith("text/csv"):
        rows = list(csv.reader(io.StringIO(body.decode())))
        return {"rows": len(rows) - 1, "header": rows[0] if rows else [],
                "first": rows[1:3], "last": rows[-2:]}
    if ctype.startswith("image/png"):
        width, height = struct.unpack(">II", body[16:24])
        return {"width": width, "height": height}
    return {}


def check_expect(case: dict, status: int, summary: dict) -> list[str]:
    """Check the branch each case is meant to hit (spec 16 §3.1 '預期分支')."""
    exp = case.get("expect", {})
    problems = []
    if status != exp.get("status"):
        problems.append(f"status {status} != {exp.get('status')}")
        return problems
    if "rows" in exp and summary.get("rows") != exp["rows"]:
        problems.append(f"rows {summary.get('rows')} != {exp['rows']}")
    if "months" in exp and summary.get("months") != exp["months"]:
        problems.append(f"months {summary.get('months')} != {exp['months']}")
    if "columns" in exp and summary.get("columns") != exp["columns"]:
        problems.append(f"columns {summary.get('columns')} != {exp['columns']}")
    if "first" in exp:
        first = (summary.get("first") or [{}])[0]
        for key, value in exp["first"].items():
            if first.get(key) != value:
                problems.append(f"first.{key} {first.get(key)} != {value}")
    detail = (summary.get("object") or {}).get("detail")
    if "detail" in exp and detail != exp["detail"]:
        problems.append(f"detail {detail!r} != {exp['detail']!r}")
    if "detail_prefix" in exp and not str(detail).startswith(exp["detail_prefix"]):
        problems.append(f"detail {detail!r} does not start with {exp['detail_prefix']!r}")
    return problems


def fetch(client, case: dict) -> dict:
    resp = client.get(case["url"])
    headers = {k: resp.headers.get(k) for k in compared_headers(case, resp.status_code)}
    body = resp.content
    summary = summarize(case, resp.status_code, headers, body)
    if case["kind"] == "json" and "all_level" in case.get("expect", {}) and resp.status_code == 200:
        levels = {r.get("level") for r in json.loads(body)}
        summary["level_values"] = sorted(levels)
    return {"status": resp.status_code, "headers": headers, "body": body,
            "sha256": hashlib.sha256(body).hexdigest(), "size": len(body), "summary": summary}


def expect_problems(case: dict, result: dict) -> list[str]:
    problems = check_expect(case, result["status"], result["summary"])
    want = case.get("expect", {}).get("all_level")
    if want is not None and result["summary"].get("level_values") != [want]:
        problems.append(f"level values {result['summary'].get('level_values')} != [{want}]")
    return problems


def make_client():
    from tests.harness import isolate_env, patch_app

    isolate_env()
    app = patch_app()
    from fastapi.testclient import TestClient

    return TestClient(app)


def git_commit() -> str:
    out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=API_DIR, capture_output=True, text=True)
    dirty = subprocess.run(["git", "status", "--porcelain", "--", "src", "mhw_app.py"],
                           cwd=API_DIR, capture_output=True, text=True).stdout.strip()
    return out.stdout.strip() + ("+dirty(src)" if dirty else "")


def cmd_run_once(args) -> int:
    from tests.golden.cases import CASES
    from tests.harness import assert_no_unexpected_client_attempts, describe_env, has_default_client, dask_client_calls

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    only = set(args.only.split(",")) if args.only else None
    results = {}
    with make_client() as client:
        assert not has_default_client(), "a distributed default client exists"
        env = describe_env()
        for case in CASES:
            if only and case["id"] not in only:
                continue
            res = fetch(client, case)
            body_path(out, case, res["status"]).write_bytes(res["body"])
            res.pop("body")
            res["expect_problems"] = expect_problems(case, res) + header_problems(case, res)
            results[case["id"]] = res
            print(f"{case['id']:5s} {res['status']} {res['size']:>10d} {res['sha256'][:12]}"
                  f" {'OK' if not res['expect_problems'] else res['expect_problems']}", flush=True)
        assert not has_default_client(), "a distributed default client appeared during the run"
    assert_no_unexpected_client_attempts()
    env["get_dask_client_calls"] = dask_client_calls["count"]
    (out / "results.json").write_text(json.dumps({"env": env, "results": results}, indent=1, default=str))
    return 0


def zarr_fingerprint() -> dict:
    zpath = (API_DIR / "data" / "mhw.zarr").resolve()
    zmeta = (zpath / ".zmetadata").read_bytes()
    return {"path": str(zpath), "zmetadata_md5": hashlib.md5(zmeta).hexdigest()}


def cmd_record(args) -> int:
    from tests.golden.cases import CASES

    tag_dir = GOLDEN_ROOT / args.tag
    if tag_dir.exists() and not args.force:
        print(f"{tag_dir} exists; use --force to replace", file=sys.stderr)
        return 2
    runs = []
    for i in range(1, args.runs + 1):
        run_dir = tag_dir / f"run{i}"
        print(f"--- run {i} -> {run_dir}", flush=True)
        cmd = [sys.executable, "-m", "tests.golden.runner", "run-once", "--out", str(run_dir)]
        if args.only:
            cmd += ["--only", args.only]
        subprocess.run(cmd, cwd=API_DIR, check=True)
        runs.append(json.loads((run_dir / "results.json").read_text()))

    cases_out = {}
    unstable = []
    for case in CASES:
        cid = case["id"]
        if cid not in runs[0]["results"]:
            continue
        per_run = [r["results"][cid] for r in runs]
        shas = [r["sha256"] for r in per_run]
        stab = stability_problems(per_run)
        stable = not stab
        if not stable:
            unstable.append(cid)
        first = runs[0]["results"][cid]
        src = body_path(tag_dir / "run1", case, first["status"])
        shutil.copy2(src, tag_dir / src.name)
        cases_out[cid] = {"url": case["url"], "kind": case["kind"], "why": case.get("why"),
                          "status": first["status"], "headers": first["headers"],
                          "sha256": first["sha256"], "size": first["size"],
                          "stable": stable, "run_sha256": shas if not stable else None,
                          "stability_problems": stab or None,
                          "expect_problems": first["expect_problems"], "summary": first["summary"]}
    manifest = {
        "tag": args.tag,
        "commit": git_commit(),
        "recorded_at": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "runs": args.runs,
        "zarr": zarr_fingerprint(),
        "env": runs[0]["env"],
        "unstable": unstable,
        "approved_unstable": [],
        "cases": cases_out,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=1, ensure_ascii=False, default=str) + "\n")
    print(f"manifest -> {MANIFEST}; unstable: {unstable or 'none'}")
    bad = {cid: c["expect_problems"] for cid, c in cases_out.items() if c["expect_problems"]}
    if bad:
        print(f"EXPECTATION PROBLEMS: {bad}")
    return 1 if (unstable or bad) else 0


def cmd_amend_headers(args) -> int:
    """Add header evidence to the manifest without touching golden bodies."""
    manifest = load_manifest()
    manifest.setdefault("approved_unstable", [])
    tag = manifest["tag"]
    stamp = _dt.datetime.now().strftime("%Y%m%dT%H%M%S")
    work = GOLDEN_ROOT / f"{tag}_amend_headers_{stamp}"
    runs = []
    for i in range(1, args.runs + 1):
        run_dir = work / f"run{i}"
        print(f"--- amend run {i} -> {run_dir}", flush=True)
        subprocess.run([sys.executable, "-m", "tests.golden.runner", "run-once", "--out", str(run_dir)],
                       cwd=API_DIR, check=True)
        runs.append(json.loads((run_dir / "results.json").read_text()))
    problems = []
    new_headers = {}
    for cid, g in manifest["cases"].items():
        per_run = [r["results"][cid] for r in runs]
        for i, r in enumerate(per_run, start=1):
            if r["sha256"] != g["sha256"]:
                problems.append(f"{cid}: run{i} body differs from the recorded golden")
            if r["status"] != g["status"]:
                problems.append(f"{cid}: run{i} status {r['status']} != {g['status']}")
            if r["expect_problems"]:
                problems.append(f"{cid}: run{i} {r['expect_problems']}")
        problems += [f"{cid}: {p}" for p in stability_problems(per_run)]
        for k, v in g["headers"].items():
            if per_run[0]["headers"].get(k) != v:
                problems.append(f"{cid}: existing header {k} changed ({v!r} -> {per_run[0]['headers'].get(k)!r})")
        new_headers[cid] = per_run[0]["headers"]
    if problems:
        print("REFUSING to amend:\n  " + "\n  ".join(problems))
        return 1
    added = {}
    for cid, h in new_headers.items():
        extra = sorted(set(h) - set(manifest["cases"][cid]["headers"]))
        if extra:
            added[cid] = extra
        manifest["cases"][cid]["headers"] = h
    manifest.setdefault("amendments", []).append({
        "at": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "what": "header evidence added; bodies, statuses and existing headers verified unchanged",
        "runs": args.runs, "commit": git_commit(), "work_dir": str(work), "added": added,
        "env_client_attempts": runs[0]["env"].get("client_attempts"),
    })
    MANIFEST.write_text(json.dumps(manifest, indent=1, ensure_ascii=False, default=str) + "\n")
    print(f"amended {len(added)} case(s): {added}")
    return 0


def main(argv=None) -> int:
    os.chdir(API_DIR)
    if str(API_DIR) not in sys.path:
        sys.path.insert(0, str(API_DIR))
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p1 = sub.add_parser("run-once")
    p1.add_argument("--out", required=True)
    p1.add_argument("--only")
    p2 = sub.add_parser("record")
    p2.add_argument("--tag", required=True)
    p2.add_argument("--runs", type=int, default=3)
    p2.add_argument("--only")
    p2.add_argument("--force", action="store_true")
    p3 = sub.add_parser("amend-headers")
    p3.add_argument("--runs", type=int, default=3)
    args = parser.parse_args(argv)
    return {"run-once": cmd_run_once, "record": cmd_record, "amend-headers": cmd_amend_headers}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
