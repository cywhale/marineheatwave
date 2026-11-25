"""P1-2 V-c: isolated API-only validation on vm24 (owner 2026-11-25).

Benchmark tooling only; never imported by the API.  Every server is a fresh
candidate started by the stage A runner (tools/bench/candidate.run_step: proven-
identity stop, exit / port check, watchdog), bound to 127.0.0.1, with an
explicit worker class and the HTTP readiness probe (r20).

    # one environment (A = plain polars, B = polars[rtcompat]); stops at the first
    # step that is not exit 0
    python -B -m tools.p1_2.vc_server steps --dir D --port P --python PY --gunicorn G \\
        --label B --steps echo,golden,smoke,j08,j11,plot,combo [--extra-env POLARS_FORCE_PKG=compat] \\
        [--step-suffix _forcecompat]

    # paired latency (doc 29 r3 §5.1): every candidate started (nested, each its
    # own run_step), then the same light queries alternated over production
    # (L2 loopback) and the candidates
    python -B -m tools.p1_2.vc_server paired --spec A:D_A:PORT_A:PY_A:G_A --spec B:... \\
        --prod https://127.0.0.1:8030 --rounds 30 --out F

Gate before every server start: MemAvailable >= 8 GiB, port free, no
update_mhw process (doc 24 r6 §3.3 / 26).  Thresholds: worker RSS observation
1.5 GB, abort 3.5 GB or MemAvailable < 6 GiB.

Each step also records which polars runtime every worker actually mapped
(/proc/<worker>/maps, Linux) and counts the polars "Missing required CPU
features" warning in the step's logs.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
import os
import random
import statistics
import sys
import time
import uuid
from pathlib import Path

import httpx
import psutil

from tools.bench import candidate as cand
from tools.bench import load
from tools.bench import procident as pi
from tools.p1_2 import vb_server as vb

GIB = 1024 ** 3
GATE_MIN_AVAILABLE = 8 * GIB
LIMITS = dict(observe_mb=1500, abort_rss_mb=3500, abort_avail_mb=6144)
CPU_WARNING = "Missing required CPU features"
EXIT_GATE = 9
STAGE_B = ("smoke", "j08", "j11", "plot", "combo")
# light point queries of doc 27 §4.2 / latency_layers (Q1 = golden J01)
QUERIES = {"Q1": load.J01, "Q3": "/api/mhw?lon0=135.3&lat0=15.3&start=1982-01-01&end=2026-10-01"}


# ------------------------------------------------------------------ helpers


def gate(port: int) -> list[str]:
    """Reasons not to start a server now (empty = go)."""
    probs = []
    avail = psutil.virtual_memory().available
    if avail < GATE_MIN_AVAILABLE:
        probs.append(f"MemAvailable {avail / GIB:.1f} GiB < 8 GiB")
    if not pi.port_free(port):
        probs.append(f"port {port} in use")
    me = os.getpid()
    for p in psutil.process_iter(["pid", "cmdline"]):
        try:
            if p.info["pid"] != me and any("update_mhw" in x for x in (p.info["cmdline"] or [])):
                probs.append(f"update_mhw running (pid {p.info['pid']})")
        except psutil.Error:
            pass
    return probs


def polars_maps(pids) -> dict:
    """pid -> sorted polars runtime libraries mapped in that process (Linux)."""
    out = {}
    for pid in pids:
        try:
            with open(f"/proc/{pid}/maps", encoding="utf-8", errors="replace") as f:
                libs = {line.split()[-1] for line in f if "_polars_runtime" in line and line.rstrip().endswith(".so")}
            out[str(pid)] = sorted(libs)
        except OSError as exc:
            out[str(pid)] = f"unreadable: {type(exc).__name__}"
    return out


def runtime_of(libs) -> str | None:
    if not isinstance(libs, list) or len(libs) != 1:
        return None
    for name in ("compat", "64", "32"):
        if f"/_polars_runtime_{name}/" in libs[0]:
            return name
    return None


def count_in(path: Path, text: str) -> int:
    try:
        return path.read_text(errors="replace").count(text)
    except FileNotFoundError:
        return 0


def with_maps(fn, expect_runtime: str | None):
    """Wrap a load scenario: after it, record the polars runtime of every worker."""
    def scenario(base, sampler, timeout=load.DEFAULT_TIMEOUT):
        out = fn(base, sampler, timeout=timeout)
        workers = [p.pid for p in psutil.Process(sampler.master.pid).children()
                   if p.pid in sampler.observed_workers()]
        maps = polars_maps(workers)
        rts = {pid: runtime_of(v) for pid, v in maps.items()}
        out = dict(out or {})
        out["polars_maps"], out["polars_runtime"] = maps, rts
        if expect_runtime and (not rts or any(r != expect_runtime for r in rts.values())):
            raise load.LoadError("polars runtime differs from the expected one", [{"expected": expect_runtime, **rts}])
        return out
    return scenario


# -------------------------------------------------------------------- steps


def run_one(D: Path, step: str, a, extra_env: dict, expect_runtime: str | None) -> dict:
    name = step + a.step_suffix
    problems = gate(a.port)
    rec = {"step": name, "gate": problems, "t_start": time.time()}
    if problems:
        rec["exit"] = EXIT_GATE
        return rec
    common = dict(python=a.python, gunicorn=a.gunicorn, worker_class=a.worker_class, extra_env=extra_env,
                  ready_log=True, **LIMITS)
    if step == "echo":
        argv = cand.gunicorn_argv(a.python, a.gunicorn, a.port, D, D / "run" / name, name,
                                  worker_class=a.worker_class, workers=1, app="tools.p1_2.asgi_echo:app")
        argv[4:4] = ["--pythonpath", str(Path(__file__).resolve().parents[2])]  # the echo app lives with the tools
        rc = cand.run_step(D, name, a.port, argv=argv, expected_workers=1, scenario_fn=vb.scenario_echo,
                           probe_path="/", ready_log=cand.error_log(D, name), scenario_name=f"asgi_echo{a.step_suffix}",
                           python=a.python, extra_env=extra_env, **LIMITS)
    elif step == "golden":
        rc = cand.run_step(D, name, a.port, expected_workers=2, workers=2, step_deadline=1800,
                           scenario_fn=with_maps(vb.make_scenario_golden(D, a.port), expect_runtime),
                           scenario_name=f"vc_golden{a.step_suffix}", **common)
    elif step in STAGE_B:
        fn = load.SCENARIOS[cand.STEPS[step]]
        if step == "smoke":
            fn = functools.partial(load.scenario_smoke_J01, out_dir=D / "out")
        rc = cand.run_step(D, name, a.port, expected_workers=1, workers=1, scenario_fn=with_maps(fn, expect_runtime),
                           scenario_name=f"{cand.STEPS[step]}{a.step_suffix}", **common)
    else:
        raise SystemExit(f"unknown step {step}")
    rec["exit"] = rc
    rec["t_end"] = time.time()
    rec["cpu_warning_count"] = {"error_log": count_in(cand.error_log(D, name), CPU_WARNING),
                                "candidate_out": count_in(D / "logs" / f"candidate_{name}.out", CPU_WARNING)}
    try:
        res = json.loads((D / "run" / name / "result.json").read_text())
        rec["exit_check"] = res.get("exit_check")
        rec["ready"] = res.get("ready")
        rec["peak_worker_rss_mb"] = res.get("monitor", {}).get("peak_worker_rss_mb")
        rec["min_available_mb"] = res.get("monitor", {}).get("min_available_mb")
    except (OSError, ValueError) as exc:
        rec["result_error"] = str(exc)
    return rec


def cmd_steps(a) -> int:
    D = Path(a.dir).resolve()
    extra_env = dict(kv.split("=", 1) for kv in a.extra_env)
    summary = {"label": a.label, "dir": str(D), "extra_env": extra_env, "worker_class": a.worker_class,
               "expect_runtime": a.expect_runtime, "steps": []}
    out = D / f"vc_steps_{a.label}{a.step_suffix}.json"
    status = 0
    for step in a.steps.split(","):
        rec = run_one(D, step, a, extra_env, a.expect_runtime)
        summary["steps"].append(rec)
        print(json.dumps({k: rec.get(k) for k in ("step", "exit", "gate", "exit_check", "peak_worker_rss_mb",
                                                   "cpu_warning_count")}, default=str), flush=True)
        out.write_text(json.dumps(summary, indent=1, default=str))
        if rec["exit"] != 0:
            status = rec["exit"]
            summary["stopped_after"] = rec["step"]
            break
    out.write_text(json.dumps(summary, indent=1, default=str))
    return status


# ------------------------------------------------------------------- paired


def bootstrap_ci(diffs: list[float], n: int = 10000, seed: int = 20261125) -> tuple[float, float]:
    rng = random.Random(seed)
    meds = sorted(statistics.median(rng.choices(diffs, k=len(diffs))) for _ in range(n))
    return meds[int(0.025 * n)], meds[int(0.975 * n) - 1]


def describe(xs: list[float]) -> dict:
    q = statistics.quantiles(xs, n=4)
    return {"n": len(xs), "median": statistics.median(xs), "p25": q[0], "p75": q[2], "min": min(xs), "max": max(xs)}


def measure_paired(bases: dict[str, str], rounds: int, warmup: int = 1, timeout: float = 60.0) -> dict:
    """Same queries, labels rotated every round (each label first equally often).
    Every response must be 200 and byte-identical across labels (Q1 also equals
    golden J01).  Times are client wall time to the full body (persistent
    connections; production is TLS on loopback, the candidates plain HTTP)."""
    labels = list(bases)
    clients = {lab: httpx.Client(base_url=b, timeout=timeout, verify=False) for lab, b in bases.items()}
    rows, problems = [], []

    def one(lab, q, r, phase):
        path = f"{QUERIES[q]}&_nc=vc{uuid.uuid4().hex[:10]}"
        t = time.perf_counter()
        resp = clients[lab].get(path)
        dt = time.perf_counter() - t
        row = {"round": r, "phase": phase, "query": q, "label": lab, "status": resp.status_code,
               "seconds": dt, "bytes": len(resp.content), "sha256": hashlib.sha256(resp.content).hexdigest(),
               "at": time.time()}
        rows.append(row)
        return row

    try:
        for w in range(warmup):
            for q in QUERIES:
                for lab in labels:
                    one(lab, q, w, "warmup")
        for r in range(rounds):
            for q in QUERIES:
                k = r % len(labels)
                order = labels[k:] + labels[:k]
                got = [one(lab, q, r, "measure") for lab in order]
                if any(g["status"] != 200 for g in got) or len({g["sha256"] for g in got}) != 1:
                    problems.append({"round": r, "query": q, "rows": got})
                if q == "Q1" and got[0]["sha256"] != load.GOLDEN_J01_SHA256:
                    problems.append({"round": r, "query": q, "problem": "Q1 != golden J01"})
    finally:
        for c in clients.values():
            c.close()
    stats: dict = {"per_label": {}, "pairs": {}}
    meas = [x for x in rows if x["phase"] == "measure"]
    for q in QUERIES:
        by = {lab: {x["round"]: x["seconds"] for x in meas if x["query"] == q and x["label"] == lab} for lab in labels}
        stats["per_label"][q] = {lab: describe(list(v.values())) for lab, v in by.items()}
        for i, l1 in enumerate(labels):
            for l2 in labels[i + 1:]:
                d = [by[l2][r] - by[l1][r] for r in sorted(by[l1]) if r in by[l2]]
                lo, hi = bootstrap_ci(d)
                stats["pairs"].setdefault(q, {})[f"{l2}-{l1}"] = {"n": len(d), "median_diff": statistics.median(d),
                                                                   "ci95": [lo, hi]}
    return {"labels": labels, "bases": bases, "rounds": rounds, "warmup": warmup, "rows": rows,
            "stats": stats, "problems": problems}


def cmd_paired(a) -> int:
    specs = []
    for s in a.spec:
        label, d, port, py, gun = s.split(":", 4)
        specs.append({"label": label, "D": Path(d).resolve(), "port": int(port), "python": py, "gunicorn": gun})
    out: dict = {"prod": a.prod, "specs": [dict(s, D=str(s["D"])) for s in specs], "bases": {}, "inner_rc": {},
                 "gates": {}}
    if a.prod:
        out["bases"]["prod"] = a.prod

    def start(i: int) -> int:
        s = specs[i]
        probs = gate(s["port"])
        out["gates"][s["label"]] = probs
        if probs:
            return EXIT_GATE

        def scenario(base, sampler, timeout=load.DEFAULT_TIMEOUT):
            out["bases"][s["label"]] = base
            if i + 1 < len(specs):
                rc = start(i + 1)
                out["inner_rc"][specs[i + 1]["label"]] = rc
                if rc != 0:
                    raise load.LoadError(f"candidate {specs[i + 1]['label']} exit {rc}")
            else:
                out["measure"] = measure_paired(dict(out["bases"]), a.rounds)
                if out["measure"]["problems"]:
                    raise load.LoadError("paired responses differ or failed", out["measure"]["problems"][:5])
            return {"label": s["label"]}

        return cand.run_step(s["D"], f"paired{a.step_suffix}", s["port"], python=s["python"], gunicorn=s["gunicorn"],
                             worker_class=a.worker_class, expected_workers=1, workers=1, scenario_fn=scenario,
                             scenario_name=f"vc_paired_{s['label']}", install_handlers=(i == 0),
                             step_deadline=1800, **LIMITS)

    rc = start(0)
    out["exit"] = rc
    Path(a.out).write_text(json.dumps(out, indent=1, default=str))
    m = out.get("measure") or {}
    print(json.dumps({"exit": rc, "inner_rc": out["inner_rc"], "gates": out["gates"],
                      "stats": m.get("stats"), "problems": len(m.get("problems", []))}, indent=1, default=str))
    return rc


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("steps")
    s.add_argument("--dir", required=True)
    s.add_argument("--port", type=int, required=True)
    s.add_argument("--python", required=True)
    s.add_argument("--gunicorn", required=True)
    s.add_argument("--label", required=True)
    s.add_argument("--steps", default="echo,golden," + ",".join(STAGE_B))
    s.add_argument("--extra-env", action="append", default=[])
    s.add_argument("--expect-runtime", choices=("compat", "32", "64"), default=None)
    s.add_argument("--step-suffix", default="")
    s.add_argument("--worker-class", default=cand.DEFAULT_WORKER_CLASS)
    p = sub.add_parser("paired")
    p.add_argument("--spec", action="append", required=True, help="LABEL:DIR:PORT:PYTHON:GUNICORN")
    p.add_argument("--prod", default=None, help="production L2 loopback base URL (optional)")
    p.add_argument("--rounds", type=int, default=30)
    p.add_argument("--out", required=True)
    p.add_argument("--step-suffix", default="")
    p.add_argument("--worker-class", default=cand.DEFAULT_WORKER_CLASS)
    a = ap.parse_args(argv)
    return cmd_steps(a) if a.cmd == "steps" else cmd_paired(a)


if __name__ == "__main__":
    sys.exit(main())
