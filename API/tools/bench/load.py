"""Local concurrency / memory scenarios (spec 16 r3 §5.1). Never point this at production.

    .venv/bin/python -m tools.bench.load --base http://127.0.0.1:18030 --master-pid PID \
        --scenario light|mixed|heavy2|leak|heavy_types|plot_month_large|sub_<key>_x<1|4>
                   |smoke_J01|single_J08_csv_antimeridian|single_J11_global_json|combo_p11 \
        --expected-workers N --out DIR

The last four are the per-step scenarios of the isolated candidate runner
(tools/bench/candidate.py, doc 24 r4); every request records the sha256 of its body.

RSS is sampled every 0.2 s per process.  Every process is recorded with PID,
PPID, name and cmdline (for manual cross-checking) and given a role:

* ``master``       the given PID;
* ``worker``       a direct child of the master whose cmdline equals the
                   master's (a plain forked gunicorn worker) or starts with
                   ``gunicorn: worker`` (setproctitle installed);
* ``unclassified`` any other direct child of the master;
* ``helper``       a deeper descendant (e.g. a subprocess a worker started).

The cmdline rule is a heuristic, not a proof of gunicorn identity (Codex
round 5), so the caller states ``--expected-workers`` and the memory result is
only valid when every sample shows exactly that many workers, the same worker
PIDs throughout, and no unclassified child.  PM2's 4G is a restart threshold,
not a per-worker limit (Codex round 3 #2), so roles are reported separately.

Failure handling (Codex round 4 #1, round 5): an exception in any request
thread, a non-200 status, a timeout, a connection failure or a missing request
makes the HTTP part invalid; a sampler exception, the master disappearing, the
sampler not stopping in time, no samples, or a worker-classification problem
makes the memory part invalid.  Either one makes the run fail: it is written as
``local_<scenario>.FAILED.json`` with the reasons and the process exits
non-zero; it never produces a ``local_<scenario>.json`` result.

``--timeout`` is the per-request HTTP timeout, not a deadline for the whole
scenario (a scenario issues many requests); there is no overall deadline.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import threading
import time
from pathlib import Path

import httpx
import psutil

FULL = "start=1982-01-01&end=2026-10-01"
LIGHT = f"/api/mhw?lon0=135&lat0=15&{FULL}"
HEAVY = {
    "J11_global_json": "/api/mhw?lon0=-180&lon1=180&lat0=-90&lat1=90&start=2026-10-01&end=2026-10-01",
    "plot_month_large": "/api/mhw/plot?bbox=(-90,-60,0,60)&sstvar=sst&period=1982-2011,2012-2021,2022,2023&mode=month",
    "J08_csv_antimeridian": "/api/mhw/csv?lon0=170&lon1=-170&lat0=-5&lat1=5&start=2020-01-01&end=2021-12-01&append=sst",
}
# Largest JSON queries whose guard claim stays just below a threshold, so they
# run unguarded and concurrently (P1-1c, Codex round 4).  The first two are just
# below the 4 MiB P1-0 candidate (kept for the record: they exceeded the 1.5 GB
# line); the last two are just below the 2 MiB default adopted in P1-1c.
SUBTHRESHOLD = {
    "global_1m_sst": "/api/mhw?lon0=-180&lon1=180&lat0=-90&lat1=90&start=2026-10-01&end=2026-10-01&append=sst",  # claim 3.97 MiB
    "box10_10y_3var": "/api/mhw?lon0=130&lon1=140&lat0=10&lat1=20&start=2010-01-01&end=2019-12-01&append=sst,sst_anomaly,td",  # 3.66 MiB
    "half_1m_sst": "/api/mhw?lon0=-180&lon1=-0.25&lat0=-90&lat1=90&start=2026-10-01&end=2026-10-01&append=sst",  # 1.98 MiB
    "box9_10y_2var": "/api/mhw?lon0=130&lon1=139&lat0=10&lat1=19&start=2010-01-01&end=2019-12-01&append=sst,td",  # 1.78 MiB
}
# Smoke sample for the candidate runner: golden J01 (tests/golden/manifest.json;
# tests/test_bench_candidate.py checks that the two values stay equal).
J01 = "/api/mhw?lon0=135&lat0=15&start=1982-01-01&end=2026-10-01"
GOLDEN_J01_SHA256 = "f377522aa11af1cf52d5507b8c773eec9fdf5c622cae9a098aa589efef5b133c"
LEAK = "/api/mhw/plot?bbox=(135,15,137,17)&period=1982-2011,2012-2021,2022&sstvar=sst&mode=month"
MB = 1024 * 1024
DEFAULT_TIMEOUT = 600.0


class LoadError(RuntimeError):
    """A run that must not be reported as a valid performance result."""

    def __init__(self, message: str, details: list | None = None):
        super().__init__(message)
        self.details = details or []


# --------------------------------------------------------------------------- RSS


class RssSampler(threading.Thread):
    def __init__(self, master_pid: int, interval: float = 0.2, stop_timeout: float = 5.0):
        super().__init__(daemon=True)
        self.master = psutil.Process(master_pid)
        self.master_cmdline = self._cmdline(self.master)
        self.interval = interval
        self.stop_timeout = stop_timeout
        self.samples: list[dict] = []
        self.processes: dict[int, dict] = {}
        self.error: dict | None = None
        self.started_at: float | None = None
        self.ended_at: float | None = None
        self.stop_timed_out = False
        self._halt = threading.Event()

    @staticmethod
    def _cmdline(p: psutil.Process) -> list[str]:
        try:
            return p.cmdline()
        except (psutil.AccessDenied, psutil.ZombieProcess, psutil.NoSuchProcess):
            return []

    def _role(self, pid: int, ppid: int, cmd: list[str]) -> str:
        if pid == self.master.pid:
            return "master"
        if ppid != self.master.pid:
            return "helper"
        if cmd and (cmd == self.master_cmdline or cmd[0].startswith("gunicorn: worker")):
            return "worker"
        return "unclassified"

    def _describe(self, p: psutil.Process) -> dict:
        info = self.processes.get(p.pid)
        if info is None:
            cmd = self._cmdline(p)
            ppid = p.ppid()
            info = {"pid": p.pid, "ppid": ppid, "name": p.name(), "cmdline": cmd,
                    "role": self._role(p.pid, ppid, cmd)}
            self.processes[p.pid] = info
        return info

    def snapshot(self) -> dict:
        """One RSS reading of master + descendants.

        A vanished master raises psutil.NoSuchProcess; a descendant that
        vanishes between listing and reading is skipped (the per-sample worker
        count check in memory_problems() catches a missing worker).
        """
        procs = [self.master] + self.master.children(recursive=True)
        snap: dict = {"t": time.time(), "procs": []}
        for p in procs:
            try:
                info = self._describe(p)
                snap["procs"].append({"pid": p.pid, "role": info["role"], "rss_mb": p.memory_info().rss / MB})
            except psutil.NoSuchProcess:
                if p.pid == self.master.pid:
                    raise
                continue
        return snap

    def run(self):
        self.started_at = time.time()
        try:
            while not self._halt.is_set():
                self.samples.append(self.snapshot())
                self._halt.wait(self.interval)
        except Exception as exc:  # reported through memory_problems(), never swallowed
            self.error = {"type": type(exc).__name__, "message": str(exc), "t": time.time(),
                          "samples_before_error": len(self.samples)}
        finally:
            self.ended_at = time.time()

    def stop(self) -> dict:
        self._halt.set()
        if self.is_alive():
            self.join(timeout=self.stop_timeout)
        self.stop_timed_out = self.is_alive()
        return self.peaks()

    def peaks(self) -> dict:
        peaks: dict[str, dict] = {}
        for s in list(self.samples):
            for p in s["procs"]:
                key = f"{p['role']}:{p['pid']}"
                peaks[key] = {"pid": p["pid"], "role": p["role"],
                              "rss_mb": max(peaks.get(key, {}).get("rss_mb", 0.0), p["rss_mb"])}
        return peaks

    def observed_workers(self) -> list[int]:
        return sorted(pid for pid, i in self.processes.items() if i["role"] == "worker")

    def info(self) -> dict:
        samples = list(self.samples)
        return {"started_at": self.started_at, "ended_at": self.ended_at,
                "first_sample_t": samples[0]["t"] if samples else None,
                "last_sample_t": samples[-1]["t"] if samples else None,
                "n_samples": len(samples), "error": self.error,
                "stop_timed_out": self.stop_timed_out, "alive_after_stop": self.is_alive(),
                "observed_workers": self.observed_workers()}

    def memory_problems(self, expected_workers: int | None) -> list[dict]:
        """Reasons the RSS series must not be reported as a valid memory result."""
        probs: list[dict] = []
        if self.error:
            probs.append({"kind": "sampler_error"} | self.error)
        if self.stop_timed_out or self.is_alive():
            probs.append({"kind": "sampler_stop_timeout", "timeout": self.stop_timeout})
        samples = list(self.samples)
        if not samples:
            probs.append({"kind": "no_samples"})
        unclassified = [i for i in self.processes.values() if i["role"] == "unclassified"]
        if unclassified:
            probs.append({"kind": "unclassified_child", "processes": unclassified})
        if expected_workers is None:
            probs.append({"kind": "expected_workers_missing"})
            return probs
        observed = self.observed_workers()
        if len(observed) != expected_workers:
            probs.append({"kind": "worker_count", "expected": expected_workers, "observed": observed})
        for s in samples:
            n = sum(1 for p in s["procs"] if p["role"] == "worker")
            if n != expected_workers:
                probs.append({"kind": "worker_count_in_sample", "t": s["t"], "expected": expected_workers, "got": n})
                break
        return probs


def rss_by_pid(snap: dict, role: str | None = None) -> dict[int, float]:
    return {p["pid"]: p["rss_mb"] for p in snap["procs"] if role is None or p["role"] == role}


# ---------------------------------------------------------------------- requests


def timed_get(client: httpx.Client, path: str, save_to: Path | None = None) -> dict:
    """One request with wall-clock start/end and the body's sha256; transport
    errors propagate.  The body is written to ``save_to`` only when given."""
    start = time.time()
    t0 = time.perf_counter()
    r = client.get(path)
    body = r.content
    out = {"path": path, "status": r.status_code, "start": start, "end": time.time(),
           "seconds": time.perf_counter() - t0, "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest()}
    if save_to is not None:
        Path(save_to).write_bytes(body)
        out["saved_to"] = str(save_to)
    return out


def check_results(results: list[dict], expected: int, errors: list[dict], expect_status: int = 200) -> None:
    """Raise LoadError unless exactly ``expected`` requests all returned ``expect_status``."""
    problems = list(errors)
    bad = [r for r in results if r["status"] != expect_status]
    problems += [{"kind": "status", "path": r["path"], "status": r["status"]} for r in bad]
    if len(results) != expected:
        problems.append({"kind": "count", "expected": expected, "got": len(results)})
    if problems:
        raise LoadError(f"{len(problems)} problem(s) in load run", problems)


def run_threads(base: str, jobs: list[list[str]], timeout: float = DEFAULT_TIMEOUT) -> list[dict]:
    """Each inner list runs sequentially in its own thread.

    A thread that hits an exception records it and stops; after all threads are
    joined the run is validated (count and status) and LoadError is raised on
    any problem, so no partial result can pass as valid.
    """
    out: list[dict] = []
    errors: list[dict] = []
    lock = threading.Lock()

    def worker(paths, tid):
        try:
            with httpx.Client(base_url=base, timeout=timeout) as c:
                for p in paths:
                    res = timed_get(c, p) | {"thread": tid}
                    with lock:
                        out.append(res)
        except Exception as exc:  # transported back to the main thread
            with lock:
                errors.append({"kind": "exception", "thread": tid, "type": type(exc).__name__, "error": str(exc)})

    threads = [threading.Thread(target=worker, args=(paths, i), daemon=True) for i, paths in enumerate(jobs)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    check_results(out, sum(len(j) for j in jobs), errors)
    return sorted(out, key=lambda r: r["start"])


# --------------------------------------------------------------------- scenarios


def scenario_light(base, sampler, timeout=DEFAULT_TIMEOUT):
    res = {}
    for c in (1, 2, 4):
        t0 = time.perf_counter()
        r = run_threads(base, [[LIGHT] * (8 // c) for _ in range(c)], timeout=timeout)
        res[str(c)] = {"wall": time.perf_counter() - t0, "requests": r}
    return res


def scenario_mixed(base, sampler, heavy_path=HEAVY["J11_global_json"], light_path=LIGHT,
                   timeout=DEFAULT_TIMEOUT, min_light_per_thread=3, max_light_per_thread=200):
    """1 heavy + 3 light threads that keep asking while the heavy one runs.

    ``done`` is set in ``finally`` whether the heavy request succeeds or fails;
    light threads also stop at their first error or after max_light_per_thread.
    """
    done = threading.Event()
    lock = threading.Lock()
    light: list[dict] = []
    heavy: list[dict] = []
    errors: list[dict] = []

    def heavy_job():
        try:
            with httpx.Client(base_url=base, timeout=timeout) as c:
                res = timed_get(c, heavy_path)
                with lock:
                    heavy.append(res)
        except Exception as exc:
            with lock:
                errors.append({"kind": "exception", "thread": "heavy", "type": type(exc).__name__, "error": str(exc)})
        finally:
            done.set()

    def light_job(tid):
        n = 0
        try:
            with httpx.Client(base_url=base, timeout=timeout) as c:
                while (not done.is_set() or n < min_light_per_thread) and n < max_light_per_thread:
                    during = not done.is_set()
                    res = timed_get(c, light_path) | {"thread": tid, "during_heavy": during}
                    with lock:
                        light.append(res)
                    n += 1
                    if res["status"] != 200:
                        break
        except Exception as exc:
            with lock:
                errors.append({"kind": "exception", "thread": tid, "type": type(exc).__name__, "error": str(exc)})

    ts = [threading.Thread(target=heavy_job, daemon=True)] + \
         [threading.Thread(target=light_job, args=(i,), daemon=True) for i in range(3)]
    for t in ts:
        t.start()
    for t in ts:
        t.join()
    check_results(heavy, 1, errors)
    check_results(light, len(light), [])  # every light request must be 200
    per_thread = {i: sum(1 for r in light if r["thread"] == i) for i in range(3)}
    if min(per_thread.values()) < min_light_per_thread:
        raise LoadError("light threads stopped early", [{"kind": "count", "per_thread": per_thread}])
    return {"heavy": heavy[0], "light": sorted(light, key=lambda r: r["start"])}


def scenario_heavy2(base, sampler, timeout=DEFAULT_TIMEOUT):
    """Two heavy requests racing for the same worker (run the server with -w 1)."""
    t0 = time.perf_counter()
    r = run_threads(base, [[HEAVY["J11_global_json"]], [HEAVY["J11_global_json"]]], timeout=timeout)
    return {"wall": time.perf_counter() - t0, "requests": r}


def scenario_leak(base, sampler, timeout=DEFAULT_TIMEOUT):
    """20 sequential month-mode plots (run with -w 1); RSS read ~0.5 s after each."""
    rss = []
    with httpx.Client(base_url=base, timeout=timeout) as c:
        for i in range(20):
            r = timed_get(c, LEAK)
            check_results([r], 1, [])
            time.sleep(0.5)
            snap = sampler.snapshot()
            rss.append({"i": i, "request": r, "worker_rss_mb": rss_by_pid(snap, "worker")})
    return {"rss_after_each": rss}


def _measure_one(c, sampler, path, settle=1.0):
    """One request with worker RSS before, peak during and ~settle s after."""
    time.sleep(settle)
    before = sampler.snapshot()
    n0 = len(sampler.samples)
    r = timed_get(c, path)
    check_results([r], 1, [])
    time.sleep(settle)
    during = sampler.samples[n0:]
    after = sampler.snapshot()
    b = rss_by_pid(before, "worker")
    peak = {pid: max([rss_by_pid(s, "worker").get(pid, 0.0) for s in during] + [rss_by_pid(after, "worker").get(pid, 0.0)])
            for pid in b}
    return {"request": r, "before": b, "peak": peak, "after_1s": rss_by_pid(after, "worker"),
            "peak_minus_before": {pid: peak[pid] - b[pid] for pid in b}, "n_samples_during": len(during)}


def scenario_heavy_types(base, sampler, timeout=DEFAULT_TIMEOUT):
    """One request per heavy type, sequential, in ONE worker (-w 1).

    Note: types run back to back, so later types may reuse memory freed (or not
    returned to the OS) by earlier ones; only the first type is an isolated peak.
    """
    out = {}
    with httpx.Client(base_url=base, timeout=timeout) as c:
        for order, (name, path) in enumerate(HEAVY.items()):
            out[name] = {"order": order} | _measure_one(c, sampler, path)
    return out


def scenario_plot_month_large(base, sampler, timeout=DEFAULT_TIMEOUT):
    """The large month plot alone: run against a freshly started -w 1 server so
    the peak is independent of earlier heavy requests."""
    with httpx.Client(base_url=base, timeout=timeout) as c:
        return {"plot_month_large": _measure_one(c, sampler, HEAVY["plot_month_large"])}


def _measure_concurrent(base, sampler, path, n, timeout):
    """n identical requests started together on one worker; worker RSS before / peak / after."""
    time.sleep(1.0)
    before = sampler.snapshot()
    n0 = len(sampler.samples)
    t0 = time.perf_counter()
    reqs = run_threads(base, [[path] for _ in range(n)], timeout=timeout)
    wall = time.perf_counter() - t0
    time.sleep(1.0)
    during = sampler.samples[n0:]
    after = sampler.snapshot()
    b = rss_by_pid(before, "worker")
    peak = {pid: max([rss_by_pid(s, "worker").get(pid, 0.0) for s in during] + [rss_by_pid(after, "worker").get(pid, 0.0)])
            for pid in b}
    return {"path": path, "n": n, "wall": wall, "requests": reqs, "before": b, "peak": peak,
            "after_1s": rss_by_pid(after, "worker"), "peak_minus_before": {pid: peak[pid] - b[pid] for pid in b},
            "n_samples_during": len(during)}


def scenario_smoke_J01(base, sampler, timeout=DEFAULT_TIMEOUT, out_dir: Path | None = None):
    """Golden J01 once; its sha256 must equal the golden value (a smoke sample, not a
    check of the whole store).  Only this scenario keeps a response body."""
    with httpx.Client(base_url=base, timeout=timeout) as c:
        r = timed_get(c, J01, save_to=(Path(out_dir) / "j01.json") if out_dir else None)
    check_results([r], 1, [])
    if r["sha256"] != GOLDEN_J01_SHA256:
        raise LoadError("smoke J01 sha256 differs from golden",
                        [{"kind": "sha256", "got": r["sha256"], "expected": GOLDEN_J01_SHA256}])
    return {"J01": r, "golden_sha256": GOLDEN_J01_SHA256}


def _single(key):
    def scenario(base, sampler, timeout=DEFAULT_TIMEOUT):
        with httpx.Client(base_url=base, timeout=timeout) as c:
            return {key: _measure_one(c, sampler, HEAVY[key])}
    scenario.__doc__ = f"{key} alone on a freshly started -w 1 server."
    return scenario


COMBO_LIGHT = ["half_1m_sst"] * 4


def scenario_combo_p11(base, sampler, timeout=DEFAULT_TIMEOUT):
    """1 large month plot + 4 half-globe 1-month sst JSON (just below 2 MiB) + 1 J01,
    all started together on one worker.  Measured directly, never summed from
    single runs.  Valid only if every other request's interval overlaps the plot's."""
    paths = [HEAVY["plot_month_large"]] + [SUBTHRESHOLD[k] for k in COMBO_LIGHT] + [J01]
    time.sleep(1.0)
    before = sampler.snapshot()
    n0 = len(sampler.samples)
    t0 = time.perf_counter()
    reqs = run_threads(base, [[p] for p in paths], timeout=timeout)
    wall = time.perf_counter() - t0
    time.sleep(1.0)
    during = sampler.samples[n0:]
    after = sampler.snapshot()
    b = rss_by_pid(before, "worker")
    peak = {pid: max([rss_by_pid(s, "worker").get(pid, 0.0) for s in during] + [rss_by_pid(after, "worker").get(pid, 0.0)])
            for pid in b}
    plot = next(r for r in reqs if r["path"] == paths[0])
    others = [r for r in reqs if r is not plot]
    overlap = [{"path": r["path"], "start": r["start"], "end": r["end"],
                "overlaps_plot": r["start"] < plot["end"] and r["end"] > plot["start"]} for r in others]
    if not all(o["overlaps_plot"] for o in overlap):
        raise LoadError("combo requests did not overlap the plot", [{"kind": "no_overlap", "plot": plot,
                                                                       "others": overlap}])
    return {"combo_p11": {"paths": paths, "wall": wall, "requests": reqs, "overlap": overlap,
                          "before": b, "peak": peak, "after_1s": rss_by_pid(after, "worker"),
                          "peak_minus_before": {pid: peak[pid] - b[pid] for pid in b},
                          "n_samples_during": len(during)}}


def _subthreshold(key, n):
    def scenario(base, sampler, timeout=DEFAULT_TIMEOUT):
        return {f"{key}_x{n}": _measure_concurrent(base, sampler, SUBTHRESHOLD[key], n, timeout)}
    scenario.__doc__ = f"{n} x {key} at once on a fresh -w 1 worker (sub-threshold JSON)."
    return scenario


SCENARIOS = {"light": scenario_light, "mixed": scenario_mixed, "heavy2": scenario_heavy2,
             "leak": scenario_leak, "heavy_types": scenario_heavy_types,
             "plot_month_large": scenario_plot_month_large}
SCENARIOS.update({f"sub_{k}_x{n}": _subthreshold(k, n) for k in SUBTHRESHOLD for n in (1, 4)})
SCENARIOS.update({"smoke_J01": scenario_smoke_J01,
                  "single_J08_csv_antimeridian": _single("J08_csv_antimeridian"),
                  "single_J11_global_json": _single("J11_global_json"),
                  "combo_p11": scenario_combo_p11})


def run_scenario(name: str, base: str, master_pid: int, out_dir: Path, timeout: float = DEFAULT_TIMEOUT,
                 scenario_fn=None, expected_workers: int | None = None, sampler_kwargs: dict | None = None) -> int:
    """Run one scenario with RSS sampling; the sampler is always stopped.

    Returns 0 and writes local_<name>.json only when both the HTTP part and the
    memory part are valid; otherwise returns 1 and writes
    local_<name>.FAILED.json with every reason.
    """
    fn = scenario_fn or SCENARIOS[name]
    out_dir.mkdir(parents=True, exist_ok=True)
    sampler = RssSampler(master_pid, **(sampler_kwargs or {}))
    sampler.start()
    doc = {"scenario": name, "base": base, "started": time.time(), "expected_workers": expected_workers,
           "http_timeout_s": timeout, "note": "http_timeout_s is per request, not a total scenario deadline"}
    try:
        time.sleep(0.3)
        doc["result"] = fn(base, sampler, timeout=timeout)
        doc["http_valid"] = True
    except Exception as exc:
        doc["http_valid"] = False
        doc["error"] = {"type": type(exc).__name__, "message": str(exc), "details": getattr(exc, "details", [])}
    finally:
        doc["rss_peak_mb"] = sampler.stop()
        doc["ended"] = time.time()
        doc["sampler"] = sampler.info()
        doc["processes"] = list(sampler.processes.values())
        doc["rss_series_mb"] = list(sampler.samples)
        doc["sampler_alive_after_stop"] = sampler.is_alive()
        doc["memory_problems"] = sampler.memory_problems(expected_workers)
        doc["memory_valid"] = not doc["memory_problems"]
    doc["valid"] = doc["http_valid"] and doc["memory_valid"]
    path = out_dir / (f"local_{name}.json" if doc["valid"] else f"local_{name}.FAILED.json")
    path.write_text(json.dumps(doc, indent=1, default=str))
    summary = {"scenario": name, "valid": doc["valid"], "http_valid": doc["http_valid"],
               "memory_valid": doc["memory_valid"], "file": str(path),
               "observed_workers": doc["sampler"]["observed_workers"],
               "rss_peak_mb": {k: round(v["rss_mb"]) for k, v in doc["rss_peak_mb"].items()}}
    if not doc["http_valid"]:
        summary["error"] = doc["error"]
    if doc["memory_problems"]:
        summary["memory_problems"] = [p["kind"] for p in doc["memory_problems"]]
    print(json.dumps(summary, default=str))
    return 0 if doc["valid"] else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default="http://127.0.0.1:18030")
    ap.add_argument("--master-pid", type=int, required=True)
    ap.add_argument("--scenario", choices=SCENARIOS, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--expected-workers", type=int, required=True,
                    help="number of gunicorn workers (-w); a mismatch invalidates the memory result")
    ap.add_argument("--timeout", type=float, default=DEFAULT_TIMEOUT, help="per-request HTTP timeout (s)")
    args = ap.parse_args(argv)
    return run_scenario(args.scenario, args.base, args.master_pid, Path(args.out), timeout=args.timeout,
                        expected_workers=args.expected_workers)


if __name__ == "__main__":
    sys.exit(main())
