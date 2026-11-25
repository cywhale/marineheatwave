"""Isolated candidate runner for the vm24 validation (doc 24 r4 §1.2-§1.5).

Benchmark tooling only; never imported by the API.  One ``run`` = one freshly
started single-worker candidate gunicorn, one measured step, a watchdog, and a
stop + exit check that only signals processes whose identity is proven
(tools/bench/procident.py).

    python -B -m tools.bench.candidate run --dir D --step smoke|j08|j11|plot|combo --port 18031 \\
        --python PY --gunicorn GUNICORN_SCRIPT \\
        --observe-worker-rss-mb 1500 --abort-worker-rss-mb 3500 --abort-min-available-mb 6144
    python -B -m tools.bench.candidate recover --dir D --mode list|stop --port 18031

D layout: app/API (cwd of the candidate: mhw_app.py, src/, data/mhw.zarr),
run/<step>/ (launcher.json, candidate.id, candidate.pids, gunicorn.pid,
result.json; a step directory that already exists is refused), logs/, out/,
tmp/, home/, xdg/{cache,config,data,state}/, mpl/.  The candidate environment
is built from scratch (nothing inherited), so no DASK_* or other variable leaks.

Exit status: 0 ok, 1 HTTP/result invalid, 2 usage/state-file format, 3 observation
line crossed (do not start the next step), 4 watchdog abort, 5 exit check failed
(candidate process left or port still bound), 6 identity undetermined (nothing
signalled), 7 a recorded number now belongs to another process (nothing
signalled for it).  Precedence: 6, 5, 4, 2, 1, 3, 7.
"""
from __future__ import annotations

import argparse
import ctypes
import functools
import json
import os
import signal
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path

import psutil

from tools.bench import load
from tools.bench import procident as pi

STEPS = {"smoke": "smoke_J01", "j08": "single_J08_csv_antimeridian", "j11": "single_J11_global_json",
         "plot": "plot_month_large", "combo": "combo_p11"}
EXIT_OK, EXIT_FAILED, EXIT_FORMAT, EXIT_OBSERVE, EXIT_ABORT, EXIT_RESIDUAL, EXIT_UNDETERMINED, EXIT_MISMATCH = range(8)
PRECEDENCE = [EXIT_UNDETERMINED, EXIT_RESIDUAL, EXIT_ABORT, EXIT_FORMAT, EXIT_FAILED, EXIT_OBSERVE, EXIT_MISMATCH]
SUBDIRS = ("app/API", "run", "logs", "out", "tmp", "home", "xdg/cache", "xdg/config", "xdg/data", "xdg/state", "mpl")
MB = 1024 * 1024
_SIGS = (signal.SIGTERM, signal.SIGHUP, signal.SIGINT)
AVAILABLE_SOURCE = "psutil.virtual_memory().available (Linux: /proc/meminfo MemAvailable)"


def worst(codes) -> int:
    codes = set(codes)
    for c in PRECEDENCE:
        if c in codes:
            return c
    return EXIT_OK


class Interrupted(BaseException):
    """SIGTERM / SIGHUP / SIGINT delivered to the runner (cleanup still runs)."""

    def __init__(self, signum: int):
        super().__init__(signal.Signals(signum).name)
        self.signum = signum


def candidate_env(D: Path, python: str) -> dict[str, str]:
    d = str(D)
    return {"PATH": f"{Path(python).parent}:/usr/bin:/bin", "LANG": "C.UTF-8",
            "HOME": f"{d}/home", "XDG_CACHE_HOME": f"{d}/xdg/cache", "XDG_CONFIG_HOME": f"{d}/xdg/config",
            "XDG_DATA_HOME": f"{d}/xdg/data", "XDG_STATE_HOME": f"{d}/xdg/state",
            "TMPDIR": f"{d}/tmp", "MPLCONFIGDIR": f"{d}/mpl", "MPLBACKEND": "Agg",
            "PYTHONDONTWRITEBYTECODE": "1", "PYTHONNOUSERSITE": "1",
            "MHW_HEAVY_SLOTS": "1", "MHW_HEAVY_BYTES": "2097152"}


def gunicorn_argv(python: str, gunicorn: str, port: int, D: Path, run_dir: Path, step: str) -> list[str]:
    return [python, "-B", gunicorn, "mhw_app:app", "-w", "1", "-k", "uvicorn.workers.UvicornWorker",
            "-b", f"127.0.0.1:{port}", "--timeout", "600", "--graceful-timeout", "30",
            "--pid", str(run_dir / "gunicorn.pid"), "--worker-tmp-dir", str(D / "tmp"),
            "--error-logfile", str(D / "logs" / f"gunicorn_{step}.err"),
            "--access-logfile", str(D / "logs" / f"access_{step}.log")]


# ------------------------------------------------------------------------ launcher


class Launcher:
    """Starts the candidate in its own session; the child records its own identity
    (candidate.id) before exec, so no candidate code ever runs without it."""

    def __init__(self, argv: list[str], cwd: Path, env: dict, run_dir: Path, *, step: str, port: int,
                 table=None, log_path: Path | None = None, child_pause_before_id: float = 0.0):
        self.argv, self.cwd, self.env, self.run_dir = argv, Path(cwd), env, Path(run_dir)
        self.step, self.port = step, port
        self.table = table or pi.platform_table()
        self.log_path = log_path
        self.child_pause_before_id = child_pause_before_id  # tests only (V8)
        self.proc: subprocess.Popen | None = None
        self.cid: pi.CandidateId | None = None
        self.id_problem: str | None = None
        self.launcher_pid = os.getpid()
        self._libc = ctypes.CDLL(None, use_errno=True) if sys.platform.startswith("linux") else None

    # runs in the child after setsid(), before exec: os-level calls only
    def _child(self):
        for s in _SIGS:
            signal.signal(s, signal.SIG_DFL)
        if self.child_pause_before_id:  # tests only (V8): announce the pause, then wait
            pi.write_atomic(str(self.run_dir / ".test_child_paused"), f"{os.getpid()}\n")
            time.sleep(self.child_pause_before_id)
        me = self.table.read(os.getpid())
        if me is None or not (me.pgid == me.sid == me.pid):
            os._exit(71)
        cid = pi.CandidateId(self.table.unit, me.pgid, me.sid, me.pid, me.start)
        pi.write_atomic(str(self.run_dir / pi.ID_NAME), cid.line())
        if self._libc is not None and self._libc.prctl(1, int(signal.SIGTERM), 0, 0, 0) != 0:  # PR_SET_PDEATHSIG
            os._exit(72)
        if os.getppid() != self.launcher_pid:  # the runner is already gone: never exec
            os._exit(70)
        signal.pthread_sigmask(signal.SIG_UNBLOCK, _SIGS)

    def start(self) -> None:
        me = self.table.read(os.getpid())
        info = {"pid": me.pid, "start": me.start, "unit": self.table.unit, "step": self.step, "port": self.port,
                "cwd": str(self.cwd), "argv": self.argv, "t": time.time()}
        pi.write_atomic(str(self.run_dir / pi.LAUNCHER_NAME), json.dumps(info, indent=1) + "\n")
        pi.write_atomic(str(self.run_dir / pi.PIDS_NAME), pi.pids_header(self.table.unit))
        out = open(self.log_path, "ab") if self.log_path else subprocess.DEVNULL
        old = signal.pthread_sigmask(signal.SIG_BLOCK, _SIGS)
        try:
            self.proc = subprocess.Popen(self.argv, cwd=str(self.cwd), env=self.env, start_new_session=True,
                                         preexec_fn=self._child, stdin=subprocess.DEVNULL, stdout=out,
                                         stderr=subprocess.STDOUT, close_fds=True)
            # identity is settled before any pending signal can interrupt us
            try:
                with open(self.run_dir / pi.ID_NAME, encoding="utf-8") as f:
                    cid = pi.parse_id(f.read())
                if cid.master_pid != self.proc.pid or cid.unit != self.table.unit:
                    raise pi.FormatError(f"candidate.id pid {cid.master_pid} != Popen pid {self.proc.pid}")
                self.cid = cid
            except (OSError, pi.FormatError) as exc:
                # the master is our unreaped child, so its pid cannot be reused yet
                self.id_problem = str(exc)
                p = self.table.read(self.proc.pid)
                if p is not None and p.pgid == p.sid == p.pid:
                    self.cid = pi.CandidateId(self.table.unit, p.pgid, p.sid, p.pid, p.start)
        finally:
            signal.pthread_sigmask(signal.SIG_SETMASK, old)
            if out is not subprocess.DEVNULL:
                out.close()

    def reap(self) -> None:
        if self.proc is not None:
            self.proc.poll()


# ------------------------------------------------------------------------- monitor


class Monitor(threading.Thread):
    """Records every candidate process in candidate.pids, samples worker RSS and
    available memory, and fires the abort callback once when an abort line is crossed."""

    def __init__(self, table, cid: pi.CandidateId, run_dir: Path, *, abort_rss_mb: float, abort_avail_mb: float,
                 on_abort, fast: float = 0.05, slow: float = 0.2):
        super().__init__(daemon=True)
        self.table, self.cid, self.run_dir = table, cid, Path(run_dir)
        self.abort_rss_mb, self.abort_avail_mb, self.on_abort = abort_rss_mb, abort_avail_mb, on_abort
        self.fast, self.slow = fast, slow
        self.fast_mode = True
        self._lock = threading.Lock()
        self._halt = threading.Event()
        self.recorded: list[pi.Proc] = []
        self._seen: set[tuple[int, int]] = set()
        self.samples: list[dict] = []
        self.peak_worker_mb = 0.0
        self.min_available_mb: float | None = None
        self.abort: dict | None = None
        self.errors: list[dict] = []
        self.undetermined_seen: list[dict] = []

    def recorded_list(self) -> list[pi.Proc]:
        with self._lock:
            return list(self.recorded)

    def workers(self, res: pi.Resolution) -> list[pi.Proc]:
        return [p for p in res.targets.values() if p.ppid == self.cid.master_pid]

    def record(self, res: pi.Resolution) -> None:
        new = [p for p in res.targets.values() if (p.pid, p.start) not in self._seen]
        if not new:
            return
        with open(self.run_dir / pi.PIDS_NAME, "a", encoding="utf-8") as f:
            for p in new:
                f.write(pi.pids_line(p))
            f.flush()
            os.fsync(f.fileno())
        with self._lock:
            for p in new:
                self._seen.add((p.pid, p.start))
                self.recorded.append(p)

    def trigger_abort(self, info: dict) -> None:
        with self._lock:
            if self.abort is not None:
                return
            self.abort = info | {"t": time.time()}
        self.on_abort(self.abort)

    def tick(self) -> None:
        res = pi.scan_resolve(self.table, self.cid, self.recorded_list())
        self.record(res)
        rss = {}
        for p in self.workers(res):
            try:
                rss[p.pid] = psutil.Process(p.pid).memory_info().rss / MB
            except psutil.Error:
                pass
        avail = psutil.virtual_memory().available / MB
        top = max(rss.values(), default=0.0)
        self.peak_worker_mb = max(self.peak_worker_mb, top)
        self.min_available_mb = avail if self.min_available_mb is None else min(self.min_available_mb, avail)
        self.samples.append({"t": time.time(), "worker_rss_mb": rss, "available_mb": avail,
                             "targets": len(res.targets), "state": res.state})
        if res.state == "undetermined":
            self.undetermined_seen.append({"t": time.time(), "pids": sorted(res.undetermined),
                                           "unreadable": [u["pid"] for u in res.unreadable]})
        if top > self.abort_rss_mb:
            self.trigger_abort({"reason": "worker_rss", "worker_rss_mb": top, "limit_mb": self.abort_rss_mb})
        elif avail < self.abort_avail_mb:
            self.trigger_abort({"reason": "available_memory", "available_mb": avail, "limit_mb": self.abort_avail_mb})

    def run(self):
        while not self._halt.is_set():
            try:
                self.tick()
            except Exception as exc:  # recorded; the watchdog keeps running
                if len(self.errors) < 20:
                    self.errors.append({"t": time.time(), "type": type(exc).__name__, "message": str(exc)})
            self._halt.wait(self.fast if self.fast_mode else self.slow)

    def halt(self, timeout: float = 10.0) -> None:
        self._halt.set()
        if self.is_alive() and threading.current_thread() is not self:
            self.join(timeout)

    def summary(self) -> dict:
        return {"peak_worker_rss_mb": self.peak_worker_mb, "min_available_mb": self.min_available_mb,
                "available_source": AVAILABLE_SOURCE, "abort": self.abort, "n_samples": len(self.samples),
                "errors": self.errors, "undetermined_seen": self.undetermined_seen[:20],
                "recorded": [vars(p) for p in self.recorded_list()]}


# -------------------------------------------------------------------------- run


def wait_ready(launcher: Launcher, port: int, expected_workers: int, timeout: float) -> str | None:
    """None when ready, otherwise the reason."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if launcher.proc.poll() is not None:
            return f"master exited with {launcher.proc.returncode} before ready"
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=1.0):
                pass
            if len(psutil.Process(launcher.proc.pid).children()) >= expected_workers:
                return None
        except (OSError, psutil.Error):
            pass
        time.sleep(0.1)
    return f"not ready within {timeout} s"


def _raise_interrupted(signum, frame):
    for s in _SIGS:
        signal.signal(s, signal.SIG_IGN)
    raise Interrupted(signum)


def run_step(D: Path, step: str, port: int, *, observe_mb: float, abort_rss_mb: float, abort_avail_mb: float,
             python: str | None = None, gunicorn: str | None = None, argv: list[str] | None = None,
             scenario_fn=None, expected_workers: int = 1, startup_timeout: float = 120.0,
             step_deadline: float = 1200.0, http_timeout: float = load.DEFAULT_TIMEOUT,
             term_wait: float = 30.0, abort_term_wait: float = 10.0, table=None,
             install_handlers: bool = True) -> int:
    """One step on a fresh candidate.  ``argv``/``scenario_fn`` replace gunicorn and
    the load scenario (tests).  Returns the exit status (module docstring)."""
    D = Path(D)
    run_dir = D / "run" / step
    if run_dir.exists():
        print(json.dumps({"step": step, "error": f"{run_dir} exists; a step directory is never reused"}))
        return EXIT_FORMAT
    for sub in SUBDIRS:
        (D / sub).mkdir(parents=True, exist_ok=True)
    run_dir.mkdir(parents=True)
    table = table or pi.platform_table()
    if argv is None:
        argv = gunicorn_argv(python, gunicorn, port, D, run_dir, step)
    env = candidate_env(D, python or sys.executable)
    launcher = Launcher(argv, D / "app" / "API", env, run_dir, step=step, port=port, table=table,
                        log_path=D / "logs" / f"candidate_{step}.out")
    doc: dict = {"step": step, "scenario": STEPS.get(step), "port": port, "argv": argv, "env": env,
                 "thresholds": {"observe_worker_rss_mb": observe_mb, "abort_worker_rss_mb": abort_rss_mb,
                                "abort_min_available_mb": abort_avail_mb},
                 "started": time.time()}
    codes: list[int] = []
    stops: list[dict] = []
    stop_lock = threading.Lock()
    monitor: Monitor | None = None
    timer: threading.Timer | None = None
    old_handlers = {s: signal.getsignal(s) for s in _SIGS} if install_handlers else {}

    def do_stop(wait: float) -> pi.StopResult | None:
        with stop_lock:
            if launcher.cid is None:
                return None
            rec = monitor.recorded_list if monitor else (lambda: [])
            res = pi.stop_candidate(table, launcher.cid, rec, term_wait=wait, port=port, reap=launcher.reap)
            stops.append(res.as_dict())
            return res

    if install_handlers:
        for s in _SIGS:
            signal.signal(s, _raise_interrupted)
    final: pi.StopResult | None = None
    try:
        launcher.start()
        doc["candidate_id"] = vars(launcher.cid) if launcher.cid else None
        if launcher.id_problem:
            doc["id_problem"] = launcher.id_problem
            codes.append(EXIT_FORMAT)
        if launcher.cid is None:
            raise RuntimeError("candidate identity could not be established")
        monitor = Monitor(table, launcher.cid, run_dir, abort_rss_mb=abort_rss_mb, abort_avail_mb=abort_avail_mb,
                          on_abort=lambda info: do_stop(abort_term_wait))
        monitor.start()
        timer = threading.Timer(step_deadline, monitor.trigger_abort, args=({"reason": "step_deadline",
                                                                             "limit_s": step_deadline},))
        timer.daemon = True
        timer.start()
        problem = wait_ready(launcher, port, expected_workers, startup_timeout)
        monitor.fast_mode = False
        if problem:
            doc["startup_problem"] = problem
            codes.append(EXIT_FAILED)
        else:
            fn = scenario_fn
            if fn is None:
                fn = load.SCENARIOS[STEPS[step]]
                if step == "smoke":
                    fn = functools.partial(load.scenario_smoke_J01, out_dir=D / "out")
            rc = load.run_scenario(STEPS.get(step, step), f"http://127.0.0.1:{port}", launcher.proc.pid, D / "out",
                                   timeout=http_timeout, scenario_fn=fn, expected_workers=expected_workers)
            doc["load_exit"] = rc
            if rc != 0:
                codes.append(EXIT_FAILED)
    except Interrupted as exc:
        doc["interrupted"] = str(exc)
        codes.append(EXIT_FAILED)
    except pi.FormatError as exc:
        doc["error"] = {"type": "FormatError", "message": str(exc)}
        codes.append(EXIT_FORMAT)
    except Exception as exc:
        doc["error"] = {"type": type(exc).__name__, "message": str(exc)}
        codes.append(EXIT_FAILED)
    finally:
        if install_handlers:
            for s in _SIGS:
                signal.signal(s, signal.SIG_IGN)  # cleanup must not be interrupted
        if timer:
            timer.cancel()
        try:
            final = do_stop(term_wait)
        except Exception as exc:  # the result is still written; nothing is reported clean
            doc["stop_error"] = {"type": type(exc).__name__, "message": str(exc)}
            codes.append(EXIT_UNDETERMINED)
        if monitor:
            monitor.halt()
        launcher.reap()
        if final is None and launcher.proc is not None and launcher.proc.poll() is None:
            # no identity at all, but the master is our unreaped child: stop it through Popen
            launcher.proc.terminate()
            try:
                launcher.proc.wait(term_wait)
            except subprocess.TimeoutExpired:
                launcher.proc.kill()
                launcher.proc.wait(10)
            codes.append(EXIT_UNDETERMINED)
            doc["no_identity_stop"] = "master stopped via Popen; workers/helpers could not be checked"
        if install_handlers:
            for s, h in old_handlers.items():
                signal.signal(s, h)

    if monitor:
        doc["monitor"] = monitor.summary()
        doc["monitor_series"] = monitor.samples
        if monitor.abort:
            codes.append(EXIT_ABORT)
        elif monitor.peak_worker_mb > observe_mb:
            codes.append(EXIT_OBSERVE)
            doc["observation_line_crossed"] = True
    doc["stops"] = stops
    if final is not None:
        doc["exit_check"] = {"state": final.state, "port_free": final.port_free, "mismatch": final.mismatch}
        if final.state == "undetermined":
            codes.append(EXIT_UNDETERMINED)
        elif final.state == "residual":
            codes.append(EXIT_RESIDUAL)
        if final.mismatch:
            codes.append(EXIT_MISMATCH)
    doc["exit_status"] = worst(codes)
    doc["ended"] = time.time()
    text = json.dumps(doc, indent=1, default=str)
    (run_dir / "result.json").write_text(text)
    (D / "out" / f"candidate_{step}.json").write_text(text)
    print(json.dumps({"step": step, "exit_status": doc["exit_status"],
                      "exit_check": doc.get("exit_check"), "load_exit": doc.get("load_exit"),
                      "peak_worker_rss_mb": round(monitor.peak_worker_mb) if monitor else None,
                      "min_available_mb": round(monitor.min_available_mb) if monitor and monitor.min_available_mb else None,
                      "abort": monitor.abort if monitor else None, "interrupted": doc.get("interrupted"),
                      "error": doc.get("error"), "startup_problem": doc.get("startup_problem")}, default=str))
    return doc["exit_status"]


# ------------------------------------------------------------------------ recover


def _recover_one(rd: Path, mode: str, table, app: str, entry: dict, signaller, term_wait, kill_wait, poll) -> list[int]:
    codes: list[int] = []
    if not (rd / pi.ID_NAME).exists():
        # by construction gunicorn never ran here; list processes whose cwd is the app dir
        snap = table.scan()
        hits = [pid for pid in snap if (table.cwd(pid) and os.path.realpath(table.cwd(pid)) == app)]
        entry["no_candidate_id"] = True
        entry["cwd_matches"] = hits
        entry["unreadable"] = list(table.bad)
        entry["state"] = "undetermined" if hits or table.bad else "none"
        return [EXIT_UNDETERMINED] if entry["state"] == "undetermined" else []
    try:
        cid, recorded = pi.read_state(str(rd), table.unit)
    except (pi.FormatError, OSError) as exc:
        entry["state"] = "format_error"
        entry["error"] = str(exc)
        return [EXIT_FORMAT]
    if mode == "list":
        res = pi.scan_resolve(table, cid, recorded)
        entry.update(res.as_dict())
        codes.append({"undetermined": EXIT_UNDETERMINED, "alive": EXIT_FAILED}.get(res.state, EXIT_OK))
        if res.mismatch:
            codes.append(EXIT_MISMATCH)
    else:
        sr = pi.stop_candidate(table, cid, lambda: recorded, term_wait=term_wait, kill_wait=kill_wait,
                               poll=poll, port=None, signaller=signaller)
        entry.update(sr.as_dict())
        codes.append({"undetermined": EXIT_UNDETERMINED, "residual": EXIT_RESIDUAL}.get(sr.state, EXIT_OK))
        if sr.mismatch:
            codes.append(EXIT_MISMATCH)
    return codes


def recover(D: Path, mode: str, port: int | None, *, table=None, signaller=None,
            term_wait: float = 30.0, kill_wait: float = 5.0, poll: float = 0.2) -> int:
    """Identify (list) or stop (stop) the candidates of every step directory in D/run."""
    D = Path(D)
    table = table or pi.platform_table()
    codes: list[int] = []
    report: list[dict] = []
    run_root = D / "run"
    dirs = sorted(d for d in run_root.iterdir() if d.is_dir()) if run_root.is_dir() else []
    app = os.path.realpath(D / "app" / "API")
    for rd in dirs:
        entry: dict = {"run_dir": str(rd)}
        report.append(entry)
        try:
            codes.extend(_recover_one(rd, mode, table, app, entry, signaller, term_wait, kill_wait, poll))
        except Exception as exc:  # never report a step clean when it could not be checked
            entry["state"] = "error"
            entry["error"] = {"type": type(exc).__name__, "message": str(exc)}
            codes.append(EXIT_UNDETERMINED)
    pf = pi.port_free(port) if port else None
    # a bound port only counts when no candidate is known to be alive (it is
    # expected then); same order as candidate_recover.sh
    if pf is False and not ({EXIT_FAILED, EXIT_RESIDUAL, EXIT_UNDETERMINED} & set(codes)):
        codes.append(EXIT_RESIDUAL)
    status = worst(codes)
    msg = {EXIT_OK: "no candidate process", EXIT_FAILED: "CANDIDATE PROCESSES ALIVE",
           EXIT_FORMAT: "STATE FILE FORMAT ERROR", EXIT_RESIDUAL: "CANDIDATE STILL ALIVE OR PORT STILL BOUND",
           EXIT_UNDETERMINED: "IDENTITY UNDETERMINED: nothing signalled for those processes",
           EXIT_MISMATCH: "no candidate process; a recorded number now belongs to another process (not signalled)"}
    print(json.dumps({"mode": mode, "exit_status": status, "message": msg.get(status, ""), "port_free": pf,
                      "runs": report}, indent=1, default=str))
    return status


# ---------------------------------------------------------------------------- cli


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="start a fresh candidate, measure one step, stop and check")
    r.add_argument("--dir", required=True)
    r.add_argument("--step", choices=STEPS, required=True)
    r.add_argument("--port", type=int, required=True)
    r.add_argument("--python", required=True, help="interpreter of the environment to test (run with -B)")
    r.add_argument("--gunicorn", required=True, help="gunicorn script of that environment")
    r.add_argument("--observe-worker-rss-mb", type=float, required=True)
    r.add_argument("--abort-worker-rss-mb", type=float, required=True)
    r.add_argument("--abort-min-available-mb", type=float, required=True)
    r.add_argument("--startup-timeout", type=float, default=120.0)
    r.add_argument("--step-deadline", type=float, default=1200.0)
    r.add_argument("--timeout", type=float, default=load.DEFAULT_TIMEOUT, help="per-request HTTP timeout (s)")
    c = sub.add_parser("recover", help="identify or stop candidates left by a runner that died")
    c.add_argument("--dir", required=True)
    c.add_argument("--mode", choices=("list", "stop"), required=True)
    c.add_argument("--port", type=int, required=True)
    a = ap.parse_args(argv)
    if a.cmd == "run":
        return run_step(Path(a.dir), a.step, a.port, observe_mb=a.observe_worker_rss_mb,
                        abort_rss_mb=a.abort_worker_rss_mb, abort_avail_mb=a.abort_min_available_mb,
                        python=a.python, gunicorn=a.gunicorn, startup_timeout=a.startup_timeout,
                        step_deadline=a.step_deadline, http_timeout=a.timeout)
    return recover(Path(a.dir), a.mode, a.port)


if __name__ == "__main__":
    sys.exit(main())
