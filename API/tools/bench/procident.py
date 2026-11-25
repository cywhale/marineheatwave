"""Process identity for isolated candidate runs (doc 24 r4 §1.3-§1.5).

Benchmark tooling only; never imported by the API.  The same rules are
implemented for Linux in ``candidate_recover.sh`` (shell fallback).

Start values
------------
``linux-ticks``  field 22 (starttime) of ``/proc/PID/stat``, written verbatim
                 (clock ticks since boot).  Read directly from /proc, never via
                 psutil, so Python and the shell fallback compare the same string.
``darwin-usec``  round(psutil create_time * 1e6); macOS, local tests only.
A file whose unit differs from the reader's table is a format error.

Identity (``resolve``)
----------------------
V  = recorded processes (master from candidate.id + candidate.pids) that are
     alive now with the recorded start.
Anchored sessions / groups = the current sid / pgid of every process in V.  A
process's sid only changes through its own setsid() (to its own pid), groups
can only be joined within a session and sessions are only inherited, so every
member of an anchored session or group descends from a candidate process.
T  = V + members of anchored sessions/groups + (recursively) children of T.
F  = "foreign": not in T, its sid/pgid equals a recorded number whose leader
     (pid == that number) is alive with a different start than recorded.  While
     any of our processes used that number the kernel could not have handed it
     out again, so these processes are provably not ours.
U  = "undetermined": not in T, sid/pgid equals a recorded number, not provably
     foreign.  Any U blocks every signal (state "undetermined").

Signals are sent one process at a time after re-checking its start; on Linux
through a pidfd opened before the check, so a reused PID cannot be hit.  No
killpg.
"""
from __future__ import annotations

import errno
import os
import re
import signal
import socket
import sys
import time
from dataclasses import dataclass, field
from typing import Callable

UNIT_LINUX = "linux-ticks"
UNIT_DARWIN = "darwin-usec"
UNITS = (UNIT_LINUX, UNIT_DARWIN)
ID_NAME = "candidate.id"
PIDS_NAME = "candidate.pids"
LAUNCHER_NAME = "launcher.json"
_NUM = r"(0|[1-9][0-9]*)"
_ID_RE = re.compile(rf"^mhw-candidate-id v1 ({UNIT_LINUX}|{UNIT_DARWIN}) {_NUM} {_NUM} {_NUM} {_NUM}$")
_PIDS_HEADER_RE = re.compile(rf"^# mhw-candidate-pids v1 ({UNIT_LINUX}|{UNIT_DARWIN})$")
_PIDS_LINE_RE = re.compile(rf"^{_NUM} {_NUM} {_NUM} {_NUM} {_NUM}$")


class FormatError(ValueError):
    """A state file or /proc entry that does not match the expected format."""


@dataclass(frozen=True)
class Proc:
    pid: int
    ppid: int
    pgid: int
    sid: int
    start: int


# ------------------------------------------------------------------ process tables


def parse_stat(pid: int, text: str) -> tuple[str, Proc]:
    """Parse one /proc/PID/stat.  comm may contain spaces and parentheses, so
    the fields are split after the LAST ')'.  Returns (state, Proc)."""
    head = f"{pid} ("
    if not text.startswith(head):
        raise FormatError(f"stat of {pid} does not start with {head!r}")
    i = text.rfind(")")
    if i < len(head) - 1:
        raise FormatError(f"stat of {pid}: no closing ')'")
    rest = text[i + 1:].split()
    if len(rest) < 20:
        raise FormatError(f"stat of {pid}: {len(rest)} fields after comm, need >= 20")
    state = rest[0]
    vals = [rest[1], rest[2], rest[3], rest[19]]  # ppid, pgrp, session, starttime (field 22)
    if not all(v.isdigit() for v in vals):
        raise FormatError(f"stat of {pid}: non-numeric field in {vals}")
    ppid, pgid, sid, start = (int(v) for v in vals)
    return state, Proc(pid, ppid, pgid, sid, start)


class LinuxProcTable:
    """Reads /proc (or a fixture tree with the same layout)."""

    unit = UNIT_LINUX

    def __init__(self, root: str = "/proc"):
        self.root = root
        self.bad: list[str] = []  # stat entries that could not be parsed in the last scan

    def read(self, pid: int) -> Proc | None:
        """The live, non-zombie process with this pid, or None."""
        try:
            with open(f"{self.root}/{pid}/stat", encoding="utf-8", errors="surrogateescape") as f:
                text = f.read()
        except (FileNotFoundError, ProcessLookupError, NotADirectoryError):
            return None
        state, proc = parse_stat(pid, text)
        return None if state in ("Z", "X", "x") else proc

    def scan(self) -> dict[int, Proc]:
        out: dict[int, Proc] = {}
        self.bad = []
        for name in os.listdir(self.root):
            if not name.isdigit():
                continue
            try:
                p = self.read(int(name))
            except FormatError as exc:
                self.bad.append(str(exc))
                continue
            except OSError:
                continue
            if p is not None:
                out[p.pid] = p
        return out

    def cwd(self, pid: int) -> str | None:
        try:
            return os.readlink(f"{self.root}/{pid}/cwd")
        except OSError:
            return None


class DarwinProcTable:
    """macOS via psutil; local tests only (start = create_time in integer microseconds)."""

    unit = UNIT_DARWIN

    def __init__(self):
        import psutil  # local tests only
        self._psutil = psutil
        self.bad: list[str] = []

    def read(self, pid: int) -> Proc | None:
        ps = self._psutil
        try:
            p = ps.Process(pid)
            if p.status() == ps.STATUS_ZOMBIE:
                return None
            return Proc(pid, p.ppid(), os.getpgid(pid), os.getsid(pid), round(p.create_time() * 1e6))
        except (ps.NoSuchProcess, ps.ZombieProcess, ProcessLookupError):
            return None
        except (ps.AccessDenied, PermissionError, OSError):
            return None

    def scan(self) -> dict[int, Proc]:
        out = {}
        for pid in self._psutil.pids():
            p = self.read(pid)
            if p is not None:
                out[pid] = p
        return out

    def cwd(self, pid: int) -> str | None:
        try:
            return self._psutil.Process(pid).cwd()
        except Exception:
            return None


def platform_table():
    if sys.platform.startswith("linux"):
        return LinuxProcTable()
    if sys.platform == "darwin":
        return DarwinProcTable()
    raise RuntimeError(f"unsupported platform {sys.platform}")


# ---------------------------------------------------------------------- state files


@dataclass(frozen=True)
class CandidateId:
    unit: str
    pgid: int
    sid: int
    master_pid: int
    master_start: int

    def line(self) -> str:
        return f"mhw-candidate-id v1 {self.unit} {self.pgid} {self.sid} {self.master_pid} {self.master_start}\n"


def parse_id(text: str) -> CandidateId:
    body = text[:-1] if text.endswith("\n") else text
    if "\n" in body:
        raise FormatError("candidate.id must be exactly one line")
    m = _ID_RE.match(body)
    if not m:
        raise FormatError(f"candidate.id has an unexpected format: {body[:120]!r}")
    unit = m.group(1)
    pgid, sid, pid, start = (int(m.group(i)) for i in range(2, 6))
    if not (pgid == sid == pid) or pid <= 1:
        raise FormatError(f"candidate.id: pgid={pgid} sid={sid} pid={pid} must be equal and > 1")
    return CandidateId(unit, pgid, sid, pid, start)


def pids_header(unit: str) -> str:
    return f"# mhw-candidate-pids v1 {unit}\n"


def pids_line(p: Proc) -> str:
    return f"{p.pid} {p.start} {p.ppid} {p.pgid} {p.sid}\n"


def parse_pids(text: str, unit: str) -> list[Proc]:
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        raise FormatError("candidate.pids is empty (no header)")
    m = _PIDS_HEADER_RE.match(lines[0])
    if not m:
        raise FormatError(f"candidate.pids header: {lines[0][:120]!r}")
    if m.group(1) != unit:
        raise FormatError(f"candidate.pids unit {m.group(1)} != {unit}")
    out = []
    for n, ln in enumerate(lines[1:], start=2):
        if not _PIDS_LINE_RE.match(ln):
            raise FormatError(f"candidate.pids line {n}: {ln[:120]!r}")
        pid, start, ppid, pgid, sid = (int(x) for x in ln.split())
        out.append(Proc(pid, ppid, pgid, sid, start))
    return out


def read_state(run_dir: str, unit: str) -> tuple[CandidateId, list[Proc]]:
    """Read and validate candidate.id and candidate.pids (FormatError on any problem)."""
    try:
        with open(os.path.join(run_dir, ID_NAME), encoding="utf-8") as f:
            cid = parse_id(f.read())
    except FileNotFoundError:
        raise
    except UnicodeDecodeError as exc:
        raise FormatError(f"candidate.id is not UTF-8: {exc}") from None
    if cid.unit != unit:
        raise FormatError(f"candidate.id unit {cid.unit} != table unit {unit}")
    try:
        with open(os.path.join(run_dir, PIDS_NAME), encoding="utf-8") as f:
            recorded = parse_pids(f.read(), unit)
    except FileNotFoundError:
        raise FormatError("candidate.pids is missing") from None
    except UnicodeDecodeError as exc:
        raise FormatError(f"candidate.pids is not UTF-8: {exc}") from None
    return cid, recorded


def write_atomic(path: str, text: str) -> None:
    """os-level write + rename (also used inside preexec_fn, so no Python file objects)."""
    tmp = f"{path}.tmp.{os.getpid()}"
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        data = text.encode()
        while data:
            data = data[os.write(fd, data):]
    finally:
        os.close(fd)
    os.rename(tmp, path)


# ---------------------------------------------------------------------- resolution


@dataclass
class Resolution:
    targets: dict[int, Proc] = field(default_factory=dict)
    how: dict[int, str] = field(default_factory=dict)  # recorded | session | group | child
    foreign: dict[int, Proc] = field(default_factory=dict)
    undetermined: dict[int, Proc] = field(default_factory=dict)
    master: str = "gone"  # alive | gone | reused

    @property
    def state(self) -> str:
        if self.undetermined:
            return "undetermined"
        return "alive" if self.targets else "none"

    @property
    def mismatch(self) -> bool:
        """A recorded number/pid now provably belongs to another process."""
        return self.master == "reused" or bool(self.foreign)

    def ordered_targets(self, master_pid: int) -> list[Proc]:
        """Master first, then by pid (stable output for logs and tests)."""
        return sorted(self.targets.values(), key=lambda p: (p.pid != master_pid, p.pid))

    def as_dict(self) -> dict:
        d = lambda m: {str(k): vars(v) | ({"how": self.how[k]} if k in self.how else {}) for k, v in sorted(m.items())}
        return {"state": self.state, "master": self.master, "mismatch": self.mismatch,
                "targets": d(self.targets), "foreign": d(self.foreign), "undetermined": d(self.undetermined)}


def resolve(snapshot: dict[int, Proc], cid: CandidateId, recorded: list[Proc]) -> Resolution:
    res = Resolution()
    pairs = [(cid.master_pid, cid.master_start)] + [(r.pid, r.start) for r in recorded]
    cur = snapshot.get(cid.master_pid)
    res.master = "gone" if cur is None else ("alive" if cur.start == cid.master_start else "reused")

    for pid, start in pairs:
        p = snapshot.get(pid)
        if p is not None and p.start == start:
            res.targets[pid] = p
            res.how[pid] = "recorded"
    anchored_s = {p.sid for p in res.targets.values()}
    anchored_g = {p.pgid for p in res.targets.values()}
    for pid, p in snapshot.items():
        if pid in res.targets:
            continue
        if p.sid in anchored_s:
            res.targets[pid], res.how[pid] = p, "session"
        elif p.pgid in anchored_g:
            res.targets[pid], res.how[pid] = p, "group"
    changed = True
    while changed:
        changed = False
        for pid, p in snapshot.items():
            if pid not in res.targets and p.ppid in res.targets:
                res.targets[pid], res.how[pid] = p, "child"
                changed = True

    numbers = {cid.sid, cid.pgid}
    leader_starts: dict[int, set[int]] = {cid.sid: {cid.master_start}}
    for r in recorded:
        numbers.update((r.sid, r.pgid))
        if r.pid in (r.sid, r.pgid):
            leader_starts.setdefault(r.pid, set()).add(r.start)
    for pid, p in snapshot.items():
        if pid in res.targets:
            continue
        hits = {p.sid, p.pgid} & numbers
        if not hits:
            continue
        foreign = any(n in leader_starts and n in snapshot and snapshot[n].start not in leader_starts[n]
                      for n in hits)
        (res.foreign if foreign else res.undetermined)[pid] = p
    return res


# ---------------------------------------------------------------------- signalling


class Signaller:
    """Send a signal to one recorded process after re-checking its start.

    Linux: pidfd_open first, then the start check, then pidfd_send_signal, so
    the signal can only reach the checked process.  Without pidfd (macOS tests,
    or ENOSYS) the check is followed by os.kill and a tiny window remains;
    ``pidfd_unavailable`` records that.  ``kill`` (tests) replaces the delivery.
    """

    def __init__(self, table, kill: Callable[[int, int], None] | None = None):
        self.table = table
        self.kill = kill
        self.use_pidfd = kill is None and isinstance(table, LinuxProcTable) and table.root == "/proc" \
            and hasattr(os, "pidfd_open") and hasattr(signal, "pidfd_send_signal")
        self.pidfd_unavailable = not self.use_pidfd and kill is None
        self.log: list[dict] = []

    def _record(self, p: Proc, sig: int, outcome: str) -> str:
        self.log.append({"t": time.time(), "pid": p.pid, "start": p.start, "signal": signal.Signals(sig).name,
                         "outcome": outcome})
        return outcome

    def send(self, p: Proc, sig: int) -> str:
        fd = None
        if self.use_pidfd:
            try:
                fd = os.pidfd_open(p.pid)
            except ProcessLookupError:
                return self._record(p, sig, "gone")
            except OSError as exc:
                if exc.errno != errno.ENOSYS:
                    return self._record(p, sig, f"error:{exc}")
                self.use_pidfd, self.pidfd_unavailable = False, True
        try:
            try:
                cur = self.table.read(p.pid)
            except FormatError as exc:
                return self._record(p, sig, f"unreadable:{exc}")
            if cur is None:
                return self._record(p, sig, "gone")
            if cur.start != p.start:
                return self._record(p, sig, "changed")
            try:
                if fd is not None:
                    signal.pidfd_send_signal(fd, sig)
                elif self.kill is not None:
                    self.kill(p.pid, sig)
                else:
                    os.kill(p.pid, sig)
            except ProcessLookupError:
                return self._record(p, sig, "gone")
            except OSError as exc:
                return self._record(p, sig, f"error:{exc}")
            return self._record(p, sig, "sent")
        finally:
            if fd is not None:
                os.close(fd)


# ------------------------------------------------------------------ stop and check


def port_free(port: int, host: str = "127.0.0.1") -> bool:
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        s.bind((host, port))
        return True
    except OSError:
        return False
    finally:
        s.close()


@dataclass
class StopResult:
    state: str  # clean | residual | undetermined
    final: Resolution
    port_free: bool | None
    signals: list[dict]
    rounds: list[dict]
    pidfd_unavailable: bool

    @property
    def mismatch(self) -> bool:
        return self.final.mismatch

    def as_dict(self) -> dict:
        return {"state": self.state, "port_free": self.port_free, "mismatch": self.mismatch,
                "pidfd_unavailable": self.pidfd_unavailable, "signals": self.signals,
                "rounds": self.rounds, "final": self.final.as_dict()}


def stop_candidate(table, cid: CandidateId, recorded: Callable[[], list[Proc]], *,
                   term_wait: float = 30.0, kill_wait: float = 5.0, poll: float = 0.2,
                   port: int | None = None, reap: Callable[[], None] | None = None,
                   signaller: Signaller | None = None) -> StopResult:
    """TERM -> wait -> KILL -> wait, re-resolving identity before every signal round.

    Any undetermined process stops the procedure at once (no signal at all in
    that round or later); the result is never reported as clean then.  Every
    process proven to be a target in one round is carried into later rounds as
    recorded (pid + start), so stopping the master first cannot turn its
    unrecorded workers into undetermined processes.
    """
    sg = signaller or Signaller(table)
    rounds: list[dict] = []
    proven: dict[tuple[int, int], Proc] = {}

    def check() -> Resolution:
        if reap:
            reap()
        r = resolve(table.scan(), cid, recorded() + list(proven.values()))
        for p in r.targets.values():
            proven[(p.pid, p.start)] = p
        rounds.append({"t": time.time(), "state": r.state, "targets": sorted(r.targets),
                       "undetermined": sorted(r.undetermined), "foreign": sorted(r.foreign)})
        return r

    def wait_gone(limit: float) -> Resolution:
        deadline = time.monotonic() + limit
        r = check()
        while r.state == "alive" and time.monotonic() < deadline:
            time.sleep(poll)
            r = check()
        return r

    def finish(r: Resolution) -> StopResult:
        if r.state == "undetermined":
            state, pf = "undetermined", None
        else:
            pf = port_free(port) if port else None
            state = "residual" if r.state == "alive" or pf is False else "clean"
        return StopResult(state, r, pf, sg.log, rounds, sg.pidfd_unavailable)

    r = check()
    for sig, limit in ((signal.SIGTERM, term_wait), (signal.SIGKILL, kill_wait)):
        if r.state != "alive":
            break
        for p in r.ordered_targets(cid.master_pid):
            sg.send(p, sig)
        r = wait_gone(limit)
    return finish(r)
