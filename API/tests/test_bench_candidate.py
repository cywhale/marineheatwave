"""Isolated candidate runner: identity, stop and recovery (doc 24 r4 §1.6, V1-V17).

Fixture tests build a fake /proc tree (Linux layout) and run the same cases
through procident (Python) and candidate_recover.sh (shell), with the kill
command replaced.  Live tests start small fake process trees (master -> worker
-> helper, no MHW data) on this machine.  Linux-only behaviour (pidfd,
PR_SET_PDEATHSIG, the real /proc) is not exercised here (doc 24 r4 L1).
"""
from __future__ import annotations

import json
import os
import signal
import socket
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import psutil
import pytest

from tools.bench import candidate, load
from tools.bench import procident as pi

API = Path(__file__).resolve().parents[1]
LINUX = sys.platform.startswith("linux")
linux_only = pytest.mark.skipif(not LINUX, reason="Linux-only path (doc 24 L1)")
SHELL = API / "tools" / "bench" / "candidate_recover.sh"
pytestmark = pytest.mark.skipif(sys.platform not in ("darwin",) and not sys.platform.startswith("linux"),
                                reason="POSIX only")


# =========================================================== fixture /proc tree


def stat_line(pid, ppid, pgid, sid, start, state="S", comm="python3"):
    rest = [state, str(ppid), str(pgid), str(sid)] + ["0"] * 15 + [str(start)] + ["0"] * 5
    return f"{pid} ({comm}) " + " ".join(rest) + "\n"


class FakeProc:
    """A fixture /proc tree plus the state files of one step directory."""

    def __init__(self, base: Path):
        self.root = base / "proc"
        self.root.mkdir()
        self.D = base / "D"
        self.run = self.D / "run" / "s1"
        self.run.mkdir(parents=True)
        (self.D / "app" / "API").mkdir(parents=True)
        self.killlog = base / "kill.log"

    def add(self, pid, ppid, pgid, sid, start, state="S", comm="python3"):
        d = self.root / str(pid)
        d.mkdir(exist_ok=True)
        (d / "stat").write_text(stat_line(pid, ppid, pgid, sid, start, state, comm))

    def remove(self, pid):
        d = self.root / str(pid)
        for f in d.iterdir():
            f.unlink()
        d.rmdir()

    def set_id(self, text):
        (self.run / "candidate.id").write_text(text)

    def set_master(self, pid, start):
        self.set_id(f"mhw-candidate-id v1 linux-ticks {pid} {pid} {pid} {start}\n")

    def set_pids(self, rows, header="# mhw-candidate-pids v1 linux-ticks"):
        body = header + "\n" + "".join(f"{pid} {start} {ppid} {pgid} {sid}\n" for pid, start, ppid, pgid, sid in rows)
        (self.run / "candidate.pids").write_text(body)

    def table(self, euid=None):
        return pi.LinuxProcTable(str(self.root), euid=euid)

    def py_kill(self, removes: bool):
        def kill(pid, sig):
            with open(self.killlog, "a") as f:
                f.write(f"{signal.Signals(sig).name.replace('SIG', '')} {pid}\n")
            if removes and (self.root / str(pid)).exists():
                self.remove(pid)
        return kill

    def kills(self):
        if not self.killlog.exists():
            return []
        return [tuple(l.split()) for l in self.killlog.read_text().splitlines()]

    # ---- runners
    def python(self, mode, removes=True, euid=None):
        table = self.table(euid)
        sg = pi.Signaller(table, kill=self.py_kill(removes))
        return candidate.recover(self.D, mode, None, table=table, signaller=sg, term_wait=0, kill_wait=0, poll=0.01)

    def shell(self, mode, removes=True, port_out="", extra="", euid=None, statcmd=None):
        stub = self.killlog.parent / "killstub.sh"
        gone = self.killlog.parent / "gone"
        gone.mkdir(exist_ok=True)
        # "exit" of a fixture process = move its directory out of the fixture /proc
        stub.write_text(textwrap.dedent(f"""\
            #!/bin/bash
            # called as: killstub.sh -s SIG PID
            echo "$2 $3" >> {self.killlog}
            {extra}
            {'mv ' + str(self.root) + '/"$3" ' + str(gone) + '/"$3.$RANDOM$RANDOM"' if removes else 'true'}
            """))
        stub.chmod(0o755)
        portstub = self.killlog.parent / "portstub.sh"
        portstub.write_text(f"#!/bin/bash\nprintf '%s' '{port_out}'\n")
        portstub.chmod(0o755)
        env = {"PATH": os.environ["PATH"], "MHW_RECOVER_PROC": str(self.root), "MHW_RECOVER_KILL": str(stub),
               "MHW_RECOVER_PORTCHECK": str(portstub), "MHW_RECOVER_TERM_WAIT": "0",
               "MHW_RECOVER_KILL_WAIT": "0", "MHW_RECOVER_POLL": "1", "MHW_RECOVER_TEST_MODE": "1"}
        if euid is not None:
            env["MHW_RECOVER_EUID"] = str(euid)
        if statcmd is not None:
            env["MHW_RECOVER_STAT"] = str(statcmd)
        p = subprocess.run(["bash", str(SHELL), mode, str(self.run), "18031"], env=env, capture_output=True,
                           text=True, timeout=60)
        return p.returncode, p.stdout


def shell_sets(out):
    sets = {"candidate": set(), "foreign": set(), "UNDETERMINED": set()}
    for line in out.splitlines():
        parts = line.split()
        if parts and parts[0] in sets:
            sets[parts[0]].add(int(parts[1]))
    return sets


def py_sets(fp):
    cid, rec = pi.read_state(str(fp.run), pi.UNIT_LINUX)
    r = pi.scan_resolve(fp.table(), cid, rec)
    return {"candidate": set(r.targets), "foreign": set(r.foreign), "UNDETERMINED": set(r.undetermined)}, r


@pytest.fixture()
def fp(tmp_path):
    return FakeProc(tmp_path)


def base_tree(fp, *, master_start=1000):
    fp.add(100, 1, 100, 100, master_start, comm="python3.11")
    fp.add(101, 100, 100, 100, 1001, comm="gunicorn: worker [mhw_app:app]")
    fp.add(102, 101, 100, 100, 1002, comm="python3 (resource tracker) x")
    fp.add(500, 1, 500, 500, 50, comm="sshd")  # unrelated


# ----------------------------------------------------------------- V14 parsing


@pytest.mark.parametrize("comm", ["python3", "gunicorn: worker [mhw_app:app]", "a) b (c", "x)", "((", "sp ace) S 1 2 3"])
def test_v14_parse_stat_comm_with_spaces_and_parentheses(comm):
    state, p = pi.parse_stat(42, stat_line(42, 7, 8, 9, 123456, comm=comm))
    assert (state, p) == ("S", pi.Proc(42, 7, 8, 9, 123456))


@pytest.mark.parametrize("text", ["42 python3 S 1 2 3\n", "43 (x) S 1 2 3 4\n", "42 (x) S a 2 3 " + "0 " * 16 + "\n",
                                  "42 (x) S 1 2 3 " + "0 " * 15 + "-5 0\n"])
def test_v14_parse_stat_rejects_bad_lines(text):
    with pytest.raises(pi.FormatError):
        pi.parse_stat(42, text)


def test_v14_shell_and_python_read_the_same_start_field(fp):
    fp.add(100, 1, 100, 100, 98765, comm="odd ) name ( with ) parens")
    fp.add(101, 100, 100, 100, 98766, comm="w")
    fp.set_master(100, 98765)
    fp.set_pids([])
    sets, r = py_sets(fp)
    assert sets["candidate"] == {100, 101} and r.master == "alive"
    rc, out = fp.shell("list")
    assert rc == 1 and shell_sets(out)["candidate"] == {100, 101}


# ---------------------------------------------------------------- V15 formats

BAD_IDS = ["", "\n", "mhw-candidate-id v2 linux-ticks 100 100 100 1000\n",
           "mhw-candidate-id v1 linux-ticks 100 100 101 1000\n", "mhw-candidate-id v1 linux-ticks 100 100 100\n",
           "mhw-candidate-id v1 linux-ticks 100 100 100 1000 7\n", "mhw-candidate-id v1 linux-ticks 0100 0100 0100 1\n",
           "mhw-candidate-id v1 linux-ticks 1 1 1 5\n", "mhw-candidate-id v1 linux-ticks 100 100 100 x\n",
           "mhw-candidate-id v1 linux-ticks 100 100 100 1000\nextra\n", "mhw-candidate-id v1 seconds 100 100 100 1000\n"]


@pytest.mark.parametrize("text", BAD_IDS)
def test_v15_bad_candidate_id_is_a_format_error_in_python_and_shell(fp, text):
    base_tree(fp)
    fp.set_id(text)
    fp.set_pids([])
    with pytest.raises(pi.FormatError):
        pi.parse_id(text)
    for mode in ("list", "stop"):
        assert fp.python(mode) == candidate.EXIT_FORMAT
        rc, out = fp.shell(mode)
        assert rc == 2 and "FORMAT" in out, out
    assert fp.kills() == []


def test_v15_darwin_unit_is_rejected_by_the_linux_reader_and_the_shell(fp):
    base_tree(fp)
    fp.set_id("mhw-candidate-id v1 darwin-usec 100 100 100 1000\n")
    fp.set_pids([], header="# mhw-candidate-pids v1 darwin-usec")
    assert fp.python("list") == candidate.EXIT_FORMAT
    assert fp.shell("list")[0] == 2


@pytest.mark.parametrize("body", ["", "# wrong header\n", "# mhw-candidate-pids v1 linux-ticks\n101 1001 100 100\n",
                                  "# mhw-candidate-pids v1 linux-ticks\n101 1001 100 100 x\n",
                                  "# mhw-candidate-pids v1 darwin-usec\n"])
def test_v15_bad_candidate_pids(fp, body):
    base_tree(fp)
    fp.set_master(100, 1000)
    (fp.run / "candidate.pids").write_text(body)
    assert fp.python("stop") == candidate.EXIT_FORMAT
    assert fp.shell("stop")[0] == 2
    assert fp.kills() == []


def test_v15_missing_pids_file(fp):
    base_tree(fp)
    fp.set_master(100, 1000)
    assert fp.python("list") == candidate.EXIT_FORMAT
    assert fp.shell("list")[0] == 2


# ------------------------------------------- V16/V17 resolution, Python == shell


def check_parity(fp, expect: dict, list_rc: int):
    sets, _ = py_sets(fp)
    assert sets == expect
    rc, out = fp.shell("list")
    assert shell_sets(out) == expect, out
    assert fp.python("list") == list_rc
    assert rc == list_rc, out


def test_v17_normal_tree(fp):
    base_tree(fp)
    fp.set_master(100, 1000)
    fp.set_pids([(101, 1001, 100, 100, 100)])
    check_parity(fp, {"candidate": {100, 101, 102}, "foreign": set(), "UNDETERMINED": set()}, 1)


def test_v16_v17_stop_signals_master_first_and_ends_clean(fp):
    base_tree(fp)
    fp.set_master(100, 1000)
    fp.set_pids([(101, 1001, 100, 100, 100)])
    assert fp.python("stop") == candidate.EXIT_OK
    py = fp.kills()
    fp.killlog.unlink()
    base_tree(fp)
    rc, out = fp.shell("stop")
    assert rc == 0 and "no candidate process" in out, out
    assert py == fp.kills() == [("TERM", "100"), ("TERM", "101"), ("TERM", "102")]
    assert (fp.root / "500").exists()


def test_v10_stale_state_file_points_at_other_processes(fp):
    base_tree(fp, master_start=777)  # pid 100 is now another process (start differs)
    fp.set_master(100, 1000)
    fp.set_pids([(101, 9999, 100, 100, 100)])  # recorded worker pid reused as well
    check_parity(fp, {"candidate": set(), "foreign": {100, 101, 102}, "UNDETERMINED": set()}, candidate.EXIT_MISMATCH)
    assert fp.python("stop") == candidate.EXIT_MISMATCH
    rc, out = fp.shell("stop")
    assert rc == candidate.EXIT_MISMATCH and "no candidate process" in out
    assert fp.kills() == []
    assert "no candidate process; port" not in out  # never reported as plainly clean


def test_v9b_master_gone_unrecorded_worker_is_undetermined(fp):
    fp.add(101, 1, 100, 100, 1001)  # worker orphaned, never recorded
    fp.set_master(100, 1000)
    fp.set_pids([])
    check_parity(fp, {"candidate": set(), "foreign": set(), "UNDETERMINED": {101}}, candidate.EXIT_UNDETERMINED)
    assert fp.python("stop") == candidate.EXIT_UNDETERMINED
    rc, out = fp.shell("stop")
    assert rc == 6 and "no candidate process" not in out
    assert fp.kills() == []


def test_v9a_master_gone_recorded_worker_is_cleaned(fp):
    fp.add(101, 1, 100, 100, 1001)
    fp.add(103, 101, 100, 100, 1003)  # unrecorded child of the worker: found through the session
    fp.set_master(100, 1000)
    fp.set_pids([(101, 1001, 100, 100, 100)])
    check_parity(fp, {"candidate": {101, 103}, "foreign": set(), "UNDETERMINED": set()}, 1)
    assert fp.python("stop") == candidate.EXIT_OK
    assert sorted(fp.kills()) == [("TERM", "101"), ("TERM", "103")]


def test_v5_fixture_unrecorded_child_of_dead_helper_found_through_session(fp):
    base_tree(fp)
    fp.add(201, 1, 200, 200, 2001)  # B: unrecorded, helper 200 is gone
    fp.add(202, 1, 200, 200, 2002)  # A: recorded anchor in the helper's session
    fp.set_master(100, 1000)
    fp.set_pids([(101, 1001, 100, 100, 100), (200, 2000, 101, 200, 200), (202, 2002, 200, 200, 200)])
    sets, r = py_sets(fp)
    assert sets["candidate"] == {100, 101, 102, 201, 202} and r.how[201] == "session"
    check_parity(fp, sets, 1)
    rc, out = fp.shell("list")
    assert "candidate 201 start=2001 via=session" in out


def test_v5b_fixture_dead_helper_without_anchor_is_undetermined(fp):
    base_tree(fp)
    fp.add(201, 1, 200, 200, 2001)
    fp.set_master(100, 1000)
    fp.set_pids([(101, 1001, 100, 100, 100), (200, 2000, 101, 200, 200)])
    check_parity(fp, {"candidate": {100, 101, 102}, "foreign": set(), "UNDETERMINED": {201}},
                 candidate.EXIT_UNDETERMINED)
    assert fp.python("stop") == candidate.EXIT_UNDETERMINED
    assert fp.shell("stop")[0] == 6
    assert fp.kills() == []  # undetermined blocks every signal, also to provable candidates


def test_v11_fixture_session_number_reused_while_detached_helper_lives(fp):
    fp.add(100, 1, 100, 100, 5000)   # someone else now owns pid/session 100
    fp.add(103, 100, 100, 100, 5001)  # ... and its child
    fp.add(200, 1, 200, 200, 2000)   # our recorded detached helper, still alive
    fp.set_master(100, 1000)
    fp.set_pids([(101, 1001, 100, 100, 100), (200, 2000, 101, 200, 200)])
    check_parity(fp, {"candidate": {200}, "foreign": {100, 103}, "UNDETERMINED": set()}, 1)
    assert fp.python("stop") == candidate.EXIT_MISMATCH
    assert fp.kills() == [("TERM", "200")]
    fp.killlog.unlink()
    fp.add(200, 1, 200, 200, 2000)
    rc, out = fp.shell("stop")
    assert rc == candidate.EXIT_MISMATCH, out
    assert fp.kills() == [("TERM", "200")]
    assert (fp.root / "100").exists() and (fp.root / "103").exists()


def test_v16_zombie_master_counts_as_gone(fp):
    fp.add(100, 1, 100, 100, 1000, state="Z")
    fp.set_master(100, 1000)
    fp.set_pids([])
    check_parity(fp, {"candidate": set(), "foreign": set(), "UNDETERMINED": set()}, 0)


def test_v16_residual_after_kill_is_never_reported_clean(fp):
    base_tree(fp)
    fp.set_master(100, 1000)
    fp.set_pids([])
    assert fp.python("stop", removes=False) == candidate.EXIT_RESIDUAL
    rc, out = fp.shell("stop", removes=False)
    assert rc == 5 and "CANDIDATE STILL ALIVE" in out and "no candidate process" not in out
    assert ("KILL", "100") in fp.kills()


def test_v16_shell_port_still_bound(fp):
    fp.set_master(100, 1000)
    fp.set_pids([])
    rc, out = fp.shell("list", port_out="LISTEN 0 128 127.0.0.1:18031")
    assert rc == 5 and "STILL BOUND" in out


def test_v16_no_candidate_id_lists_cwd_matches_without_signals(fp):
    fp.add(300, 1, 300, 300, 3000)
    os.symlink(str(fp.D / "app" / "API"), fp.root / "300" / "cwd")
    rc, out = fp.shell("stop")
    assert rc == 6 and "300" in out and fp.kills() == []
    assert fp.python("stop") == candidate.EXIT_UNDETERMINED
    fp.remove(300)
    assert fp.shell("list")[0] == 0 and fp.python("list") == 0


# ------------------------------- round 14 #1: an incomplete scan is never clean


@pytest.fixture()
def unreadable(fp):
    """Restore permissions so tmp_path cleanup works whatever the test did."""
    changed = []
    yield changed
    for f in changed:
        if f.exists():
            f.chmod(0o644)


@pytest.mark.parametrize("how", ["chmod000", "garbage", "short"])
@pytest.mark.parametrize("mode", ["list", "stop"])
def test_r14_unreadable_entry_of_ours_blocks_python_and_shell(fp, unreadable, how, mode):
    base_tree(fp)
    fp.set_master(100, 1000)
    fp.set_pids([(101, 1001, 100, 100, 100)])
    fp.add(103, 101, 100, 100, 1003)  # a leftover helper the scan cannot read
    stat = fp.root / "103" / "stat"
    if how == "chmod000":
        stat.chmod(0o000)
        unreadable.append(stat)
    elif how == "garbage":
        stat.write_text("garbage\n")
    else:
        stat.write_text("103 (x) S 101 100 100\n")
    table = fp.table()
    r = pi.scan_resolve(table, *pi.read_state(str(fp.run), pi.UNIT_LINUX))
    assert r.state == "undetermined" and [u["pid"] for u in r.unreadable] == [103]
    assert fp.python(mode) == candidate.EXIT_UNDETERMINED
    rc, out = fp.shell(mode)
    assert rc == 6 and "UNREADABLE 103" in out and "no candidate process" not in out, out
    assert fp.kills() == []  # nothing signalled, also not the provable candidates


def test_r14_unreadable_entry_without_candidate_id(fp, unreadable):
    fp.add(103, 1, 103, 103, 1003)
    stat = fp.root / "103" / "stat"
    stat.chmod(0o000)
    unreadable.append(stat)
    assert fp.python("list") == candidate.EXIT_UNDETERMINED
    rc, out = fp.shell("list")
    assert rc == 6 and "unreadable" in out


def test_r15_entry_present_without_stat_is_undetermined(fp):
    """V21 (r15 policy): an entry that is still there but has no readable stat is
    not "gone" (only a confirmed ENOENT on the entry is); Python and shell agree."""
    base_tree(fp)
    (fp.root / "777").mkdir()
    fp.set_master(100, 1000)
    fp.set_pids([])
    _, r = py_sets(fp)
    assert [u["pid"] for u in r.unreadable] == [777] and r.state == "undetermined"
    assert fp.python("list") == candidate.EXIT_UNDETERMINED
    rc, out = fp.shell("list")
    assert rc == 6 and "UNREADABLE 777" in out, out


def test_r15_confirmed_enoent_on_the_entry_is_gone(fp, unreadable, monkeypatch):
    base_tree(fp)
    fp.add(777, 1, 777, 777, 7)
    (fp.root / "777" / "stat").chmod(0o000)
    unreadable.append(fp.root / "777" / "stat")
    target = str(fp.root / "777")
    real = os.lstat

    def lstat(path, *a, **kw):
        if str(path) == target:
            raise FileNotFoundError(2, "No such file or directory", target)
        return real(path, *a, **kw)

    monkeypatch.setattr(pi.os, "lstat", lstat)
    table = fp.table()
    assert table.entry(777) == (pi.GONE, None)
    table.scan()
    assert table.bad == [] and table.other_unreadable == []


@pytest.mark.parametrize("err", [PermissionError(13, "denied"), OSError(5, "I/O error")])
def test_r15_unknown_owner_is_fail_closed_python(fp, unreadable, monkeypatch, err):
    base_tree(fp)
    fp.set_master(100, 1000)
    fp.set_pids([])
    fp.add(777, 1, 777, 777, 7)
    (fp.root / "777" / "stat").chmod(0o000)
    unreadable.append(fp.root / "777" / "stat")
    target = str(fp.root / "777")
    real = os.lstat

    def lstat(path, *a, **kw):
        if str(path) == target:
            raise err
        return real(path, *a, **kw)

    monkeypatch.setattr(pi.os, "lstat", lstat)
    table = fp.table(euid=4242)  # even with "another uid" configured, an unknown owner blocks
    table.scan()
    assert [b["pid"] for b in table.bad] == [777] and table.bad[0]["owner"] is None
    assert table.other_unreadable == []


def test_r15_unknown_owner_is_fail_closed_shell(fp, unreadable):
    base_tree(fp)
    fp.set_master(100, 1000)
    fp.set_pids([])
    fp.add(777, 1, 777, 777, 7)
    (fp.root / "777" / "stat").chmod(0o000)
    unreadable.append(fp.root / "777" / "stat")
    stub = fp.killlog.parent / "statstub.sh"
    stub.write_text(f"""#!/bin/bash
for a in "$@"; do last=$a; done
if [ "$last" = "{fp.root}/777" ]; then echo "stat: $last: Permission denied" >&2; exit 1; fi
exec stat "$@"
""")
    stub.chmod(0o755)
    rc, out = fp.shell("stop", statcmd=stub, euid=4242)
    assert rc == 6 and "UNREADABLE 777" in out and fp.kills() == [], out


def test_r15_other_uid_unreadable_is_listed_not_blocking(fp, unreadable):
    base_tree(fp)
    fp.set_master(100, 1000)
    fp.set_pids([(101, 1001, 100, 100, 100)])
    fp.add(888, 1, 888, 888, 8)
    (fp.root / "888" / "stat").chmod(0o000)
    unreadable.append(fp.root / "888" / "stat")
    me = os.geteuid()
    table = fp.table(euid=me + 4242)  # the fixture files are ours, so "another uid" from this view
    r = pi.scan_resolve(table, *pi.read_state(str(fp.run), pi.UNIT_LINUX))
    assert r.state == "alive" and r.unreadable == []
    assert [(o["pid"], o["owner"]) for o in r.other_unreadable] == [(888, me)]
    assert r.as_dict()["other_unreadable"][0]["pid"] == 888  # in the recorded result
    assert fp.python("list", euid=me + 4242) == 1
    rc, out = fp.shell("list", euid=me + 4242)
    assert rc == 1 and f"OTHER-UID-UNREADABLE 888 uid={me}" in out, out
    assert shell_sets(out)["candidate"] == {100, 101, 102}
    # the same entry from our own view blocks
    assert fp.python("list") == candidate.EXIT_UNDETERMINED and fp.shell("list")[0] == 6


def test_r14_signaller_does_not_signal_unreadable_process():
    class T:
        unit = pi.UNIT_LINUX

        def read(self, pid):
            raise PermissionError(13, "denied")

    sent = []
    sg = pi.Signaller(T(), kill=lambda pid, sig: sent.append(pid))
    assert sg.send(pi.Proc(100, 1, 100, 100, 1000), signal.SIGTERM).startswith("unreadable")
    assert sent == []


def test_r14_final_stop_error_still_writes_result(live, monkeypatch):
    real = pi.stop_candidate
    calls = {"n": 0}

    def flaky(*a, **kw):
        calls["n"] += 1
        raise OSError(5, "simulated /proc read failure")

    monkeypatch.setattr(candidate.pi, "stop_candidate", flaky)
    rc = run(live)
    assert rc == candidate.EXIT_UNDETERMINED and calls["n"] >= 1
    doc = result(live["D"])
    assert doc["stop_error"]["message"].endswith("simulated /proc read failure")
    assert "exit_check" not in doc  # never reported clean
    monkeypatch.setattr(candidate.pi, "stop_candidate", real)
    assert candidate.recover(live["D"], "stop", live["port"]) == candidate.EXIT_OK  # recovery still works
    assert no_candidate_left(live["D"], live["port"])


def test_r14_shell_refuses_overrides_without_test_mode(fp):
    base_tree(fp)
    fp.set_master(100, 1000)
    fp.set_pids([])
    env = {"PATH": os.environ["PATH"], "MHW_RECOVER_PROC": str(fp.root), "MHW_RECOVER_KILL": "/bin/echo"}
    p = subprocess.run(["bash", str(SHELL), "stop", str(fp.run), "18031"], env=env, capture_output=True,
                       text=True, timeout=30)
    assert p.returncode == 2 and "REFUSED" in p.stdout


# ------------------------------------------------------- V12 identity re-check


def test_v12_pid_reused_between_resolution_and_kill_python_and_shell(fp):
    """TERM to the master makes pid 101 belong to another process (new session)
    before 101's own kill: both implementations must re-read its start and skip it."""
    reused = stat_line(101, 1, 101, 101, 9101, comm="other")

    def setup():
        base_tree(fp)
        fp.set_master(100, 1000)
        fp.set_pids([(101, 1001, 100, 100, 100)])

    setup()
    table = fp.table()
    plain = fp.py_kill(True)

    def kill(pid, sig):
        plain(pid, sig)
        if pid == 100:
            (fp.root / "101" / "stat").write_text(reused)

    rc = candidate.recover(fp.D, "stop", None, table=table, signaller=pi.Signaller(table, kill=kill),
                           term_wait=0, kill_wait=0, poll=0.01)
    assert ("TERM", "101") not in fp.kills() and ("TERM", "102") in fp.kills()
    assert (fp.root / "101").exists() and rc == candidate.EXIT_OK
    fp.killlog.unlink()
    fp.remove(101)
    setup()
    extra = f'[ "$3" = 100 ] && printf %s "{reused.strip()}" > {fp.root}/101/stat'
    rc, out = fp.shell("stop", extra=extra)
    assert ("TERM", "101") not in fp.kills() and ("TERM", "102") in fp.kills(), out
    assert "skip TERM 101" in out and rc == 0


class _SeqTable:
    unit = pi.UNIT_LINUX

    def __init__(self, reads):
        self.reads = list(reads)

    def read(self, pid):
        return self.reads.pop(0)


def test_v12_signaller_rechecks_start_before_every_signal():
    sent = []
    p = pi.Proc(100, 1, 100, 100, 1000)
    table = _SeqTable([p, pi.Proc(100, 1, 100, 100, 4242), None])
    sg = pi.Signaller(table, kill=lambda pid, sig: sent.append((pid, sig)))
    assert sg.send(p, signal.SIGTERM) == "sent"
    assert sg.send(p, signal.SIGKILL) == "changed"  # pid reused between TERM and KILL
    assert sg.send(p, signal.SIGKILL) == "gone"
    assert sent == [(100, signal.SIGTERM)]


def test_v12_fixture_pid_reused_between_term_and_kill(fp):
    base_tree(fp)
    fp.set_master(100, 1000)
    fp.set_pids([(101, 1001, 100, 100, 100)])
    table = fp.table()

    def kill(pid, sig):
        with open(fp.killlog, "a") as f:
            f.write(f"{signal.Signals(sig).name} {pid}\n")
        if sig == signal.SIGTERM and pid == 100:
            # the whole candidate exits after TERM and pid 100 is handed to a new session
            for q in (100, 101, 102):
                fp.remove(q)
            fp.add(100, 1, 100, 100, 7777)

    res = pi.stop_candidate(table, pi.read_state(str(fp.run), pi.UNIT_LINUX)[0],
                            lambda: pi.read_state(str(fp.run), pi.UNIT_LINUX)[1],
                            term_wait=0, kill_wait=0, signaller=pi.Signaller(table, kill=kill))
    log = fp.killlog.read_text().split("\n")
    assert "SIGKILL 100" not in log and res.state == "clean" and res.mismatch
    assert set(res.final.foreign) == {100}


# ==================================================================== live tests

FAKE = r'''
import os, signal, socket, sys, time
args = sys.argv[1:]
def opt(name, default=None):
    return args[args.index(name) + 1] if name in args else default
flag = lambda n: n in args
sock = socket.socket()
sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
sock.bind(("127.0.0.1", int(opt("--port"))))
sock.listen(16)

def sleeper():
    pid = os.fork()
    if pid == 0:
        signal.signal(signal.SIGTERM, signal.SIG_DFL)
        while True:
            time.sleep(1)
    return pid

def helper():
    if os.fork():
        return
    sock.close()
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    if flag("--helper-setsid"):
        os.setsid()
    go = opt("--go")
    if go:
        if not flag("--no-anchor"):
            sleeper()                      # A: recorded while the helper lives
        while not os.path.exists(go):
            time.sleep(0.05)
        b = sleeper()                      # B: created after recording stopped
        with open(go + ".b", "w") as f:
            f.write(str(b))
        os._exit(0)
    while True:
        time.sleep(1)

def worker():
    if flag("--worker-ignore-term"):
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
    else:
        signal.signal(signal.SIGTERM, lambda *a: os._exit(0))
    if flag("--helper"):
        helper()
    hold = bytearray(int(opt("--alloc-mb", "0")) * 1024 * 1024)
    for i in range(0, len(hold), 4096):
        hold[i] = 1
    # readiness modes: the port is open (bound + listening) from the start, but
    # the app does not answer until --not-ready-for has passed (never with
    # --never-ready): "close" drops the connection without a response, "503"
    # answers 503, "hang" reads and never answers
    ready_at = time.time() + float(opt("--not-ready-for", "0"))
    never, mode = flag("--never-ready"), opt("--not-ready-mode", "close")
    marker_log, marker_at, marker_done = opt("--marker-log"), time.time() + float(opt("--marker-after", "0")), False
    req_log = opt("--req-log")
    sock.settimeout(0.2)
    while True:
        if marker_log and not marker_done and time.time() >= marker_at:
            with open(marker_log, "a") as f:
                f.write("[INFO] Application startup complete.\n")
            marker_done = True
        try:
            c, _ = sock.accept()
        except socket.timeout:
            continue
        with c:
            ready = not never and time.time() >= ready_at
            if not ready and mode == "close":
                continue
            c.settimeout(2)
            try:
                req = c.recv(65536)
            except OSError:
                req = b""
            if req_log:
                with open(req_log, "a") as f:
                    f.write(f"{time.time()} {ready} {req.split(b' ')[1].decode() if req.count(b' ') > 1 else '?'}\n")
            if not ready and mode == "hang":
                time.sleep(3)
                continue
            if not ready:
                c.sendall(b"HTTP/1.1 503 Service Unavailable\r\nContent-Length: 0\r\nConnection: close\r\n\r\n")
                continue
            c.sendall(b"HTTP/1.1 200 OK\r\nContent-Type: text/plain\r\nContent-Length: 2\r\nConnection: close\r\n\r\nok")

wpid = os.fork()
if wpid == 0:
    worker()
    os._exit(0)
if flag("--master-helper"):  # a direct child of the master that is not a gunicorn worker
    if os.fork() == 0:
        sock.close()
        os.execv("/bin/sleep", ["sleep", "600"])
if opt("--master-exit-after"):
    time.sleep(float(opt("--master-exit-after")))
    os._exit(0)
def term(*a):
    if opt("--term-marker"):
        with open(opt("--term-marker"), "w") as f:
            f.write(str(os.getpid()))
    try:
        os.kill(wpid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        os.waitpid(wpid, 0)
    except ChildProcessError:
        pass
    os._exit(0)
signal.signal(signal.SIGTERM, term)
while True:
    time.sleep(1)
'''

DRIVER = r'''
import sys, time
sys.path.insert(0, sys.argv[5])
from pathlib import Path
from tools.bench import candidate, load
D, port, fake, mode = Path(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
argv = [sys.executable, "-B", fake, "--port", str(port), "--helper"] + sys.argv[6:]
def blocking(base, sampler, timeout=5):
    with load.httpx.Client(base_url=base, timeout=5) as c:
        while True:
            c.get("/")
            time.sleep(0.2)
if mode == "pause_before_id":
    for sub in candidate.SUBDIRS:
        (D / sub).mkdir(parents=True, exist_ok=True)
    (D / "run" / "s").mkdir()
    l = candidate.Launcher(argv, D / "app" / "API", candidate.candidate_env(D, sys.executable), D / "run" / "s",
                           step="s", port=port,
                           child_pause=(str(D / "run" / "s" / ".test_child_paused"), str(D / "release")))
    print("starting", flush=True)
    l.start()
    print("started", flush=True)
    time.sleep(60)
else:
    sys.exit(candidate.run_step(D, "s", port, observe_mb=100000, abort_rss_mb=100000, abort_avail_mb=0,
                                python=sys.executable, argv=argv, scenario_fn=blocking, startup_timeout=20,
                                term_wait=5))
'''


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    port = s.getsockname()[1]
    s.close()
    return port


@pytest.fixture()
def live(tmp_path):
    fake = tmp_path / "fakesrv.py"
    fake.write_text(FAKE)
    D = tmp_path / "D"
    started: list[int] = []
    yield {"D": D, "fake": str(fake), "port": free_port(), "tmp": tmp_path, "started": started}
    # safety net: never leave fake processes behind, whatever a test did
    app = os.path.realpath(D / "app" / "API")
    for p in psutil.process_iter():
        try:
            if p.pid in started or os.path.realpath(p.cwd()) == app:
                p.kill()
        except (psutil.Error, OSError):
            pass


def fake_argv(lv, *extra):
    return [sys.executable, "-B", lv["fake"], "--port", str(lv["port"]), *extra]


def quick(base, sampler, timeout=5):
    return {"n": len(load.run_threads(base, [["/"], ["/"]], timeout=5))}


def until(cond, timeout=10.0, step=0.05):
    end = time.monotonic() + timeout
    while time.monotonic() < end:
        v = cond()
        if v:
            return v
        time.sleep(step)
    raise AssertionError("condition not met in time")


def no_candidate_left(D, port):
    app = os.path.realpath(Path(D) / "app" / "API")
    left = []
    for p in psutil.process_iter():
        try:
            if os.path.realpath(p.cwd()) == app and p.status() != psutil.STATUS_ZOMBIE:
                left.append(p.pid)
        except (psutil.Error, OSError):
            pass
    return not left and pi.port_free(port)


def result(D, step="s"):
    return json.loads((Path(D) / "run" / step / "result.json").read_text())


def run(lv, *extra, step="s", **kw):
    params = dict(observe_mb=100000, abort_rss_mb=100000, abort_avail_mb=0, python=sys.executable,
                  argv=fake_argv(lv, *extra), scenario_fn=quick, startup_timeout=20, term_wait=5,
                  install_handlers=False)
    params.update(kw)
    return candidate.run_step(lv["D"], step, lv["port"], **params)


def test_v1_normal_stop_leaves_nothing(live):
    assert run(live, "--helper") == candidate.EXIT_OK
    doc = result(live["D"])
    assert doc["exit_check"] == {"state": "clean", "port_free": True, "mismatch": False}
    assert len(doc["monitor"]["recorded"]) >= 3  # master, worker, helper
    assert no_candidate_left(live["D"], live["port"])
    assert (live["D"] / "run" / "s" / "candidate.id").read_text().startswith("mhw-candidate-id v1 ")


def test_v1_step_directory_is_never_reused(live):
    assert run(live) == candidate.EXIT_OK
    assert run(live) == candidate.EXIT_FORMAT


def test_v2_worker_ignoring_term_is_killed(live):
    assert run(live, "--worker-ignore-term", term_wait=1) == candidate.EXIT_OK
    sigs = [s["signal"] for s in result(live["D"])["stops"][-1]["signals"]]
    assert "SIGKILL" in sigs
    assert no_candidate_left(live["D"], live["port"])


def test_v3_watchdog_worker_rss_abort(live):
    rc = run(live, "--alloc-mb", "300", abort_rss_mb=150, scenario_fn=_hammer)
    assert rc == candidate.EXIT_ABORT
    doc = result(live["D"])
    assert doc["monitor"]["abort"]["reason"] == "worker_rss"
    assert doc["exit_check"]["state"] == "clean" and doc["exit_check"]["port_free"]
    assert no_candidate_left(live["D"], live["port"])


def test_v3_watchdog_available_memory_abort(live):
    rc = run(live, "--helper", abort_avail_mb=10 ** 9, scenario_fn=_hammer)
    assert rc == candidate.EXIT_ABORT
    assert result(live["D"])["monitor"]["abort"]["reason"] == "available_memory"
    assert no_candidate_left(live["D"], live["port"])


def _hammer(base, sampler, timeout=5):
    with load.httpx.Client(base_url=base, timeout=5) as c:
        for _ in range(300):
            c.get("/")
            time.sleep(0.1)
    raise load.LoadError("watchdog did not fire")


def test_v4_detached_helper_recorded_and_cleaned(live):
    assert run(live, "--helper", "--helper-setsid") == candidate.EXIT_OK
    doc = result(live["D"])
    sid = doc["candidate_id"]["sid"]
    assert any(p["sid"] != sid for p in doc["monitor"]["recorded"])  # the detached helper was recorded
    assert no_candidate_left(live["D"], live["port"])


def test_observation_line_only_marks(live):
    rc = run(live, "--alloc-mb", "120", observe_mb=50)
    assert rc == candidate.EXIT_OBSERVE
    assert result(live["D"])["observation_line_crossed"] is True
    assert no_candidate_left(live["D"], live["port"])


def _launch(lv, *extra, step="s"):
    D = lv["D"]
    for sub in candidate.SUBDIRS:
        (D / sub).mkdir(parents=True, exist_ok=True)
    rd = D / "run" / step
    rd.mkdir()
    table = pi.platform_table()
    l = candidate.Launcher(fake_argv(lv, *extra), D / "app" / "API", candidate.candidate_env(D, sys.executable),
                           rd, step=step, port=lv["port"], table=table)
    l.start()
    lv["started"].append(l.proc.pid)
    return l, table, rd


@pytest.mark.parametrize("anchor", [True, False])
def test_v5_dead_helper_leaves_unrecorded_session_child(live, anchor):
    go = live["tmp"] / "go"
    extra = ["--helper", "--helper-setsid", "--go", str(go)] + ([] if anchor else ["--no-anchor"])
    l, table, rd = _launch(live, *extra)
    m = candidate.Monitor(table, l.cid, rd, abort_rss_mb=1e9, abort_avail_mb=0, on_abort=lambda i: None)
    m.start()
    want = 2 if anchor else 1  # helper (+ A) in a session other than the candidate's
    until(lambda: sum(1 for p in m.recorded_list() if p.sid != l.cid.sid) >= want)
    m.halt()
    helper = next(p for p in m.recorded_list() if p.sid != l.cid.sid and p.pid == p.sid)
    go.write_text("go")
    b = int(until(lambda: (go.parent / "go.b").exists() and (go.parent / "go.b").read_text()))
    live["started"].append(b)
    until(lambda: table.read(helper.pid) is None)  # the helper has exited
    cid, recorded = pi.read_state(str(rd), table.unit)
    assert b not in {p.pid for p in recorded}  # B was never recorded
    res = pi.resolve(table.scan(), cid, recorded)
    if anchor:
        assert res.how[b] == "session"  # found only through the anchored session
        sr = pi.stop_candidate(table, cid, lambda: recorded, term_wait=5, port=live["port"], reap=l.reap)
        assert sr.state == "clean"
        assert no_candidate_left(live["D"], live["port"]) and table.read(b) is None
    else:
        assert b in res.undetermined
        sr = pi.stop_candidate(table, cid, lambda: recorded, term_wait=5, port=live["port"], reap=l.reap)
        assert sr.state == "undetermined" and sr.signals == []
        assert table.read(b) is not None and table.read(l.cid.master_pid) is not None  # nothing signalled
        assert candidate.recover(live["D"], "stop", live["port"]) == candidate.EXIT_UNDETERMINED


def _driver(lv, mode, *extra):
    p = subprocess.Popen([sys.executable, "-B", "-c", DRIVER, str(lv["D"]), str(lv["port"]), lv["fake"], mode,
                          str(API), *extra], cwd=str(API), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    lv["started"].append(p.pid)
    return p


def _wait_candidate_ready(lv):
    idf = lv["D"] / "run" / "s" / "candidate.id"
    until(lambda: idf.exists(), timeout=20)
    until(lambda: not pi.port_free(lv["port"]), timeout=20)
    rec = lv["D"] / "run" / "s" / "candidate.pids"
    until(lambda: len(rec.read_text().splitlines()) >= 3, timeout=20)  # header + master + worker


@pytest.mark.parametrize("sig", [signal.SIGTERM, signal.SIGHUP])
def test_v6_runner_interrupted_cleans_up(live, sig):
    p = _driver(live, "run")
    _wait_candidate_ready(live)
    time.sleep(0.5)
    p.send_signal(sig)
    out, _ = p.communicate(timeout=60)
    assert p.returncode == candidate.EXIT_FAILED, out
    doc = result(live["D"])
    assert doc["interrupted"] == signal.Signals(sig).name
    assert doc["exit_check"]["state"] == "clean"
    assert no_candidate_left(live["D"], live["port"])


def test_v7_runner_killed_before_ready_then_recover(live):
    p = _driver(live, "run")
    _wait_candidate_ready(live)
    cid = pi.parse_id((live["D"] / "run" / "s" / "candidate.id").read_text())
    p.kill()
    p.wait(10)
    if LINUX:
        # PR_SET_PDEATHSIG: the master gets SIGTERM when the runner dies; recover cleans the rest
        until(lambda: pi.platform_table().read(cid.master_pid) is None, timeout=30)
        assert candidate.recover(live["D"], "list", live["port"]) in (candidate.EXIT_OK, candidate.EXIT_FAILED)
    else:
        assert not no_candidate_left(live["D"], live["port"])  # no PDEATHSIG on macOS: the candidate outlives it
        assert candidate.recover(live["D"], "list", live["port"]) == candidate.EXIT_FAILED
    assert candidate.recover(live["D"], "stop", live["port"]) == candidate.EXIT_OK
    assert no_candidate_left(live["D"], live["port"])


def test_v8_runner_killed_before_identity_written_never_execs(live):
    p = _driver(live, "pause_before_id")
    assert p.stdout.readline().strip() == "starting"
    idf = live["D"] / "run" / "s" / "candidate.id"
    marker = live["D"] / "run" / "s" / ".test_child_paused"
    release = live["D"] / "release"
    # explicit handshake (Codex rounds 14-15): the child announces the pause and then
    # waits for the release file, which is created only after the runner is dead
    child = int(until(lambda: marker.exists() and marker.read_text().strip(), timeout=20))
    live["started"].append(child)
    assert not idf.exists()  # the child is parked before candidate.id
    assert pi.platform_table().read(child) is not None
    p.kill()
    p.wait(10)
    assert pi.platform_table().read(p.pid) is None  # the runner is gone ...
    assert not idf.exists()  # ... and the child still has not passed the pause
    release.write_text("go")
    until(lambda: idf.exists(), timeout=10)  # written by the child itself after the release
    assert idf.stat().st_mtime_ns >= release.stat().st_mtime_ns
    cid = pi.parse_id(idf.read_text())
    assert cid.master_pid == child
    table = pi.platform_table()
    until(lambda: table.read(cid.master_pid) is None, timeout=10)  # it noticed the runner was gone and exited
    assert pi.port_free(live["port"])  # gunicorn/fake server never started
    assert candidate.recover(live["D"], "list", live["port"]) == candidate.EXIT_OK


@pytest.mark.parametrize("recorded", [True, False])
def test_v9_master_exited_worker_alive(live, recorded):
    l, table, rd = _launch(live, "--master-exit-after", "1.5")
    m = None
    if recorded:
        m = candidate.Monitor(table, l.cid, rd, abort_rss_mb=1e9, abort_avail_mb=0, on_abort=lambda i: None)
        m.start()
        until(lambda: any(p.ppid == l.cid.master_pid for p in m.recorded_list()))
    until(lambda: l.proc.poll() is not None)  # master exited (and reaped)
    if m:
        m.halt()
    worker = until(lambda: [p for p in table.scan().values() if p.sid == l.cid.sid])
    live["started"].extend(p.pid for p in worker)
    if recorded:
        assert candidate.recover(live["D"], "stop", live["port"]) == candidate.EXIT_OK
        assert no_candidate_left(live["D"], live["port"])
    else:
        assert candidate.recover(live["D"], "stop", live["port"]) == candidate.EXIT_UNDETERMINED
        assert all(table.read(p.pid) is not None for p in worker)  # not signalled


def _foreign_sleeper(lv):
    p = subprocess.Popen(["sleep", "60"], start_new_session=True)
    lv["started"].append(p.pid)
    return p


def test_v10_stale_files_point_at_live_unrelated_processes(live):
    table = pi.platform_table()
    a, b = _foreign_sleeper(live), _foreign_sleeper(live)
    until(lambda: table.read(a.pid) and table.read(b.pid))
    rd = live["D"] / "run" / "s"
    rd.mkdir(parents=True)
    pa, pb = table.read(a.pid), table.read(b.pid)
    (rd / "candidate.id").write_text(pi.CandidateId(table.unit, a.pid, a.pid, a.pid, pa.start + 1).line())
    (rd / "candidate.pids").write_text(pi.pids_header(table.unit) + pi.pids_line(pi.Proc(b.pid, 1, b.pid, b.pid, pb.start - 1)))
    assert candidate.recover(live["D"], "stop", None) == candidate.EXIT_MISMATCH
    assert a.poll() is None and b.poll() is None  # untouched


def test_v11_number_reused_while_detached_helper_lives(live):
    table = pi.platform_table()
    foreign = _foreign_sleeper(live)
    helper = _foreign_sleeper(live)  # plays our recorded detached helper (own session)
    until(lambda: table.read(foreign.pid) and table.read(helper.pid))
    pf, ph = table.read(foreign.pid), table.read(helper.pid)
    rd = live["D"] / "run" / "s"
    rd.mkdir(parents=True)
    (rd / "candidate.id").write_text(pi.CandidateId(table.unit, pf.pid, pf.pid, pf.pid, pf.start + 1).line())
    (rd / "candidate.pids").write_text(pi.pids_header(table.unit) + pi.pids_line(ph))
    assert candidate.recover(live["D"], "stop", None) == candidate.EXIT_MISMATCH
    helper.wait(10)
    assert foreign.poll() is None  # the process now holding the old number was not signalled


# ============================================================ configuration


# ------------------------------------------------- readiness (Codex r20, P1-2 V-c)


def _probe_and_scenario_times(log):
    rows = [l.split() for l in Path(log).read_text().splitlines()]
    probes = [float(t) for t, ready, path in rows if path == candidate.DEFAULT_PROBE_PATH]
    scen = [float(t) for t, ready, path in rows if path == "/"]
    return rows, probes, scen


@pytest.mark.parametrize("mode", ["close", "503", "hang"])
def test_open_port_without_app_answer_is_not_ready(live, mode):
    log = live["tmp"] / "req.log"
    t0 = time.time()
    assert run(live, "--not-ready-for", "2.5", "--not-ready-mode", mode, "--req-log", str(log)) == candidate.EXIT_OK
    doc = result(live["D"])
    rd = doc["ready"]
    assert rd["ready_after_s"] >= 2.0 and rd["probe_last"]["status"] == 200
    assert rd["probe_attempts"] >= (1 if mode == "hang" else 2)
    rows, probes, scen = _probe_and_scenario_times(log)
    # the scenario (the measurement) starts only after the successful probe, and
    # every probe request happened before it: probes are never measured
    assert scen and probes and max(probes) < min(scen)
    assert all(r[1] == "True" for r in rows if r[2] == "/")
    if mode != "close":
        assert any(r[1] == "False" for r in rows)  # the port answered TCP before the app was ready
    assert min(scen) - t0 >= 2.0
    assert doc["exit_check"]["state"] == "clean" and no_candidate_left(live["D"], live["port"])


def test_never_ready_times_out_and_is_cleaned_up(live):
    rc = run(live, "--never-ready", "--not-ready-mode", "503", "--helper", startup_timeout=3)
    assert rc == candidate.EXIT_FAILED
    doc = result(live["D"])
    assert doc["startup_problem"].startswith("not ready within 3")
    assert "load_exit" not in doc and doc["ready"]["ready_after_s"] is None
    assert doc["ready"]["probe_last"]["status"] == 503
    assert doc["exit_check"] == {"state": "clean", "port_free": True, "mismatch": False}
    assert no_candidate_left(live["D"], live["port"])


def test_ready_waits_for_the_startup_marker_of_every_worker(live):
    log = live["tmp"] / "err.log"
    assert run(live, "--marker-log", str(log), "--marker-after", "2", ready_log=log) == candidate.EXIT_OK
    rd = result(live["D"])["ready"]
    assert rd["markers"] == 1 and rd["ready_after_s"] >= 1.8
    # two workers expected, one marker: never ready, then stopped and checked
    lv2 = dict(live, D=live["tmp"] / "D2")
    log2 = live["tmp"] / "err2.log"
    rc = run(lv2, "--marker-log", str(log2), ready_log=log2, expected_workers=2, workers=2, startup_timeout=3)
    assert rc == candidate.EXIT_FAILED and no_candidate_left(lv2["D"], live["port"])


def test_helper_child_of_master_is_not_counted_as_a_worker(live):
    rc = run(live, "--master-helper", expected_workers=2, startup_timeout=3)
    doc = result(live["D"])
    assert rc == candidate.EXIT_FAILED and len(doc["ready"]["workers"]) == 1
    assert "workers 1/2" in doc["startup_problem"]
    assert doc["exit_check"]["state"] == "clean" and no_candidate_left(live["D"], live["port"])
    # one worker expected: readiness counts only the worker, so the step starts;
    # the load sampler then still rejects the unclassified child (fail closed)
    lv2 = dict(live, D=live["tmp"] / "D2")
    assert run(lv2, "--master-helper", step="t") == candidate.EXIT_FAILED
    doc = result(lv2["D"], "t")
    assert doc.get("startup_problem") is None and doc["ready"]["workers"] and len(doc["ready"]["workers"]) == 1
    assert doc["load_exit"] == 1 and len(doc["monitor"]["recorded"]) >= 3  # master, worker, sleep helper
    assert no_candidate_left(lv2["D"], live["port"])


def test_golden_j01_constant_matches_manifest():
    manifest = json.loads((API / "tests" / "golden" / "manifest.json").read_text())
    assert manifest["cases"]["J01"]["sha256"] == load.GOLDEN_J01_SHA256
    assert manifest["cases"]["J01"]["url"] == load.J01


def test_candidate_env_is_built_from_scratch(tmp_path):
    env = candidate.candidate_env(tmp_path, "/opt/py/bin/python3.11")
    assert not any(k.startswith("DASK") for k in env)
    for k in ("HOME", "XDG_CACHE_HOME", "XDG_CONFIG_HOME", "XDG_DATA_HOME", "XDG_STATE_HOME", "TMPDIR", "MPLCONFIGDIR"):
        assert env[k].startswith(str(tmp_path) + "/")
    assert env["PYTHONDONTWRITEBYTECODE"] == "1" and env["PYTHONNOUSERSITE"] == "1"
    assert env["MHW_HEAVY_SLOTS"] == "1" and env["MHW_HEAVY_BYTES"] == "2097152"
    assert env["PATH"].split(":")[0] == "/opt/py/bin"


def test_gunicorn_argv_single_worker_loopback(tmp_path):
    a = candidate.gunicorn_argv("/py/bin/python3.11", "/py/bin/gunicorn", 18031, tmp_path, tmp_path / "run" / "x", "x")
    assert a[:3] == ["/py/bin/python3.11", "-B", "/py/bin/gunicorn"]
    assert a[a.index("-w") + 1] == "1" and a[a.index("-b") + 1] == "127.0.0.1:18031"
    assert "--reload" not in a and "--keyfile" not in a and "--certfile" not in a
    # P1-2: the worker class is explicit and the deprecated uvicorn.workers is not the default
    assert a[a.index("-k") + 1] == "uvicorn_worker.UvicornWorker" == candidate.DEFAULT_WORKER_CLASS
    assert a[a.index("--error-logfile") + 1] == str(candidate.error_log(tmp_path, "x"))


def test_gunicorn_argv_worker_class_and_workers_are_configurable(tmp_path):
    a = candidate.gunicorn_argv("/p", "/g", 1, tmp_path, tmp_path, "x", worker_class="uvicorn.workers.UvicornWorker",
                                workers=2, app="tools.p1_2.asgi_echo:app")
    assert a[a.index("-k") + 1] == "uvicorn.workers.UvicornWorker" and a[a.index("-w") + 1] == "2"
    assert a[3] == "tools.p1_2.asgi_echo:app"
    with pytest.raises(ValueError):
        candidate.gunicorn_argv("/p", "/g", 1, tmp_path, tmp_path, "x", worker_class="")


def test_cli_passes_worker_class_and_probe(monkeypatch):
    got = {}
    monkeypatch.setattr(candidate, "run_step", lambda *a, **k: got.update(k) or 0)
    base = ["run", "--dir", "/x", "--step", "smoke", "--port", "1", "--python", "/p", "--gunicorn", "/g",
            "--observe-worker-rss-mb", "1", "--abort-worker-rss-mb", "2", "--abort-min-available-mb", "3"]
    assert candidate.main(base) == 0
    assert got["worker_class"] == "uvicorn_worker.UvicornWorker" and got["probe_path"] == candidate.DEFAULT_PROBE_PATH
    assert got["workers"] == 1 and got["expected_workers"] == 1
    assert candidate.main(base + ["--worker-class", "x.Y", "--workers", "2", "--probe-path", "/p"]) == 0
    assert (got["worker_class"], got["workers"], got["expected_workers"], got["probe_path"]) == ("x.Y", 2, 2, "/p")


def test_exit_precedence():
    assert candidate.worst([candidate.EXIT_OBSERVE, candidate.EXIT_RESIDUAL]) == candidate.EXIT_RESIDUAL
    assert candidate.worst([candidate.EXIT_ABORT, candidate.EXIT_UNDETERMINED]) == candidate.EXIT_UNDETERMINED
    assert candidate.worst([candidate.EXIT_FAILED, candidate.EXIT_OBSERVE]) == candidate.EXIT_FAILED
    assert candidate.worst([]) == candidate.EXIT_OK


def test_combo_overlap_rule(monkeypatch):
    plot = {"path": load.HEAVY["plot_month_large"], "start": 10.0, "end": 12.0}
    reqs = [plot] + [{"path": load.SUBTHRESHOLD["half_1m_sst"], "start": 10.1, "end": 10.5}] * 4 + \
           [{"path": load.J01, "start": 12.5, "end": 12.6}]  # J01 after the plot: no overlap
    monkeypatch.setattr(load, "run_threads", lambda base, jobs, timeout: reqs)
    monkeypatch.setattr(load.time, "sleep", lambda s: None)

    class S:
        samples: list = []

        def snapshot(self):
            return {"t": 0, "procs": [{"pid": 1, "role": "worker", "rss_mb": 1.0}]}

    with pytest.raises(load.LoadError, match="overlap"):
        load.scenario_combo_p11("http://x", S())
    reqs[-1] = {"path": load.J01, "start": 11.0, "end": 11.1}
    out = load.scenario_combo_p11("http://x", S())["combo_p11"]
    assert all(o["overlaps_plot"] for o in out["overlap"]) and len(out["paths"]) == 6


# ============================================================ Linux only (L1)


@linux_only
def test_linux_real_proc_scan_is_complete_for_our_uid():
    table = pi.LinuxProcTable()
    snap = table.scan()
    me = table.read(os.getpid())
    assert me is not None and snap[os.getpid()] == me
    assert me.sid == os.getsid(0) and me.pgid == os.getpgid(0)
    assert table.bad == []  # every entry of our uid is readable
    # record (not assert) how many entries of other uids are unreadable (hidepid etc.)
    print(json.dumps({"entries": len(snap), "other_unreadable": len(table.other_unreadable)}))


@linux_only
def test_linux_pidfd_is_used_for_signals(live):
    funcs = pi.pidfd_functions()
    assert funcs is not None  # os wrappers or the ctypes syscalls
    print(json.dumps({"pidfd_source": funcs[2], "os_has_pidfd_open": hasattr(os, "pidfd_open")}))
    assert run(live, "--helper") == candidate.EXIT_OK
    last = result(live["D"])["stops"][-1]
    sent = [s for st in result(live["D"])["stops"] for s in st["signals"] if s["outcome"] == "sent"]
    assert sent and all(s["via"] == "pidfd" for s in sent)
    assert last["pidfd_unavailable"] is False and last["pidfd_source"] == funcs[2]
    assert no_candidate_left(live["D"], live["port"])


@linux_only
@pytest.mark.parametrize("prefer_os", [True, False])
def test_linux_pidfd_syscalls_signal_only_the_opened_process(prefer_os):
    funcs = pi.pidfd_functions(prefer_os=prefer_os)
    if funcs is None or (prefer_os and funcs[2] != "os"):
        pytest.skip("os.pidfd_open not in this interpreter build")
    popen, psend, source = funcs
    assert source == ("os" if prefer_os else "ctypes-syscall")
    p = subprocess.Popen(["sleep", "60"])
    fd = popen(p.pid)
    try:
        psend(fd, signal.SIGTERM)
    finally:
        os.close(fd)
    assert p.wait(10) == -signal.SIGTERM
    with pytest.raises(ProcessLookupError):
        popen(p.pid)  # reaped: the pid is gone


@linux_only
def test_linux_pdeathsig_terms_master_when_runner_is_killed(live):
    marker = live["tmp"] / "master_got_term"
    p = _driver(live, "run", "--term-marker", str(marker))
    _wait_candidate_ready(live)
    cid = pi.parse_id((live["D"] / "run" / "s" / "candidate.id").read_text())
    p.kill()
    p.wait(10)
    until(lambda: marker.exists(), timeout=30)  # SIGTERM reached the master only through PR_SET_PDEATHSIG
    assert int(marker.read_text()) == cid.master_pid
    until(lambda: pi.platform_table().read(cid.master_pid) is None, timeout=30)
    assert candidate.recover(live["D"], "stop", live["port"]) == candidate.EXIT_OK
    assert no_candidate_left(live["D"], live["port"])


@linux_only
@pytest.mark.parametrize("setsid_helper", [False, True])
def test_linux_shell_recover_on_real_proc(live, setsid_helper):
    extra = ["--helper"] + (["--helper-setsid"] if setsid_helper else [])
    l, table, rd = _launch(live, *extra)
    m = candidate.Monitor(table, l.cid, rd, abort_rss_mb=1e9, abort_avail_mb=0, on_abort=lambda i: None)
    m.start()
    until(lambda: len(m.recorded_list()) >= 3)  # master, worker, helper
    m.halt()
    pids = {p.pid for p in m.recorded_list()}
    clean = {"PATH": "/usr/bin:/bin"}  # real use: env -i
    cmd = ["bash", str(SHELL), "list", str(rd), str(live["port"])]
    p = subprocess.run(cmd, env=clean, capture_output=True, text=True, timeout=60)
    assert p.returncode == 1 and pids <= shell_sets(p.stdout)["candidate"], p.stdout
    cmd[2] = "stop"
    p = subprocess.run(cmd, env=clean, capture_output=True, text=True, timeout=120)
    l.reap()
    assert p.returncode == 0 and "no candidate process" in p.stdout, p.stdout
    assert all(table.read(x) is None for x in pids)
    assert no_candidate_left(live["D"], live["port"])
