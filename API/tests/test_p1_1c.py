"""P1-1c targeted checks (spec 17 P1-1c; r2 after Codex round 9): heavy-request lane.

Heavy work waits for a heavy-lane token in the event loop and then runs in a
thread of that lane; light work uses the default threadpool.  Mutual exclusion
and non-starvation are proven with explicit synchronisation points (guard
hooks and events), never inferred from timing.  Every wait has a timeout and
every blocked thread is released and joined in ``finally``.  Responses are
compared with the serial golden manifest.
"""
import ast
import asyncio
import hashlib
import struct
import inspect
import json
import threading
import time
from pathlib import Path

import pytest

from src import guard as G

API = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((API / "tests/golden/manifest.json").read_text())
CASES = MANIFEST["cases"]
WAIT = 60  # generous upper bound for one local request; only reached on failure


def url(cid):
    return CASES[cid]["url"]


PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def assert_golden(cid, r):
    """Data cases (JSON / CSV / errors): body byte-identical to golden.  PNG cases:
    status, headers, PNG signature and image size only -- PNG bytes are not compared
    across matplotlib versions (owner decision 2026-11-25, doc 29); the data behind
    the plots is covered by tests/test_p1_1d.py."""
    c = CASES[cid]
    assert r.status_code == c["status"], (cid, r.status_code)
    for k, v in c["headers"].items():
        assert r.headers.get(k) == v, (cid, k)
    if c["headers"].get("content-type") == "image/png":
        assert r.content[:8] == PNG_SIGNATURE, cid
        assert struct.unpack(">II", r.content[16:24]) == (c["summary"]["width"], c["summary"]["height"]), cid
    else:
        assert hashlib.sha256(r.content).hexdigest() == c["sha256"], cid


def wait_until(cond, timeout=WAIT, what="condition"):
    deadline = time.monotonic() + timeout
    while not cond():
        assert time.monotonic() < deadline, f"timed out waiting for {what}"
        time.sleep(0.01)


@pytest.fixture()
def guard_cfg():
    """Set slots / threshold for one test; restore defaults and clear the hook after."""
    def set_(slots=G.DEFAULT_SLOTS, heavy_bytes=G.DEFAULT_HEAVY_BYTES):
        G.configure(slots=slots, heavy_bytes=heavy_bytes)
        G.reset_stats()
    yield set_
    G.hook = None
    G.configure(slots=G.DEFAULT_SLOTS, heavy_bytes=G.DEFAULT_HEAVY_BYTES)
    G.reset_stats()


class Runner:
    """Run client requests in threads; collect responses or exceptions."""

    def __init__(self, client):
        self.client, self.threads, self.results = client, [], {}

    def start(self, key, path):
        def run():
            try:
                self.results[key] = self.client.get(path)
            except BaseException as e:  # surfaced by join()
                self.results[key] = e
        t = threading.Thread(target=run, daemon=True)
        t.key = key
        self.threads.append(t)
        t.start()
        return t

    def join(self, timeout=WAIT):
        deadline = time.monotonic() + timeout
        for t in self.threads:
            t.join(max(0.0, deadline - time.monotonic()))
        stuck = [t.key for t in self.threads if t.is_alive()]
        assert not stuck, f"requests did not finish in {timeout}s: {stuck}"
        for k, v in self.results.items():
            if isinstance(v, BaseException):
                raise AssertionError(f"request {k} raised") from v
        return self.results


def borrowed(client, limiter_fn):
    """Tokens currently borrowed from a limiter, read inside the app's event loop."""
    return client.portal.call(lambda: limiter_fn().borrowed_tokens)


# ------------------------------------------------------------------ structure


def test_data_and_plot_work_is_synchronous_and_endpoints_only_dispatch():
    from src import mhw_plot, mhw_utils

    for fn in (mhw_utils.process_mhw_data, mhw_utils.select_mhw_data, mhw_utils.finish_mhw_data,
               mhw_plot.month_climatology, mhw_plot.region_climatology):
        assert not inspect.iscoroutinefunction(fn), fn.__name__
    tree = ast.parse((API / "mhw_app.py").read_text())
    heavy_calls = {"finish_mhw_data", "process_mhw_data", "month_climatology", "region_climatology",
                   "select_mhw_data"}
    eps = [n for n in tree.body if isinstance(n, ast.AsyncFunctionDef)
           and n.name in ("read_mhw", "read_mhw_csv", "climatology")]
    assert len(eps) == 3
    for ep in eps:
        # awaits only hand work to a thread
        awaited = {n.value.func.id for n in ast.walk(ep) if isinstance(n, ast.Await)}
        assert awaited <= {"run_in_threadpool", "run_guarded"}, (ep.name, awaited)
        # data / plot functions are called only inside the nested work() (a thread),
        # or passed to run_in_threadpool, never directly in the event loop
        nested = {id(x) for f in ast.walk(ep) if isinstance(f, ast.FunctionDef) for x in ast.walk(f)}
        direct = [n.func.id for n in ast.walk(ep) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                  and n.func.id in heavy_calls and id(n) not in nested]
        assert direct == [], (ep.name, direct)
    for rel in ("src/mhw_utils.py", "src/mhw_plot.py"):
        t = ast.parse((API / rel).read_text())
        assert not [n for n in ast.walk(t) if isinstance(n, (ast.Await, ast.AsyncFunctionDef))], rel


def test_select_reads_no_data_and_matches_process(client):
    """select builds a lazy selection only; finish(select(..)) == process(..)."""
    from datetime import date

    from polars.testing import assert_frame_equal

    from src.mhw_utils import finish_mhw_data, process_mhw_data, select_mhw_data
    from tests.test_p1_1a import COORDS
    from tests.zarr_spy import chunk_reads, data_var_keys

    kw = dict(lon0=-5, lat0=-5, lon1=5, lat1=5, start=date(2020, 1, 1), end=date(2021, 12, 1), append="sst")
    with chunk_reads() as keys:
        sel = select_mhw_data(**kw)
    assert data_var_keys(keys, COORDS) == []
    assert len(sel.parts) == 2 and sel.nbytes == sum(p.nbytes for p in sel.parts) > 0
    a, b = finish_mhw_data(sel), process_mhw_data(**kw)
    assert a[0] == b[0]
    assert_frame_equal(a[1], b[1], check_exact=True)


def test_run_guarded_releases_on_error_and_counts_once(guard_cfg):
    guard_cfg(slots=1)

    async def scenario():
        with pytest.raises(ZeroDivisionError):
            await G.run_guarded(lambda: 1 / 0, heavy=True)
        assert G._heavy_limiter().borrowed_tokens == 0
        assert await G.run_guarded(lambda: 42, heavy=True) == 42
        assert await G.run_guarded(lambda: 7, heavy=False) == 7

    asyncio.run(scenario())
    assert G.stats == {"waiting": 0, "active": 0, "max_active": 1, "acquired": 2, "released": 2}


# ------------------------------------------------------ slot before reading


@pytest.mark.parametrize("cid", ["J07", "J09", "C03", "J08"])
def test_slot_taken_before_any_data_read(client, guard_cfg, cid):
    """J07/C03/J08 go through xr.concat (reads both parts at once); J09 is one bbox."""
    from tests.test_p1_1a import COORDS
    from tests.zarr_spy import chunk_reads, data_var_keys

    guard_cfg(slots=1, heavy_bytes=1)
    at_acquire = []
    with chunk_reads() as keys:
        G.hook = lambda ev, t: at_acquire.append(data_var_keys(list(keys), COORDS)) if ev == "acquired" else None
        r = client.get(url(cid))
    assert_golden(cid, r)
    assert len(at_acquire) == 1, "exactly one acquire per request"
    assert at_acquire[0] == [], f"data read before the slot: {at_acquire[0][:5]}"
    assert data_var_keys(keys.fetched, COORDS), "the spy saw no data fetched at all"  # data was read afterwards


def test_no_data_read_while_plot_lock_held(client, guard_cfg):
    """Requests run one at a time here, so "the plot lock is held" is exactly
    G.PLOT_LOCK.locked() -- valid whichever thread performs the read (zarr 3
    reads on its own IO thread, so a thread-local check would be vacuous)."""
    from tests.zarr_spy import chunk_reads, data_keys

    guard_cfg()
    # Serial execution is what makes PLOT_LOCK.locked() mean "this request holds
    # it": nothing else holds the lock or a heavy slot before we start, the two
    # requests are sent one after the other from this thread, and the guard's
    # own counters show that no two heavy requests overlapped.
    assert not G.PLOT_LOCK.locked() and G.stats["active"] == G.stats["waiting"] == 0
    under_lock = []
    with chunk_reads(on_read=lambda k: under_lock.append(k) if G.PLOT_LOCK.locked() else None) as reads:
        for cid in ("P03", "P01"):
            assert_golden(cid, client.get(url(cid)))
    assert G.stats == {"waiting": 0, "active": 0, "max_active": 1, "acquired": 2, "released": 2}, G.stats
    assert not G.PLOT_LOCK.locked()
    assert data_keys(reads.fetched), "the spy saw no data fetched"
    assert data_keys(under_lock) == [], under_lock[:5]


def test_plot_takes_slot_first_and_lock_after_data(client, guard_cfg):
    guard_cfg()
    events = []
    G.hook = lambda ev, t: events.append(ev)
    assert_golden("P03", client.get(url("P03")))  # two bboxes
    assert events == ["acquiring", "acquired", "plot_locked", "releasing"]


# ---------------------------------------------------- spec acceptance 1 - 4


def _two_heavy(client, slots, first_cid="J07", second_cid="P03"):
    """Start first; while it holds a token, start second and wait until it asks for one."""
    first_in, second_asking, release_first = threading.Event(), threading.Event(), threading.Event()
    order, lock, first = [], threading.Lock(), []

    def hook(ev, ticket):
        with lock:
            order.append((id(ticket), ev))
        if ev == "acquired" and not first:
            first.append(ticket)
            first_in.set()
            assert release_first.wait(WAIT), "first request was never released"
        elif ev == "acquiring" and first and ticket is not first[0]:
            second_asking.set()

    G.hook = hook
    run, err = Runner(client), None
    try:
        run.start(first_cid, url(first_cid))
        assert first_in.wait(WAIT), "first request never acquired"
        run.start(second_cid, url(second_cid))
        assert second_asking.wait(WAIT), "second request never asked for a slot"
        if slots == 1:
            time.sleep(0.3)  # give a broken guard the chance to let it in
            assert G.stats["active"] == 1 and G.stats["waiting"] == 1, G.stats
        else:
            wait_until(lambda: G.stats["active"] == 2, what="two active heavy requests")
    except BaseException as e:
        err = e
    finally:
        release_first.set()
    res = run.join()  # always join, so no request leaks into the next test
    if err is not None:
        raise err
    return res, order, first[0]


def test_slots_1_two_heavy_requests_never_overlap(client, guard_cfg):
    guard_cfg(slots=1, heavy_bytes=1)
    res, order, first = _two_heavy(client, slots=1)
    assert G.stats["max_active"] == 1 and G.stats["active"] == 0 and G.stats["waiting"] == 0
    second = next(t for t, ev in order if ev == "acquired" and t != id(first))
    # the second started only after the first finished (serialisation / savefig done)
    assert order.index((id(first), "releasing")) < order.index((second, "acquired"))
    assert_golden("J07", res["J07"])
    assert_golden("P03", res["P03"])


def test_slots_2_two_heavy_requests_run_together(client, guard_cfg):
    guard_cfg(slots=2, heavy_bytes=1)
    res, _, _ = _two_heavy(client, slots=2)
    assert G.stats["max_active"] == 2
    assert_golden("J07", res["J07"])
    assert_golden("P03", res["P03"])


def test_error_after_slot_releases_it(client, guard_cfg):
    guard_cfg(slots=1, heavy_bytes=1)
    r = client.get(url("E02b"))  # land bbox + sst: 400 raised while holding the slot
    assert_golden("E02b", r)
    assert G.stats["acquired"] == 1 and G.stats["active"] == 0
    assert borrowed(client, G._heavy_limiter) == 0, "slot still held after the 400"
    run = Runner(client)
    run.start("J07", url("J07"))
    assert_golden("J07", run.join()["J07"])
    assert G.stats["acquired"] == 2


def test_light_request_not_blocked_by_held_slot(client, guard_cfg):
    guard_cfg()  # default threshold and slots
    held, release = threading.Event(), threading.Event()

    def hook(ev, t):
        if ev == "acquired":
            held.set()
            assert release.wait(WAIT)

    G.hook = hook
    run, light, err = Runner(client), Runner(client), None
    try:
        run.start("P01", url("P01"))
        assert held.wait(WAIT)
        light.start("J01", url("J01"))  # same app, while the only slot is held
        assert_golden("J01", light.join()["J01"])  # bounded: a blocked J01 fails here
        assert G.stats["acquired"] == 1 and G.stats["active"] == 1  # J01 took no slot
    except BaseException as e:
        err = e
    finally:
        release.set()
    res = run.join()
    light.join()
    if err is not None:
        raise err
    assert_golden("P01", res["P01"])


# ------------------------------------ Codex round 9: no threadpool starvation


def test_queued_heavy_requests_do_not_starve_the_threadpool(client, guard_cfg):
    """Shrink the default threadpool to 3 threads, hold the only heavy slot and
    queue 8 more heavy requests: they wait in the event loop without taking a
    default thread, so a light request still runs (r1 would block here)."""
    import anyio.to_thread

    guard_cfg(slots=1)
    default = anyio.to_thread.current_default_thread_limiter
    old_total = client.portal.call(lambda: default().total_tokens)
    client.portal.call(lambda: setattr(default(), "total_tokens", 3))
    held, release = threading.Event(), threading.Event()

    def hook(ev, t):
        if ev == "acquired" and not held.is_set():
            held.set()
            assert release.wait(WAIT)

    G.hook = hook
    heavy, light, err, res = Runner(client), Runner(client), None, {}
    queued = ["P01", "P06", "P02", "P04"] * 2
    try:
        heavy.start("hold-P01", url("P01"))
        assert held.wait(WAIT)
        for i, cid in enumerate(queued):
            heavy.start(f"q{i}-{cid}", url(cid))
        wait_until(lambda: G.stats["waiting"] == len(queued), what="8 queued heavy requests")
        # the queued ones hold no default thread (the slot holder runs in the heavy lane)
        assert borrowed(client, default) == 0
        light.start("J01", url("J01"))
        assert_golden("J01", light.join()["J01"])
        assert G.stats["active"] == 1 and G.stats["waiting"] == len(queued)
    except BaseException as e:
        err = e
    finally:
        release.set()
        try:
            res = heavy.join(timeout=WAIT * 2)
            light.join()
        finally:
            client.portal.call(lambda: setattr(default(), "total_tokens", old_total))
    if err is not None:
        raise err
    for key, r in res.items():
        assert_golden(key.split("-")[-1], r)
    assert G.stats["max_active"] == 1 and G.stats["acquired"] == 1 + len(queued)


# ------------------------------------------------- concurrent golden


@pytest.mark.parametrize("round_", range(3))
def test_concurrent_golden_mix(client, guard_cfg, round_):
    guard_cfg()
    run = Runner(client)
    for cid in ("J01", "J09", "P01", "P03"):
        run.start(cid, url(cid))
    res = run.join()
    for cid, r in res.items():
        assert_golden(cid, r)
    assert G.stats["active"] == 0 and G.stats["max_active"] == 1


@pytest.mark.parametrize("round_", range(3))
def test_slots_1_concurrent_multi_bbox_plots(client, guard_cfg, round_):
    guard_cfg(slots=1)
    run = Runner(client)
    for cid in ("P03", "P05"):
        run.start(cid, url(cid))
    res = run.join()
    assert_golden("P03", res["P03"])
    assert_golden("P05", res["P05"])
    assert G.stats["max_active"] == 1 and G.stats["acquired"] == 2


def test_concurrent_csv_keeps_cleanup(client, guard_cfg, tmp_path, monkeypatch):
    import os
    import tempfile

    d = tmp_path / "t"
    d.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(d))
    guard_cfg(slots=1, heavy_bytes=1)
    run = Runner(client)
    for i, cid in enumerate(("C01", "C03", "C05")):
        run.start(f"{cid}-{i}", url(cid))
    for key, r in run.join().items():
        assert_golden(key.split("-")[0], r)
    assert os.listdir(d) == []


@pytest.mark.parametrize("path,heavy", [
    # just above / below the 2 MiB default (doc 23): 4 concurrent global 1-month
    # sst requests reached 1,872 MB unguarded, above the 1.5 GB line
    ("/api/mhw?lon0=-180&lon1=180&lat0=-90&lat1=90&start=2026-10-01&end=2026-10-01&append=sst", True),   # 3.97 MiB
    ("/api/mhw?lon0=-180&lon1=-0.25&lat0=-90&lat1=90&start=2026-10-01&end=2026-10-01&append=sst", False),  # 1.98 MiB
    ("/api/mhw?lon0=130&lon1=140&lat0=10&lat1=20&start=2010-01-01&end=2019-12-01&append=sst_anomaly", False),  # J09
    ("/api/mhw?lon0=135&lat0=15&start=1982-01-01&end=2026-10-01", False),  # J01
])
def test_default_threshold_classification(client, guard_cfg, path, heavy):
    guard_cfg()  # defaults
    assert G.HEAVY_BYTES == 2 * 1024 * 1024 and G.SLOTS == 1
    assert client.get(path).status_code == 200
    assert G.stats["acquired"] == (1 if heavy else 0)
