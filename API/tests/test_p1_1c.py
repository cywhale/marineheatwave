"""P1-1c targeted checks (spec 17 P1-1c): threadpool endpoints + heavy-request guard.

Mutual exclusion is proven with explicit synchronisation points (guard hooks
and events), never inferred from timing.  Every wait has a timeout and every
blocked thread is released in ``finally``, so a broken guard fails the test
instead of hanging it.  Responses are compared with the serial golden manifest.
"""
import ast
import hashlib
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


def assert_golden(cid, r):
    c = CASES[cid]
    assert r.status_code == c["status"], (cid, r.status_code)
    assert hashlib.sha256(r.content).hexdigest() == c["sha256"], cid
    for k, v in c["headers"].items():
        assert r.headers.get(k) == v, (cid, k)


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

    def start(self, key, path, headers=None):
        def run():
            try:
                self.results[key] = self.client.get(path, headers=headers)
            except BaseException as e:  # surfaced by join()
                self.results[key] = e
        t = threading.Thread(target=run, name=f"req-{key}", daemon=True)
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


# ------------------------------------------------------------------ structure


def test_data_and_plot_paths_are_synchronous():
    from src import mhw_plot, mhw_utils
    import mhw_app

    for fn in (mhw_utils.process_mhw_data, mhw_plot.month_climatology, mhw_plot.region_climatology,
               mhw_app.read_mhw, mhw_app.read_mhw_csv, mhw_app.climatology):
        assert not inspect.iscoroutinefunction(fn), fn.__name__
    for rel in ("mhw_app.py", "src/mhw_utils.py", "src/mhw_plot.py"):
        tree = ast.parse((API / rel).read_text())
        assert not [n for n in ast.walk(tree) if isinstance(n, ast.Await)
                    and rel != "mhw_app.py"], rel
        calls = {f"{n.func.value.id}.{n.func.attr}" for n in ast.walk(tree)
                 if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                 and isinstance(n.func.value, ast.Name)}
        assert "asyncio.run" not in calls, rel


def test_guard_acquires_once_and_releases_on_error(guard_cfg):
    guard_cfg(slots=1, heavy_bytes=10)
    with G.RequestGuard() as g:
        g.claim(5)
        assert not g.held
        g.claim(11)
        g.claim(1000)
        assert g.held and G.stats["acquired"] == 1
    assert G.stats == {"active": 0, "max_active": 1, "acquired": 1, "released": 1}
    with pytest.raises(ZeroDivisionError):
        with G.RequestGuard(heavy=True):
            1 / 0
    assert G.stats["active"] == 0 and G.stats["released"] == 2
    # the single slot is free again: a non-blocking acquire succeeds
    assert G._slots.acquire(blocking=False)
    G._slots.release()


def test_slot_while_holding_plot_lock_is_refused(guard_cfg):
    guard_cfg(slots=1)
    with G.plot_lock():
        with pytest.raises(RuntimeError, match="lock order"):
            with G.RequestGuard(heavy=True):
                pass
    assert G.stats["active"] == 0 and not G.PLOT_LOCK.locked()


# ------------------------------------------------------ claim before reading


@pytest.mark.parametrize("cid", ["J07", "J09", "C03", "J08"])
def test_claim_happens_before_any_data_read(client, guard_cfg, cid):
    """J07/C03/J08 go through xr.concat (materialises at once); J09 is a single bbox."""
    from tests.test_p1_1a import COORDS, chunk_reads

    guard_cfg(slots=1, heavy_bytes=1)
    at_acquire = []
    with chunk_reads() as keys:
        G.hook = lambda ev, g: at_acquire.append([k for k in keys if k.split("/")[1] not in COORDS]) \
            if ev == "acquired" else None
        r = client.get(url(cid))
    assert_golden(cid, r)
    assert len(at_acquire) == 1, "exactly one acquire per request"
    assert at_acquire[0] == [], f"data read before the claim: {at_acquire[0][:5]}"
    assert any(k.split("/")[1] not in COORDS for k in keys)  # data was read afterwards


def test_no_data_read_while_plot_lock_held(client, guard_cfg):
    from tests.test_p1_1a import chunk_reads

    guard_cfg()
    under_lock = []
    with chunk_reads() as keys:
        import zarr.storage
        spy = zarr.storage.DirectoryStore.__getitem__

        def check(self, key):
            if getattr(G._local, "plot_lock_depth", 0):
                under_lock.append(key)
            return spy(self, key)

        zarr.storage.DirectoryStore.__getitem__ = check
        try:
            for cid in ("P03", "P01"):
                assert_golden(cid, client.get(url(cid)))
        finally:
            zarr.storage.DirectoryStore.__getitem__ = spy
    assert keys and under_lock == []


def test_plot_takes_slot_on_entry_and_lock_after_data(client, guard_cfg):
    guard_cfg()
    events = []
    G.hook = lambda ev, g: events.append(ev)
    assert_golden("P03", client.get(url("P03")))  # two bboxes
    assert events == ["acquiring", "acquired", "plot_locked", "releasing"]


# ---------------------------------------------------- spec acceptance 1 - 4


def _two_heavy(client, slots, first_cid="J07", second_cid="P03"):
    """Start first; while it holds its slot, start second and wait until it asks for one."""
    first_in, second_asking, release_first = threading.Event(), threading.Event(), threading.Event()
    order, lock = [], threading.Lock()
    first_thread = []

    def hook(ev, g):
        me = threading.current_thread()
        with lock:
            order.append((me.ident, ev))  # AnyIO worker threads share one name
        if ev == "acquired" and not first_thread:
            first_thread.append(me)
            first_in.set()
            assert release_first.wait(WAIT), "first request was never released"
        elif ev == "acquiring" and first_thread and me is not first_thread[0]:
            second_asking.set()

    G.hook = hook
    run = Runner(client)
    err = None
    try:
        run.start(first_cid, url(first_cid))
        assert first_in.wait(WAIT), "first request never acquired"
        run.start(second_cid, url(second_cid))
        assert second_asking.wait(WAIT), "second request never asked for a slot"
        if slots == 1:
            time.sleep(0.3)  # give a broken guard the chance to let it in
            assert G.stats["active"] == 1, "second request entered while the only slot was held"
        else:
            deadline = time.monotonic() + WAIT
            while G.stats["active"] < 2 and time.monotonic() < deadline:
                time.sleep(0.01)
            assert G.stats["active"] == 2, "second request could not run alongside the first"
    except BaseException as e:
        err = e
    finally:
        release_first.set()
    res = run.join()  # always join, so no request leaks into the next test
    if err is not None:
        raise err
    return res, order


def test_slots_1_two_heavy_requests_never_overlap(client, guard_cfg):
    guard_cfg(slots=1, heavy_bytes=1)
    res, order = _two_heavy(client, slots=1)
    assert G.stats["max_active"] == 1 and G.stats["active"] == 0
    # the second acquired only after the first started releasing (i.e. after
    # its serialisation / savefig finished)
    names = [n for n, ev in order if ev == "acquired"]
    first, second = names[0], names[1]
    assert order.index((first, "releasing")) < order.index((second, "acquired"))
    assert_golden("J07", res["J07"])
    assert_golden("P03", res["P03"])


def test_slots_2_two_heavy_requests_run_together(client, guard_cfg):
    guard_cfg(slots=2, heavy_bytes=1)
    res, _ = _two_heavy(client, slots=2)
    assert G.stats["max_active"] == 2
    assert_golden("J07", res["J07"])
    assert_golden("P03", res["P03"])


def test_error_after_claim_releases_slot(client, guard_cfg):
    guard_cfg(slots=1, heavy_bytes=1)
    r = client.get(url("E02b"))  # land bbox + sst: 400 after the claim
    assert_golden("E02b", r)
    assert G.stats["acquired"] == 1 and G.stats["active"] == 0
    assert G._slots.acquire(blocking=False), "slot still held after the 400"
    G._slots.release()
    run = Runner(client)
    run.start("J07", url("J07"))
    assert_golden("J07", run.join(timeout=WAIT)["J07"])
    assert G.stats["acquired"] == 2


def test_light_request_not_blocked_by_held_slot(client, guard_cfg):
    guard_cfg()  # default threshold and slots
    held, release = threading.Event(), threading.Event()

    def hook(ev, g):
        if ev == "acquired":
            held.set()
            assert release.wait(WAIT)

    G.hook = hook
    run, light = Runner(client), Runner(client)
    err = None
    try:
        run.start("P01", url("P01"))
        assert held.wait(WAIT)
        light.start("J01", url("J01"))  # same app, while the only slot is held
        r = light.join(timeout=WAIT)["J01"]  # bounded: a blocked J01 fails here
        assert_golden("J01", r)
        assert G.stats["acquired"] == 1 and G.stats["active"] == 1  # J01 did not take a slot
    except BaseException as e:
        err = e
    finally:
        release.set()
    res = run.join()
    light.join()
    if err is not None:
        raise err
    assert_golden("P01", res["P01"])


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
