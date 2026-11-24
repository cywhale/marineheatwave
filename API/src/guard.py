"""Per-worker resource control for heavy requests (P1-1c, spec 17).

Endpoints run in FastAPI's threadpool.  Each endpoint wraps its whole handling
(data loading, JSON/CSV serialisation, plotting and savefig) in::

    with RequestGuard() as g:              # plot endpoints: RequestGuard(heavy=True)
        ... process_mhw_data(..., guard=g) ...

``process_mhw_data`` calls ``g.claim(nbytes)`` with the size of the lazy
selection *before* anything is materialised (for a 0-meridian crossing, the
sum of both subsets before ``xr.concat``).  Only a claim above
``MHW_HEAVY_BYTES`` acquires one of ``MHW_HEAVY_SLOTS`` slots; a guard acquires
at most once and always releases on exit, exceptions included.  Waiting is
unbounded: no 503, no change to responses.

Lock order: heavy slot -> read data -> ``plot_lock()`` -> draw and savefig ->
release plot lock -> release slot.  Acquiring a slot while this thread holds
the plot lock is a lock-order violation and raises immediately.

Settings come from the environment at import time (per worker process):
``MHW_HEAVY_SLOTS`` (default 1) and ``MHW_HEAVY_BYTES`` (default below).
``configure()`` replaces them, for tests.
"""
from __future__ import annotations

import os
import threading
from contextlib import contextmanager

DEFAULT_SLOTS = 1
# P1-0 proposed 4 MiB (18 §4.1).  P1-1c measured 4 concurrent unguarded
# global 1-month sst JSON requests (claim 3.97 MiB each) at a 1,872 MB worker
# peak, above the 1.5 GB observation line, so the default is 2 MiB (doc 23).
DEFAULT_HEAVY_BYTES = 2 * 1024 * 1024

_state_lock = threading.Lock()
_local = threading.local()

SLOTS = int(os.environ.get("MHW_HEAVY_SLOTS", DEFAULT_SLOTS))
HEAVY_BYTES = int(os.environ.get("MHW_HEAVY_BYTES", DEFAULT_HEAVY_BYTES))
_slots = threading.BoundedSemaphore(SLOTS)

PLOT_LOCK = threading.Lock()

# Observability for tests: number of guards holding a slot now / at most.
stats = {"active": 0, "max_active": 0, "acquired": 0, "released": 0}
# Optional test hook: called as hook(event, guard) with event in
# {"acquiring", "acquired", "releasing", "plot_locked"}; never set in production.
hook = None


def configure(slots: int | None = None, heavy_bytes: int | None = None) -> None:
    """Replace the settings (tests only; call while no request is running)."""
    global SLOTS, HEAVY_BYTES, _slots
    with _state_lock:
        if slots is not None:
            SLOTS = int(slots)
            _slots = threading.BoundedSemaphore(SLOTS)
        if heavy_bytes is not None:
            HEAVY_BYTES = int(heavy_bytes)


def reset_stats() -> None:
    with _state_lock:
        stats.update(active=0, max_active=0, acquired=0, released=0)


def _call_hook(event, guard):
    h = hook
    if h is not None:
        h(event, guard)


class RequestGuard:
    """Holds at most one heavy slot for the lifetime of one request."""

    def __init__(self, heavy: bool = False):
        self.heavy = heavy
        self.held = False
        self.claimed_bytes = 0
        self._sem = None

    def __enter__(self) -> "RequestGuard":
        if self.heavy:
            self._acquire()
        return self

    def __exit__(self, *exc) -> bool:
        self.release()
        return False

    def claim(self, nbytes: int) -> None:
        """Declare the size about to be materialised; acquire if it is heavy."""
        self.claimed_bytes += int(nbytes)
        if not self.held and nbytes > HEAVY_BYTES:
            self._acquire()

    def _acquire(self) -> None:
        if getattr(_local, "plot_lock_depth", 0):
            raise RuntimeError("lock order violation: heavy slot requested while holding the plot lock")
        sem = _slots
        _call_hook("acquiring", self)
        sem.acquire()
        self._sem, self.held = sem, True
        with _state_lock:
            stats["active"] += 1
            stats["acquired"] += 1
            stats["max_active"] = max(stats["max_active"], stats["active"])
        _call_hook("acquired", self)

    def release(self) -> None:
        if not self.held:
            return
        try:
            _call_hook("releasing", self)
        finally:
            with _state_lock:
                stats["active"] -= 1
                stats["released"] += 1
            self.held = False
            sem, self._sem = self._sem, None
            sem.release()


class _NoGuard:
    """Used when process_mhw_data is called outside an endpoint (tests, tools)."""

    held = False

    def claim(self, nbytes: int) -> None:
        pass


NO_GUARD = _NoGuard()


@contextmanager
def plot_lock():
    """Serialise Figure creation .. savefig across threads of this worker."""
    with PLOT_LOCK:
        _local.plot_lock_depth = getattr(_local, "plot_lock_depth", 0) + 1
        try:
            _call_hook("plot_locked", None)
            yield
        finally:
            _local.plot_lock_depth -= 1
