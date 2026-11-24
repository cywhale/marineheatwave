"""Per-worker resource control for heavy requests (P1-1c, spec 17; r2 after Codex round 9).

Heavy work runs in its own thread lane, an AnyIO ``CapacityLimiter`` with
``MHW_HEAVY_SLOTS`` tokens, separate from the default threadpool used by light
requests.  A heavy request waits for a token **in the event loop** and only
gets a thread once it holds one, so queued heavy requests never occupy the
default threadpool (r1 waited on a ``BoundedSemaphore`` inside a threadpool
thread, so a burst of heavy requests could starve light ones: Codex round 9).

Endpoints:
  * JSON / CSV: ``select_mhw_data`` (lazy, reads no chunks) runs in the default
    pool; ``sel.nbytes > MHW_HEAVY_BYTES`` decides the lane; data reading
    (incl. ``xr.concat`` for 0-meridian crossings), serialisation and response
    construction run in that lane via ``run_guarded``.
  * plots are always heavy: the whole plot (data, Figure, savefig) runs in the
    heavy lane.

Lock order: heavy token (event loop) -> thread: read data -> ``plot_lock()`` ->
draw and savefig -> release plot lock -> thread ends -> token released.  No
code running in a thread asks for a token, so the plot lock is never held
while waiting for one.  Waiting is unbounded: no 503, no change to responses.

Settings come from the environment at import time (per worker process):
``MHW_HEAVY_SLOTS`` (default 1) and ``MHW_HEAVY_BYTES`` (default below).
``configure()``, ``stats`` and ``hook`` exist for tests only.
"""
from __future__ import annotations

import asyncio
import os
import threading
from contextlib import contextmanager

import anyio
import anyio.to_thread
from starlette.concurrency import run_in_threadpool

DEFAULT_SLOTS = 1
# P1-0 proposed 4 MiB (18 §4.1).  P1-1c measured 4 concurrent unguarded
# global 1-month sst JSON requests (claim 3.97 MiB each) at a 1,872 MB worker
# peak, above the 1.5 GB observation line, so the default is 2 MiB (doc 23).
DEFAULT_HEAVY_BYTES = 2 * 1024 * 1024

SLOTS = int(os.environ.get("MHW_HEAVY_SLOTS", DEFAULT_SLOTS))
HEAVY_BYTES = int(os.environ.get("MHW_HEAVY_BYTES", DEFAULT_HEAVY_BYTES))

_state_lock = threading.Lock()
_local = threading.local()
_limiter = None        # anyio.CapacityLimiter, created in the worker's event loop
_limiter_loop = None

PLOT_LOCK = threading.Lock()

# Observability for tests: heavy requests waiting for / holding a token.
stats = {"waiting": 0, "active": 0, "max_active": 0, "acquired": 0, "released": 0}
# Optional test hook: hook(event, ticket) with event in
#   "acquiring"   (event loop thread, before waiting for a token; must not block)
#   "acquired"    (heavy-lane thread, token held, before the work starts)
#   "releasing"   (heavy-lane thread, after the work, before the token returns)
#   "plot_locked" (thread holding the plot lock; ticket is None)
# Never set in production.
hook = None


def configure(slots: int | None = None, heavy_bytes: int | None = None) -> None:
    """Replace the settings (tests only; call while no request is running)."""
    global SLOTS, HEAVY_BYTES, _limiter, _limiter_loop
    with _state_lock:
        if slots is not None:
            SLOTS = int(slots)
            _limiter, _limiter_loop = None, None  # rebuilt with the new size on next use
        if heavy_bytes is not None:
            HEAVY_BYTES = int(heavy_bytes)


def reset_stats() -> None:
    with _state_lock:
        stats.update(waiting=0, active=0, max_active=0, acquired=0, released=0)


def is_heavy(nbytes: int) -> bool:
    return nbytes > HEAVY_BYTES


def _call_hook(event, ticket):
    h = hook
    if h is not None:
        h(event, ticket)


def _heavy_limiter():
    """The heavy lane of the running event loop (one per worker process)."""
    global _limiter, _limiter_loop
    loop = asyncio.get_running_loop()
    with _state_lock:
        if _limiter is None or _limiter_loop is not loop:
            _limiter, _limiter_loop = anyio.CapacityLimiter(SLOTS), loop
        return _limiter


class Ticket:
    """Identifies one heavy request in hook events."""

    __slots__ = ()


async def run_guarded(work, *, heavy: bool):
    """Run the synchronous ``work()`` in a thread and return its result.

    Light work uses the default threadpool.  Heavy work first waits for a
    heavy-lane token in the event loop (no thread held while waiting), then
    runs in a thread of the heavy lane; the token is returned when the thread
    finishes, whether ``work`` returns or raises.
    """
    if not heavy:
        return await run_in_threadpool(work)

    limiter = _heavy_limiter()
    ticket = Ticket()
    started = False
    with _state_lock:
        stats["waiting"] += 1
    _call_hook("acquiring", ticket)

    def in_lane():
        nonlocal started
        started = True
        with _state_lock:
            stats["waiting"] -= 1
            stats["active"] += 1
            stats["acquired"] += 1
            stats["max_active"] = max(stats["max_active"], stats["active"])
        try:
            _call_hook("acquired", ticket)
            return work()
        finally:
            try:
                _call_hook("releasing", ticket)
            finally:
                with _state_lock:
                    stats["active"] -= 1
                    stats["released"] += 1

    try:
        # waits for a token without a thread; the thread starts once it has one
        return await anyio.to_thread.run_sync(in_lane, limiter=limiter)
    finally:
        if not started:  # cancelled while still waiting for a token
            with _state_lock:
                stats["waiting"] -= 1


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
