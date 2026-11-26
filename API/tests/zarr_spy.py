"""Chunk-read spy for zarr 2 and zarr 3 (P1-2 test infrastructure, doc 29 §4).

Which methods are patched, and why these are all the read paths
----------------------------------------------------------------
zarr 2 (2.18): ``zarr.storage.DirectoryStore.__getitem__``.  ``getitems`` (used
by ``Array._chunk_getitems``) is ``{k: self[k] for k in keys if k in self}``,
so it goes through ``__getitem__`` too; a chunk that is absent on disk is
skipped by the ``in`` test and never reaches ``__getitem__`` (no attempt is
recorded for it: the array returns the fill value).

zarr 3 (3.4.2): ``zarr.storage.LocalStore``

* ``get`` / ``get_partial_values`` (async) -- the default
  ``BatchedCodecPipeline`` and metadata reads.  They run on zarr's own IO
  event-loop thread while the requesting thread blocks.
* ``get_sync`` -- the opt-in ``FusedCodecPipeline`` fast path
  (``codec_pipeline.path`` config) and the sharding codec's sync helpers;
  ``Store.get_ranges_sync`` delegates to it.  The MHW store is v2 format, not
  sharded, and the app uses the default pipeline, so the MHW tests do not go
  through it -- ``test_spy_sees_the_sync_read_path`` proves it is recorded
  anyway, and ``test_mhw_reads_use_the_async_path`` records which path the app
  actually used.

Other ``LocalStore`` methods (``exists``, ``list*``, ``getsize``) return no
chunk bytes.  A store other than ``LocalStore``/``DirectoryStore`` (fsspec,
zip, memory, a wrapper around a non-local store) is NOT covered.

Attempt vs fetch
----------------
A read is recorded twice: as an *attempt* when the store method is called
(before the original read runs, so the order against the caller's hooks is
kept and a read that blocks or fails is still seen), and, once it returns, as
*fetched* (bytes came back) or *missing* (``None`` / ``KeyError``: zarr then
uses the fill value).

* "must NOT read X" checks use the attempts (``reads`` itself): the superset,
  so they are conservative.
* "must read X" checks use ``reads.fetched``: a missing chunk is an attempt
  without data and must not let a positive check pass.

Concurrency
-----------
The patch is on the class, process-wide, and restored in ``finally``.  Only one
spy may be active at a time (a second one raises): the suite runs serially
(pytest-xdist is not installed; under xdist each worker is its own process, so
patches still could not leak between workers).

Keys are store-relative (``anomaly/sst/12.0.3`` for this v2-format store).
Metadata documents (``.zarray``, ``.zattrs``, ``.zgroup``, ``.zmetadata``,
``zarr.json``) are told apart from data chunks: zarr 3 probes ``zarr.json``
paths, which must not count as reading a variable.
"""
from __future__ import annotations

import threading
from contextlib import contextmanager
from typing import Callable

import zarr

ZARR_MAJOR = int(zarr.__version__.split(".")[0])
METADATA_NAMES = {".zarray", ".zattrs", ".zgroup", ".zmetadata", "zarr.json"}


def is_metadata(key: str) -> bool:
    return key.rsplit("/", 1)[-1] in METADATA_NAMES


def data_keys(keys: list[str]) -> list[str]:
    """Data-chunk keys only (metadata documents removed)."""
    return [k for k in keys if not is_metadata(k)]


def variable_of(key: str, group: str = "anomaly") -> str | None:
    """The variable a data-chunk key belongs to (None for metadata or other paths)."""
    parts = key.split("/")
    if is_metadata(key) or len(parts) < 3 or parts[0] != group:
        return None
    return parts[1]


def vars_read(keys: list[str], group: str = "anomaly") -> set[str]:
    """Variables whose data chunks appear in ``keys``."""
    return {v for v in (variable_of(k, group) for k in keys) if v is not None}


def data_var_keys(keys: list[str], exclude: set[str] = frozenset(), group: str = "anomaly") -> list[str]:
    """Data-chunk keys of variables not in ``exclude`` (e.g. the coordinate arrays)."""
    return [k for k in keys if (v := variable_of(k, group)) is not None and v not in exclude]


class Reads(list):
    """Attempted keys, in call order (the list itself), plus outcomes.

    ``fetched`` / ``missing``: keys whose read returned bytes / nothing.
    ``methods``: (method name, key) per attempt -- which store path was used.
    """

    def __init__(self):
        super().__init__()
        self.fetched: list[str] = []
        self.missing: list[str] = []
        self.methods: list[tuple[str, str]] = []
        self._lock = threading.Lock()

    def _attempt(self, method: str, key: str) -> None:
        with self._lock:
            self.append(key)
            self.methods.append((method, key))

    def _outcome(self, key: str, found: bool) -> None:
        with self._lock:
            (self.fetched if found else self.missing).append(key)


_active = threading.Lock()


@contextmanager
def chunk_reads(on_read: Callable[[str], None] | None = None):
    """Record every key read from the local zarr store while the block runs.

    Yields a ``Reads``.  ``on_read(key)`` (optional) is called at the moment of
    the attempt; under zarr 3 it may run on zarr's IO thread, so it must not
    rely on thread-locals of the requesting thread.
    """
    if not _active.acquire(blocking=False):
        raise RuntimeError("chunk_reads() is already active: spies must not overlap")
    reads = Reads()

    def attempt(method: str, key: str) -> None:
        reads._attempt(method, key)
        if on_read is not None:
            on_read(key)

    try:
        if ZARR_MAJOR >= 3:
            from zarr.storage import LocalStore

            orig = {n: getattr(LocalStore, n) for n in ("get", "get_partial_values", "get_sync")}

            async def get(self, key, *args, **kwargs):
                attempt("get", key)
                out = await orig["get"](self, key, *args, **kwargs)
                reads._outcome(key, out is not None)
                return out

            async def get_partial_values(self, prototype, key_ranges):
                key_ranges = list(key_ranges)
                for key, _ in key_ranges:
                    attempt("get_partial_values", key)
                out = await orig["get_partial_values"](self, prototype, key_ranges)
                for (key, _), buf in zip(key_ranges, out):
                    reads._outcome(key, buf is not None)
                return out

            def get_sync(self, key, *args, **kwargs):
                attempt("get_sync", key)
                out = orig["get_sync"](self, key, *args, **kwargs)
                reads._outcome(key, out is not None)
                return out

            patched = {"get": get, "get_partial_values": get_partial_values, "get_sync": get_sync}
            for n, f in patched.items():
                setattr(LocalStore, n, f)
            try:
                yield reads
            finally:
                for n, f in orig.items():
                    setattr(LocalStore, n, f)
        else:
            from zarr.storage import DirectoryStore

            orig2 = DirectoryStore.__getitem__

            def spy(self, key):
                attempt("__getitem__", key)
                try:
                    out = orig2(self, key)
                except KeyError:
                    reads._outcome(key, False)
                    raise
                reads._outcome(key, True)
                return out

            DirectoryStore.__getitem__ = spy
            try:
                yield reads
            finally:
                DirectoryStore.__getitem__ = orig2
    finally:
        _active.release()


def fetched_vars(reads: Reads, exclude: set[str] = frozenset(), group: str = "anomaly"):
    """For exact "these variables were read" checks.

    Returns the set of variables with fetched data chunks when it equals the
    set of variables attempted; otherwise a tuple describing both (so an
    equality assertion against a set fails with a useful message).  A variable
    that was only attempted (every chunk missing) therefore never counts.
    """
    attempted = vars_read(reads, group) - set(exclude)
    fetched = vars_read(reads.fetched, group) - set(exclude)
    return fetched if attempted == fetched else ("attempted", attempted, "fetched", fetched)
