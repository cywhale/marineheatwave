"""Chunk-read spy for zarr 2 and zarr 3 (P1-2 test infrastructure, doc 29 §4).

zarr 2 reads a directory store synchronously through
``zarr.storage.DirectoryStore.__getitem__``.  zarr 3 reads through the async
``zarr.storage.LocalStore.get`` / ``get_partial_values``, which run on zarr's
own IO event-loop thread while the calling thread blocks; the patched
coroutines record the key when they are called, before the original read is
awaited, so a read that happens is always recorded and ordering against the
caller (e.g. the guard's "acquired" hook) is preserved.

Keys are store-relative (``anomaly/sst/12.0.3`` for this v2-format store).
Metadata documents (``.zarray``, ``.zattrs``, ``.zgroup``, ``.zmetadata``,
``zarr.json``) are told apart from data chunks: zarr 3 probes ``zarr.json``
paths, which must not count as reading a variable.
"""
from __future__ import annotations

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
    """Variables whose data chunks were read."""
    return {v for v in (variable_of(k, group) for k in keys) if v is not None}


def data_var_keys(keys: list[str], exclude: set[str] = frozenset(), group: str = "anomaly") -> list[str]:
    """Data-chunk keys of variables not in ``exclude`` (e.g. the coordinate arrays)."""
    return [k for k in keys if (v := variable_of(k, group)) is not None and v not in exclude]


@contextmanager
def chunk_reads(on_read: Callable[[str], None] | None = None):
    """Record every key read from the local zarr store while the block runs.

    ``on_read(key)`` (optional) is called at the moment of the read; under
    zarr 3 it runs on zarr's IO thread, so it must not rely on thread-locals of
    the requesting thread.
    """
    keys: list[str] = []

    def note(key: str) -> None:
        keys.append(key)
        if on_read is not None:
            on_read(key)

    if ZARR_MAJOR >= 3:
        from zarr.storage import LocalStore

        orig_get, orig_partial = LocalStore.get, LocalStore.get_partial_values

        async def get(self, key, *args, **kwargs):
            note(key)
            return await orig_get(self, key, *args, **kwargs)

        async def get_partial_values(self, prototype, key_ranges):
            key_ranges = list(key_ranges)
            for key, _ in key_ranges:
                note(key)
            return await orig_partial(self, prototype, key_ranges)

        LocalStore.get, LocalStore.get_partial_values = get, get_partial_values
        try:
            yield keys
        finally:
            LocalStore.get, LocalStore.get_partial_values = orig_get, orig_partial
    else:
        from zarr.storage import DirectoryStore

        orig = DirectoryStore.__getitem__

        def spy(self, key):
            note(key)
            return orig(self, key)

        DirectoryStore.__getitem__ = spy
        try:
            yield keys
        finally:
            DirectoryStore.__getitem__ = orig
