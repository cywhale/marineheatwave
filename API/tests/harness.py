"""Isolation shared by the golden tests and the benchmark app wrapper.

Call ``isolate_env()`` before anything imports dask, distributed or matplotlib,
then ``patch_app()`` to import the MHW app with:

* no distributed Dask client (env, config files and client creation all blocked),
* a synchronous dask scheduler,
* a fixed clock (2026-11-20 12:00) in every module that reads the time,
* the Agg matplotlib backend (the macOS GUI backend aborts in worker threads).

The production app is not modified; everything here is monkeypatching.
"""
from __future__ import annotations

import datetime as _dt
import os
import tempfile
import threading
from contextlib import contextmanager

FIXED_NOW = _dt.datetime(2026, 11, 20, 12, 0, 0)

# Modules that read the clock in 0.2.8: mhw_app.py:240, src/mhw_utils.py:79,
# src/mhw_plot.py:49 and :78 (all via a module-level ``datetime`` name).
CLOCK_MODULES = ("mhw_app", "src.mhw_utils", "src.mhw_plot")
DASK_CLIENT_MODULES = ("mhw_app", "src.mhw_utils")

dask_client_calls = {"count": 0}
# Every attempt to construct distributed.Client is counted, even if the caller
# (e.g. a resilient client manager) swallows the RuntimeError we raise.
# Attempts made on purpose by the isolation tests are counted separately.
client_attempts = {"unexpected": 0, "intentional": 0}
_intentional = threading.local()


@contextmanager
def intentional_client_attempt():
    """Mark a deliberate Client() call (used only by the isolation tests)."""
    _intentional.on = True
    try:
        yield
    finally:
        _intentional.on = False


def assert_no_unexpected_client_attempts() -> None:
    if client_attempts["unexpected"]:
        raise AssertionError(f"{client_attempts['unexpected']} unexpected distributed.Client construction attempt(s)")


class FixedDateTime(_dt.datetime):
    """datetime subclass whose today()/now() return FIXED_NOW."""

    @classmethod
    def today(cls):
        f = FIXED_NOW
        return cls(f.year, f.month, f.day, f.hour, f.minute, f.second)

    @classmethod
    def now(cls, tz=None):
        if tz is not None:
            raise NotImplementedError("0.2.8 never asks for an aware now(); add support if that changes")
        return cls.today()


def isolate_env() -> str:
    """Drop DASK_* variables, point dask config at an empty dir, force Agg."""
    for key in [k for k in os.environ if k.startswith("DASK_")]:
        del os.environ[key]
    empty = tempfile.mkdtemp(prefix="mhw-dask-config-")
    os.environ["DASK_CONFIG"] = empty
    os.environ["DASK_ROOT_CONFIG"] = empty
    os.environ["MPLBACKEND"] = "Agg"
    return empty


def _no_dask_client(*_args, **_kwargs):
    dask_client_calls["count"] += 1
    return None


def _refuse_client(self, *_args, **_kwargs):
    key = "intentional" if getattr(_intentional, "on", False) else "unexpected"
    client_attempts[key] += 1
    raise RuntimeError("distributed.Client must not be created under the MHW test harness")


def patch_app(fixed_clock: bool = True):
    """Import mhw_app with Dask and the clock isolated; return the FastAPI app."""
    import importlib

    import dask
    import distributed

    dask.config.set(scheduler="synchronous")
    distributed.Client.__init__ = _refuse_client

    app_module = importlib.import_module("mhw_app")
    for name in DASK_CLIENT_MODULES:
        mod = importlib.import_module(name)
        if hasattr(mod, "get_dask_client"):
            mod.get_dask_client = _no_dask_client
    if fixed_clock:
        for name in CLOCK_MODULES:
            mod = importlib.import_module(name)
            if hasattr(mod, "datetime"):
                mod.datetime = FixedDateTime
    return app_module.app


def has_default_client() -> bool:
    import distributed

    try:
        distributed.default_client()
    except ValueError:
        return False
    return True


def describe_env() -> dict:
    """Facts worth recording next to any golden or benchmark result."""
    import importlib
    import platform

    import matplotlib

    versions = {}
    for pkg in ("xarray", "zarr", "numpy", "pandas", "polars", "dask", "distributed",
                "fastapi", "starlette", "matplotlib", "numcodecs"):
        try:
            versions[pkg] = importlib.import_module(pkg).__version__
        except Exception as exc:  # pragma: no cover - informational only
            versions[pkg] = f"unavailable: {exc}"
    import dask

    return {
        "pid": os.getpid(),
        "python": platform.python_version(),
        "machine": platform.machine(),
        "platform": platform.platform(),
        "versions": versions,
        "fixed_now": FIXED_NOW.isoformat(),
        "mpl_backend": matplotlib.get_backend(),
        "dask_scheduler": dask.config.get("scheduler", None),
        "dask_env": sorted(k for k in os.environ if k.startswith("DASK_") and k not in ("DASK_CONFIG", "DASK_ROOT_CONFIG")),
        "default_client": has_default_client(),
        "client_attempts": dict(client_attempts),
        "zarr_path": os.path.realpath("data/mhw.zarr"),
    }
