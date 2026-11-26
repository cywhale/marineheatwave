"""The isolation the golden results depend on (spec 16 r2 §2.1)."""
import os

import pytest


def test_no_dask_env_or_client(client):
    import distributed

    from tests.harness import has_default_client

    assert not [k for k in os.environ if k.startswith("DASK_") and k not in ("DASK_CONFIG", "DASK_ROOT_CONFIG")]
    assert not has_default_client()
    from tests.harness import client_attempts, intentional_client_attempt

    before = dict(client_attempts)
    with intentional_client_attempt(), pytest.raises(RuntimeError, match="must not be created"):
        distributed.Client()
    assert client_attempts["intentional"] == before["intentional"] + 1
    assert client_attempts["unexpected"] == before["unexpected"]


def test_swallowed_client_attempt_is_still_detected(client):
    """A manager that catches the refusal must not hide the attempt."""
    import distributed

    from tests import harness

    before = harness.client_attempts["unexpected"]
    try:
        try:
            distributed.Client("tls://127.0.0.1:8786")
        except Exception:
            pass  # what dask_client_manager-style code would do
        assert harness.client_attempts["unexpected"] == before + 1
        with pytest.raises(AssertionError, match="unexpected distributed.Client"):
            harness.assert_no_unexpected_client_attempts()
    finally:
        harness.client_attempts["unexpected"] = before  # restore for the autouse check


def test_synchronous_scheduler_and_agg(client):
    import dask
    import matplotlib

    assert dask.config.get("scheduler") == "synchronous"
    assert matplotlib.get_backend().lower() == "agg"


def test_fixed_clock(client):
    import mhw_app
    import src.mhw_plot
    import src.mhw_utils

    from tests.harness import FIXED_NOW

    for mod in (mhw_app, src.mhw_utils, src.mhw_plot):
        assert mod.datetime.today() == FIXED_NOW
