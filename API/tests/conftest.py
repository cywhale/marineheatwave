"""pytest setup: isolate Dask, clock and matplotlib before the app is imported.

Strict by default (Codex round 4 #2): a missing data link fails instead of
skipping, and any unexpected distributed.Client construction attempt fails the
test in which it happened.
"""
import os
from pathlib import Path

import pytest

from tests.harness import isolate_env

API_DIR = Path(__file__).resolve().parents[1]
os.chdir(API_DIR)
isolate_env()  # must run before anything imports dask/distributed/matplotlib


@pytest.fixture(scope="session")
def client():
    if not (API_DIR / "data" / "mhw.zarr" / ".zmetadata").exists():
        pytest.fail("API/data/mhw.zarr is missing (spec 16 §2: symlink to the read-only snapshot)")
    from fastapi.testclient import TestClient

    from tests.harness import assert_no_unexpected_client_attempts, has_default_client, patch_app

    app = patch_app()
    with TestClient(app) as c:
        assert not has_default_client()
        yield c
    assert not has_default_client()
    assert_no_unexpected_client_attempts()


@pytest.fixture(autouse=True)
def _no_unexpected_client_attempts():
    yield
    from tests.harness import assert_no_unexpected_client_attempts

    assert_no_unexpected_client_attempts()
