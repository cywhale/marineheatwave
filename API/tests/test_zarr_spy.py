"""Self-tests of tests/zarr_spy.py on small synthetic stores (zarr 2 and zarr 3).

These do not touch the MHW store; the MHW-store self-check is
``test_p1_1a.test_spy_records_a_small_slice_and_nothing_for_lazy_inspection``.
"""
import numpy as np
import pytest
import zarr

from tests.zarr_spy import ZARR_MAJOR, chunk_reads, vars_read

DATA = np.arange(4, dtype="f4").reshape(2, 2) + 10


def _make_store(path, compressor_none=False):
    """anomaly/x: shape (4, 4), chunks (2, 2), fill NaN; only chunk 0.0 written."""
    kw = dict(shape=(4, 4), chunks=(2, 2), dtype="f4", fill_value=np.nan)
    if ZARR_MAJOR >= 3:
        g = zarr.open_group(str(path), mode="w", zarr_format=2)
        if compressor_none:
            kw["compressors"] = None
        a = g.require_group("anomaly").create_array("x", **kw)
    else:
        g = zarr.open_group(str(path), mode="w")
        if compressor_none:
            kw["compressor"] = None
        a = g.require_group("anomaly").create_dataset("x", **kw)
    a[0:2, 0:2] = DATA
    return path


def _open(path):
    return zarr.open_group(str(path), mode="r")["anomaly/x"]


def test_present_chunk_is_fetched_with_its_data(tmp_path):
    a = _open(_make_store(tmp_path / "s.zarr"))
    with chunk_reads() as reads:
        v = a[0:2, 0:2]
    np.testing.assert_array_equal(v, DATA)
    assert "anomaly/x/0.0" in reads.fetched
    assert reads.missing == []
    assert vars_read(reads.fetched) == {"x"}


def test_missing_chunk_cannot_pass_a_positive_check(tmp_path):
    """Chunk 1.1 is not on disk: zarr returns the fill value.  The attempt (if
    any) must not count as a fetch, so a positive "x was read" check on
    ``reads.fetched`` fails, as it should."""
    a = _open(_make_store(tmp_path / "s.zarr"))
    with chunk_reads() as reads:
        v = a[2:4, 2:4]
    assert np.isnan(v).all()
    assert reads.fetched == []
    assert vars_read(reads.fetched) == set()  # a positive check would fail here
    if ZARR_MAJOR >= 3:
        # zarr 3 calls get(), which returns None: an attempt without data
        assert list(reads) == reads.missing == ["anomaly/x/1.1"]
        assert vars_read(reads) == {"x"}  # an attempt-based positive check would pass: wrong basis
    else:
        # zarr 2 skips absent keys before __getitem__ (getitems: `if k in self`)
        assert list(reads) == [] and reads.missing == []


@pytest.mark.skipif(ZARR_MAJOR < 3, reason="the sync read path exists in zarr 3 only")
def test_spy_sees_the_sync_read_path(tmp_path):
    """FusedCodecPipeline reads through LocalStore.get_sync; the spy records it."""
    path = _make_store(tmp_path / "s.zarr", compressor_none=True)
    with zarr.config.set({"codec_pipeline.path": "zarr.core.codec_pipeline.FusedCodecPipeline"}):
        a = _open(path)
        with chunk_reads() as reads:
            v = a[0:2, 0:2]
    np.testing.assert_array_equal(v, DATA)
    assert ("get_sync", "anomaly/x/0.0") in reads.methods, reads.methods
    assert "anomaly/x/0.0" in reads.fetched


def test_spies_cannot_overlap_and_patches_are_restored(tmp_path):
    if ZARR_MAJOR >= 3:
        from zarr.storage import LocalStore as S
        names = ("get", "get_partial_values", "get_sync")
    else:
        from zarr.storage import DirectoryStore as S
        names = ("__getitem__",)
    before = {n: S.__dict__[n] for n in names}

    with pytest.raises(ZeroDivisionError):
        with chunk_reads():
            with pytest.raises(RuntimeError, match="already active"):
                with chunk_reads():
                    pass
            1 / 0
    assert {n: S.__dict__[n] for n in names} == before

    a = _open(_make_store(tmp_path / "s.zarr"))
    with chunk_reads() as reads:  # usable again after the error
        a[0:2, 0:2]
    assert reads.fetched
