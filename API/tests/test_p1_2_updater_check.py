"""Metadata rules of tools/p1_2/updater_check.py (Codex r20; owner 2026-11-25).

Only the arrays update_mhw.py appends to (sst, sst_anomaly, level, td, date) may
change shape, and only by one along date; lat / lon and any other array must keep
their shape; any other metadata change must be formatting or an explicitly
listed zarr 3 rewrite.  Synthetic v2 metadata only (no data, no zarr needed).
"""
import copy
import json

import pytest

from tools.p1_2 import updater_check as uc

N_OLD = 538


def _arr(shape, chunks, dtype="<f8", **extra):
    return {"chunks": chunks, "compressor": {"id": "blosc", "cname": "lz4", "clevel": 5, "shuffle": 1,
                                             "blocksize": 0},
            "dtype": dtype, "fill_value": "NaN", "filters": None, "order": "C", "shape": shape,
            "zarr_format": 2, **extra}


def _meta(n):
    m = {".zgroup": {"zarr_format": 2}, "anomaly/.zgroup": {"zarr_format": 2},
         "anomaly/date/.zarray": _arr([n], [501], "<i8"), "anomaly/date/.zattrs": {"_ARRAY_DIMENSIONS": ["date"]},
         "anomaly/lat/.zarray": _arr([720], [720], "<f4"), "anomaly/lat/.zattrs": {"_ARRAY_DIMENSIONS": ["lat"]},
         "anomaly/lon/.zarray": _arr([1440], [1440], "<f4"), "anomaly/lon/.zattrs": {"_ARRAY_DIMENSIONS": ["lon"]}}
    for v, dt in (("sst", "<f4"), ("sst_anomaly", "<f8"), ("level", "<i8"), ("td", "<f8")):
        m[f"anomaly/{v}/.zarray"] = _arr([n, 720, 1440], [2, 180, 360], dt)
        m[f"anomaly/{v}/.zattrs"] = {"_ARRAY_DIMENSIONS": ["date", "lat", "lon"]}
    return m


def _write(root, meta, indent=4, extra_root_attrs=False):
    for k, v in meta.items():
        p = root / k
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(v, indent=indent))
    if extra_root_attrs:
        (root / ".zattrs").write_text("{}")
    (root / ".zmetadata").write_text(json.dumps({"metadata": meta, "zarra_consolidated_format": 1}, indent=indent))


def _expected_after():
    """What zarr 3 wrote in V-b: shapes +1, known rewrites, 2-space indentation."""
    m = _meta(N_OLD + 1)
    z = copy.deepcopy(m)
    z["anomaly/.zgroup"]["consolidated_metadata"] = {"kind": "inline", "must_understand": False, "metadata": {}}
    z["anomaly/lat/.zarray"]["dimension_separator"] = "."
    z["anomaly/lon/.zarray"]["dimension_separator"] = "."
    z[".zattrs"] = {}
    return m, z


@pytest.fixture()
def stores(tmp_path):
    before, after = tmp_path / "before", tmp_path / "after"
    _write(before, _meta(N_OLD))

    def make(files, zmeta, root_attrs=True):
        for k, v in files.items():
            p = after / k
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(json.dumps(v, indent=2))
        if root_attrs:
            (after / ".zattrs").write_text("{}")
        (after / ".zmetadata").write_text(json.dumps({"metadata": zmeta, "zarra_consolidated_format": 1}, indent=2))
        return uc.metadata_problems(before, after, N_OLD)

    return make


def test_expected_update_has_no_problem(stores):
    files, zmeta = _expected_after()
    assert stores(files, zmeta) == []


@pytest.mark.parametrize("where", ["zmetadata", "file"])
@pytest.mark.parametrize("arr", ["lat", "lon"])
def test_coordinate_shape_change_fails(stores, where, arr):
    files, zmeta = _expected_after()
    target = zmeta if where == "zmetadata" else files
    target[f"anomaly/{arr}/.zarray"]["shape"] = [721 if arr == "lat" else 1441]
    probs = stores(files, zmeta)
    assert any(p.startswith(f"{where}:anomaly/{arr}/.zarray: shape changed") for p in probs), probs


@pytest.mark.parametrize("where", ["zmetadata", "file"])
def test_unknown_array_shape_change_fails(stores, where, tmp_path):
    # an array the update does not touch (added to both stores first)
    extra = _arr([10, 4], [10, 4])
    (tmp_path / "before" / "anomaly" / "mask").mkdir(parents=True)
    (tmp_path / "before" / "anomaly" / "mask" / ".zarray").write_text(json.dumps(extra))
    z0 = json.loads((tmp_path / "before" / ".zmetadata").read_text())
    z0["metadata"]["anomaly/mask/.zarray"] = extra
    (tmp_path / "before" / ".zmetadata").write_text(json.dumps(z0))
    files, zmeta = _expected_after()
    files["anomaly/mask/.zarray"] = copy.deepcopy(extra)
    zmeta["anomaly/mask/.zarray"] = copy.deepcopy(extra)
    assert stores(files, zmeta) == []
    (files if where == "file" else zmeta)["anomaly/mask/.zarray"]["shape"] = [11, 4]
    probs = stores(files, zmeta)
    assert any(p.startswith(f"{where}:anomaly/mask/.zarray: shape changed") for p in probs), probs


@pytest.mark.parametrize("where", ["zmetadata", "file"])
@pytest.mark.parametrize("bad", [[N_OLD, 720, 1440], [N_OLD + 2, 720, 1440], [N_OLD + 1, 721, 1440]])
def test_appended_array_must_grow_by_exactly_one_date(stores, where, bad):
    files, zmeta = _expected_after()
    (files if where == "file" else zmeta)["anomaly/sst/.zarray"]["shape"] = bad
    probs = stores(files, zmeta)
    assert any(f"{where}:anomaly/sst/.zarray: shape" in p for p in probs), probs


def test_date_array_must_grow(stores):
    files, zmeta = _expected_after()
    zmeta["anomaly/date/.zarray"]["shape"] = [N_OLD]
    assert any("zmetadata:anomaly/date/.zarray: shape" in p for p in stores(files, zmeta))


@pytest.mark.parametrize("key,value", [("chunks", [1, 180, 360]), ("dtype", "<f4"), ("fill_value", 0),
                                       ("compressor", None), ("dimension_separator", "/")])
def test_other_array_keys_of_appended_arrays_must_not_change(stores, key, value):
    files, zmeta = _expected_after()
    files["anomaly/sst_anomaly/.zarray"][key] = value
    assert any(f"file:anomaly/sst_anomaly/.zarray: key {key!r}" in p for p in stores(files, zmeta))


def test_known_rewrite_only_with_its_expected_value(stores):
    files, zmeta = _expected_after()
    zmeta["anomaly/lat/.zarray"]["dimension_separator"] = "/"
    assert any("zmetadata:anomaly/lat/.zarray: key 'dimension_separator'" in p for p in stores(files, zmeta))
    files, zmeta = _expected_after()
    files["anomaly/lat/.zarray"]["dimension_separator"] = "."  # known only inside .zmetadata
    assert any("file:anomaly/lat/.zarray: key 'dimension_separator'" in p for p in stores(files, zmeta))


def test_attribute_change_entry_and_file_changes_fail(stores):
    files, zmeta = _expected_after()
    files["anomaly/td/.zattrs"]["units"] = "days"
    assert any("file:anomaly/td/.zattrs: key 'units'" in p for p in stores(files, zmeta))
    files, zmeta = _expected_after()
    del zmeta["anomaly/lon/.zattrs"]
    assert "zmetadata:anomaly/lon/.zattrs: entry removed" in stores(files, zmeta)
    files, zmeta = _expected_after()
    zmeta["anomaly/extra/.zattrs"] = {}
    assert "zmetadata:anomaly/extra/.zattrs: entry added" in stores(files, zmeta)
    files, zmeta = _expected_after()
    files["anomaly/.zattrs"] = {"x": 1}
    assert "file:anomaly/.zattrs: added" in stores(files, zmeta)


def test_metadata_file_classifier():
    assert uc.is_metadata_file("anomaly/sst/.zarray") and uc.is_metadata_file(".zmetadata")
    assert not uc.is_metadata_file("anomaly/sst/269.0.0") and not uc.is_metadata_file("anomaly/date/1")
