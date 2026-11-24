"""The golden guard really fails on gaps (Codex round 4 #2). No data needed."""
import copy

import pytest

from tests.golden.cases import CASES
from tests.golden.runner import (ManifestError, header_problems, load_manifest, stability_problems,
                                 validate_manifest)

CSV_H = {"content-type": "text/csv; charset=utf-8", "content-disposition": "attachment; filename=\"x.csv\"",
         "accept-ranges": "bytes", "content-length": "10"}
JSON_H = {"content-type": "application/json", "content-disposition": None}


def _good_manifest():
    cases = {}
    for c in CASES:
        h = CSV_H if c["kind"] == "csv" and c["expect"]["status"] == 200 else JSON_H
        cases[c["id"]] = {"url": c["url"], "kind": c["kind"], "status": c["expect"]["status"],
                          "headers": dict(h), "sha256": "x", "stable": True}
    return {"cases": cases, "approved_unstable": []}


def test_good_manifest_passes():
    assert validate_manifest(_good_manifest(), CASES) == []


def test_missing_manifest_file_fails(tmp_path):
    with pytest.raises(ManifestError, match="missing"):
        load_manifest(tmp_path / "nope.json")


def test_case_missing_from_manifest_fails():
    m = _good_manifest()
    del m["cases"]["J07"]
    assert any("J07: in CASES but missing" in p for p in validate_manifest(m, CASES))


def test_extra_manifest_case_fails():
    m = _good_manifest()
    m["cases"]["ZZ9"] = copy.deepcopy(m["cases"]["J01"])
    assert any("ZZ9: in manifest but not in CASES" in p for p in validate_manifest(m, CASES))


def test_url_drift_fails():
    m = _good_manifest()
    m["cases"]["J01"]["url"] += "&append=sst"
    assert any("J01: URL differs" in p for p in validate_manifest(m, CASES))


def test_unapproved_unstable_fails_and_approved_passes():
    m = _good_manifest()
    m["cases"]["P01"]["stable"] = False
    assert any("P01: recorded as unstable" in p for p in validate_manifest(m, CASES))
    m["approved_unstable"] = ["P01"]
    assert validate_manifest(m, CASES) == []


def test_missing_csv_header_evidence_fails():
    m = _good_manifest()
    del m["cases"]["C01"]["headers"]["content-length"]
    assert any("C01: manifest lacks header evidence for 'content-length'" in p for p in validate_manifest(m, CASES))


@pytest.mark.parametrize("field,value,needle", [
    ("status", 500, "status"),
    ("headers", {"content-type": "text/plain"}, "headers"),
    ("sha256", "different", "sha256"),
])
def test_stability_checks_status_headers_and_body(field, value, needle):
    base = {"status": 200, "headers": {"content-type": "application/json"}, "sha256": "a"}
    other = dict(base, **{field: value})
    assert stability_problems([base, dict(base), dict(base)]) == []
    assert any(needle in p for p in stability_problems([base, other, dict(base)]))


def test_csv_content_length_must_match_body():
    case = {"kind": "csv"}
    ok = {"status": 200, "size": 10, "headers": dict(CSV_H)}
    assert header_problems(case, ok) == []
    bad = {"status": 200, "size": 11, "headers": dict(CSV_H)}
    assert any("content-length 10 != body length 11" in p for p in header_problems(case, bad))
    no_ranges = {"status": 200, "size": 10, "headers": {k: v for k, v in CSV_H.items() if k != "accept-ranges"}}
    assert header_problems(case, no_ranges) == ["CSV response has no accept-ranges"]
