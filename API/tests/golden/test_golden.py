"""Byte comparison against tests/golden/manifest.json (spec 16 r2 §3).

Strict regression mode (Codex round 4 #2): a missing manifest, a case missing
from the manifest, CASES/manifest drift (ids, URLs, kinds) or an unstable case
that the owner has not listed in ``approved_unstable`` all FAIL; nothing is
skipped or xfailed.  A mismatching body is written to
$MHW_GOLDEN_DIR/_actual/<case>.<ext> for diffing.  Never edit the manifest to
make a test pass.

PNG cases (P01-P06): status, headers, the PNG signature and the image size are
checked, but not the bytes -- PNG output differs between matplotlib versions
and is not part of the acceptance (owner decision 2026-11-25, doc 29); the data
behind the plots is checked by tests/test_p1_1d.py.  JSON / CSV / error cases
stay byte-compared.
"""
import struct

import pytest

from tests.golden.cases import CASES
from tests.golden.runner import (GOLDEN_ROOT, ManifestError, body_path, expect_problems, fetch,
                                 header_problems, load_manifest, validate_manifest)

try:
    MANIFEST_DATA = load_manifest()
    MANIFEST_ERROR = None
except ManifestError as exc:
    MANIFEST_DATA, MANIFEST_ERROR = None, str(exc)


def _manifest():
    if MANIFEST_DATA is None:
        pytest.fail(MANIFEST_ERROR)
    return MANIFEST_DATA


def test_manifest_matches_cases():
    assert validate_manifest(_manifest(), CASES) == []


@pytest.mark.parametrize("case", CASES, ids=[c["id"] for c in CASES])
def test_golden(client, case):
    manifest = _manifest()
    golden = manifest["cases"].get(case["id"])
    if golden is None:
        pytest.fail(f"{case['id']} is missing from the manifest")
    if golden["url"] != case["url"]:
        pytest.fail(f"{case['id']}: CASES URL differs from the recorded one")
    if not golden["stable"] and case["id"] not in manifest.get("approved_unstable", []):
        pytest.fail(f"{case['id']} was not byte-stable when recorded and is not approved (spec 16 §3.3)")
    res = fetch(client, case)
    assert expect_problems(case, res) == []
    assert header_problems(case, res) == []
    assert res["status"] == golden["status"]
    assert res["headers"] == golden["headers"]
    if case["kind"] == "png":
        assert res["body"][:8] == b"\x89PNG\r\n\x1a\n", f"{case['id']} is not a PNG"
        size = struct.unpack(">II", res["body"][16:24])
        assert size == (golden["summary"]["width"], golden["summary"]["height"]), size
        return
    if res["sha256"] != golden["sha256"]:
        out = GOLDEN_ROOT / "_actual"
        out.mkdir(parents=True, exist_ok=True)
        body_path(out, case, res["status"]).write_bytes(res["body"])
        pytest.fail(f"{case['id']} body differs from golden ({res['size']} vs {golden['size']} bytes); "
                    f"actual saved under {out}")
