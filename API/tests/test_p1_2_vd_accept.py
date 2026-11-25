"""Pure parts of tools/p1_2/vd_accept.py (V-d production data acceptance)."""
import json
import struct
import zlib

import pytest

from tools.p1_2 import updater_check as uc
from tools.p1_2 import vd_accept as va

C01 = {"size": 15083}
FULL = b"date,level\n" + b"x" * 15072


def _r(status, body=b"", **headers):
    return {"status": status, "body": body, "headers": headers}


def _png(w, h):
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)
    return va.PNG_SIG + struct.pack(">I", 13) + b"IHDR" + ihdr + struct.pack(">I", zlib.crc32(b"IHDR" + ihdr))


def test_range_checks_accept_the_owner_approved_answers():
    rc = va.range_checks(C01, _r(416, **{"content-range": "bytes */15083"}),
                         _r(206, FULL[:10], **{"content-range": "bytes 0-9/15083"}), FULL,
                         _r(400, b"Only support bytes range"))
    assert rc == {"range_416": True, "range_206": True, "range_400": True}


@pytest.mark.parametrize("which,bad", [
    ("416", _r(416, **{"content-range": "bytes */15082"})),
    ("416", _r(200, FULL)),
    ("206", _r(206, FULL[:9], **{"content-range": "bytes 0-9/15083"})),
    ("206", _r(200, FULL)),
    ("400", _r(400, b"other")),
])
def test_range_checks_fail(which, bad):
    good = {"416": _r(416, **{"content-range": "bytes */15083"}),
            "206": _r(206, FULL[:10], **{"content-range": "bytes 0-9/15083"}),
            "400": _r(400, b"Only support bytes range")}
    good[which] = bad
    rc = va.range_checks(C01, good["416"], good["206"], FULL, good["400"])
    assert rc[f"range_{which}"] is False and sum(rc.values()) == 2


def test_png_size():
    assert va.png_size(_png(1200, 800)) == [1200, 800]
    assert va.png_size(b"<html>") is None


@pytest.mark.parametrize("bad,rc", [(False, 0), (True, 3)])
def test_main_writes_json_then_exits_by_ok(tmp_path, monkeypatch, bad, rc):
    def fake_run(base, cases, timeout=120.0):
        out = {"cases": {}, "j08_csv": None, "p01": None, "range": {}}
        return uc.with_checks(out, {"case_J01": True, "j08_csv_equals_57071aa": not bad})
    monkeypatch.setattr(va, "run", fake_run)
    out = tmp_path / "a.json"
    assert va.main(["--base", "https://127.0.0.1:1", "--out", str(out)]) == rc
    doc = json.loads(out.read_text())
    assert uc.result_consistency(doc) == [] and doc["ok"] is (not bad)


def test_default_cases_are_doc27s_seventeen():
    assert len(va.DOC27_CASES) == 17 and "J08" not in va.DOC27_CASES
