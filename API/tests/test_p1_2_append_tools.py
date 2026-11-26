"""Pure parts of the chained-append tools (P1-2, owner 2026-11-25): manifest diff,
point comparison, PNG header, and the run audit (no .env, socket, subprocess,
or write outside the allowed roots)."""
import json
import os
import socket
import struct
import subprocess
import sys
import zlib
from pathlib import Path

import pytest

from tools.p1_2 import append_server as aps
from tools.p1_2 import points as P
from tools.p1_2 import updater_check as uc

API = Path(__file__).resolve().parents[1]


def test_manifest_diff_groups_added_changed_removed_with_sizes():
    before = {"anomaly/sst/268.0.0": ["a", 10], "anomaly/date/1": ["d1", 5], ".zmetadata": ["m1", 100],
              "anomaly/lat/0": ["l", 7], "anomaly/old/0": ["o", 3]}
    after = {"anomaly/sst/268.0.0": ["a", 10], "anomaly/sst/269.0.0": ["n", 20], "anomaly/date/1": ["d2", 6],
             ".zmetadata": ["m2", 110], "anomaly/lat/0": ["l", 7]}
    d = uc.manifest_diff(before, after)
    assert d["summary"] == {"added": {"anomaly/sst": 1}, "changed": {"anomaly/date": 1, ".": 1},
                            "removed": {"anomaly/old": 1}}
    assert d["bytes_written"] == 20 + 6 + 110 and d["unchanged"] == 2
    ch = {e["key"]: e for g in d["changed"].values() for e in g}
    assert ch["anomaly/date/1"]["size_before"] == 5 and ch["anomaly/date/1"]["sha256_before"] == "d1"


def test_manifest_diff_accepts_vb_sha_only_snapshots():
    d = uc.manifest_diff({"x/0": "s1"}, {"x/0": ["s2", 4]})
    assert d["summary"]["changed"] == {"x": 1}


def test_compare_one_row_per_month_with_tolerances():
    want = {"2026-11-01": {"sst": 20.0, "sst_anomaly": 0.5, "level": 1.0, "td": None},
            "2026-12-01": {"sst": None, "sst_anomaly": None, "level": 0.0, "td": None}}
    rows = [{"date": "2026-11-01", "sst": 20.004, "sst_anomaly": 0.5, "level": 1, "td": None},
            {"date": "2026-12-01", "sst": None, "sst_anomaly": None, "level": 0, "td": None}]
    assert P.compare(rows, want)["match"]
    bad = [dict(rows[0], sst_anomaly=0.5001), rows[1]]
    assert not P.compare(bad, want)["match"]
    assert not P.compare(rows[:1], want)["match"]                        # month missing
    assert not P.compare(rows + [rows[1]], want)["match"]                # duplicated month
    assert not P.compare([rows[0], dict(rows[1], td=1.0)], want)["match"]  # NaN expected
    csv_rows = [{"date": "2026-11-01", "sst": "20.003", "sst_anomaly": "0.500004", "level": "1", "td": ""},
                {"date": "2026-12-01", "sst": "", "sst_anomaly": "", "level": "0", "td": ""}]
    assert P.compare(csv_rows, want, default_tol=1e-5)["match"]
    assert not P.compare(csv_rows, want)["match"]


def _png(w, h):
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)
    return b"\x89PNG\r\n\x1a\n" + struct.pack(">I", 13) + b"IHDR" + ihdr + struct.pack(">I", zlib.crc32(b"IHDR" + ihdr))


def test_png_info():
    assert aps.png_info(_png(1200, 800)) == {"png_signature": True, "width": 1200, "height": 800}
    assert aps.png_info(b"<html>")["png_signature"] is False


_AUDIT_CHILD = r"""
import json, os, socket, subprocess, sys, tempfile
sys.path.insert(0, sys.argv[1])
from tools.p1_2 import updater_check as uc
store, cwd, outside = sys.argv[2], sys.argv[3], sys.argv[4]
events, hook = uc._audit_recorder(store, [store, cwd])
sys.addaudithook(hook)
open(os.path.join(store, "a"), "w").close()          # allowed (store)
open(os.path.join(cwd, "b"), "w").close()            # allowed (cwd)
open(os.path.join(cwd, ".env"), "w").close()         # a .env open is always reported
open(os.path.join(outside, "c"), "w").close()        # write outside
try:
    socket.create_connection(("192.0.2.1", 9), timeout=0.01)
except OSError:
    pass
subprocess.run(["true"])
print(json.dumps(events))
"""


def test_run_audit_records_dotenv_socket_subprocess_and_outside_writes(tmp_path):
    store, cwd, outside = (tmp_path / n for n in ("store", "cwd", "outside"))
    for d in (store, cwd, outside):
        d.mkdir()
    out = subprocess.run([sys.executable, "-c", _AUDIT_CHILD, str(API), str(store), str(cwd), str(outside)],
                         capture_output=True, text=True, check=True).stdout
    ev = json.loads(out.strip().splitlines()[-1])
    assert ev["dotenv_file_open"] == [str(cwd / ".env")]
    assert any(e.startswith("socket.") for e in ev["socket"])
    assert ev["subprocess"]
    assert ev["write_outside"] == [str(outside / "c")] and ev["write_opens_in_store"] == 1


# ---- Codex r22 (Minor): the response's month list must be exactly the queried range

_WANT2 = {"2026-11-01": {"sst": 20.0, "sst_anomaly": 0.5, "level": 1.0, "td": None},
          "2026-12-01": {"sst": 20.1, "sst_anomaly": 0.625, "level": 2.0, "td": None}}


def _row(month, **kw):
    w = _WANT2.get(month, _WANT2["2026-11-01"])
    return {"date": month, **w, **kw}


def test_compare_exact_month_list_passes():
    res = P.compare([_row("2026-11-01"), _row("2026-12-01")], _WANT2)
    assert res["match"] and res["date_list_match"] and res["n_rows"] == 2 == res["n_want"]
    assert res["dates_got"] == res["dates_want"] == ["2026-11-01", "2026-12-01"]


@pytest.mark.parametrize("months", [
    ["2026-10-01", "2026-11-01", "2026-12-01"],   # a month before the queried range
    ["2026-11-01", "2026-12-01", "2027-01-01"],   # a month after it
    ["2026-12-01", "2026-11-01"],                 # wrong order
    ["2026-11-01", "2026-11-01", "2026-12-01"],   # duplicate
    ["2026-11-01"],                               # missing month
    [],                                           # no rows
    ["2026-11-01", "2026-12-01", ""],             # a row without a date
])
def test_compare_rejects_a_wrong_month_list(months):
    rows = [_row(m) if m else {"date": None} for m in months]
    res = P.compare(rows, _WANT2)
    assert not res["match"] and not res["date_list_match"] and res["dates_got"] == [m or "None" for m in months]


def test_compare_extra_month_fails_even_when_every_wanted_month_matches():
    rows = [_row("2026-11-01"), _row("2026-12-01"), _row("2027-01-01")]
    res = P.compare(rows, _WANT2)
    assert all(m["match"] for m in res["months"].values()) and not res["match"]


def test_compare_csv_month_list_and_tolerances_kept():
    csv_rows = [{"date": "2026-11-01", "sst": "20.003", "sst_anomaly": "0.500004", "level": "1", "td": ""},
                {"date": "2026-12-01", "sst": "20.1", "sst_anomaly": "0.625", "level": "2", "td": ""}]
    assert P.compare(csv_rows, _WANT2, default_tol=1e-5)["match"]
    assert not P.compare(csv_rows[::-1], _WANT2, default_tol=1e-5)["match"]
    assert not P.compare([dict(csv_rows[0], td="1.0"), csv_rows[1]], _WANT2, default_tol=1e-5)["match"]  # NaN kept
