"""Golden query set (spec 16 r2 §3.1).

Every URL fixes its dates so results only depend on the zarr snapshot and the
fixed clock in tests/harness.py.  ``expect`` holds facts checked before the
byte comparison: the branch each case is meant to exercise.  ``rows`` is only
set where it was verified against the 2026-10 snapshot.
"""
from __future__ import annotations

FULL = "start=1982-01-01&end=2026-10-01"

CASES: list[dict] = [
    # --- JSON: points -------------------------------------------------------
    {"id": "J01", "kind": "json", "url": f"/api/mhw?lon0=135&lat0=15&{FULL}",
     "expect": {"status": 200, "rows": 538, "first": {"lon": 135.125, "lat": 15.125},
                "columns": ["lon", "lat", "date", "level"]},
     "why": "point branch, default level"},
    {"id": "J02", "kind": "json", "url": f"/api/mhw?lon0=135&lat0=15&{FULL}&append=sst,sst_anomaly,level,td",
     "expect": {"status": 200, "rows": 538,
                "columns": ["lon", "lat", "date", "level", "sst", "sst_anomaly", "td"]},
     "why": "4 variables, sorted columns, td NaN -> null"},
    {"id": "J03", "kind": "json", "url": f"/api/mhw?lon0=135&lat0=15&{FULL}&append=td,%20sst,td",
     "expect": {"status": 200, "rows": 538, "columns": ["lon", "lat", "date", "sst", "td"]},
     "why": "duplicates, whitespace, sort"},
    {"id": "J04", "kind": "json", "url": f"/api/mhw?lon0=100&lat0=40&{FULL}",
     "expect": {"status": 200, "rows": 538, "first": {"lon": 100.125, "lat": 40.125}},
     "why": "land point: sst/anomaly/td all NaN, level all 0"},
    {"id": "J04b", "kind": "json", "url": f"/api/mhw?lon0=100&lat0=40&{FULL}&append=sst",
     "expect": {"status": 400, "detail": "No data available after removing rows with NA values."},
     "why": "land + non-level variable: dropna removes every row"},
    {"id": "J05", "kind": "json", "url": f"/api/mhw?lon0=-1.875&lat0=-60.875&{FULL}&append=sst,sst_anomaly,level",
     "expect": {"status": 200, "rows": 538, "first": {"lon": -1.875, "lat": -60.875}},
     "why": "sea ice: anomaly NaN in 358 months, level -1 in 280"},
    {"id": "J05b", "kind": "json", "url": f"/api/mhw?lon0=-60&lat0=-65&{FULL}&append=sst_anomaly,level",
     "expect": {"status": 200, "rows": 538, "first": {"lon": -59.875, "lat": -64.875}},
     "why": "anomaly all NaN, level mixes -1 (60) and 0 (478)"},
    {"id": "J06a", "kind": "json", "url": "/api/mhw?lon0=-180&lat0=0&start=2020-01-01&end=2020-12-01",
     "expect": {"status": 200, "rows": 12, "first": {"lon": -179.875, "lat": 0.125}}, "why": "lon clamp"},
    {"id": "J06b", "kind": "json", "url": "/api/mhw?lon0=180&lat0=0&start=2020-01-01&end=2020-12-01",
     "expect": {"status": 200, "rows": 12, "first": {"lon": -179.875, "lat": 0.125}}, "why": "lon clamp"},
    {"id": "J06c", "kind": "json", "url": "/api/mhw?lon0=-179.875&lat0=0&start=2020-01-01&end=2020-12-01",
     "expect": {"status": 200, "rows": 12, "first": {"lon": -179.875, "lat": 0.125}}, "why": "lon clamp"},
    {"id": "J06d", "kind": "json", "url": "/api/mhw?lon0=179.9&lat0=0&start=2020-01-01&end=2020-12-01",
     "expect": {"status": 200, "rows": 12, "first": {"lon": 179.875, "lat": 0.125}}, "why": "lon clamp"},
    {"id": "J13a", "kind": "json", "url": f"/api/mhw?lon0=0&lat0=90&{FULL}&append=sst,level",
     "expect": {"status": 200, "rows": 538, "first": {"lon": 0.125, "lat": 89.875}}, "why": "north edge"},
    {"id": "J13b", "kind": "json", "url": f"/api/mhw?lon0=0&lat0=-90&{FULL}&append=sst,level",
     "expect": {"status": 200, "rows": 538, "first": {"lon": 0.125, "lat": -89.875}},
     "why": "south edge: sst all NaN, level all 0"},
    # --- JSON: bounding boxes ----------------------------------------------
    {"id": "J07", "kind": "json", "url": "/api/mhw?lon0=-5&lon1=5&lat0=-5&lat1=5&start=2020-01-01&end=2021-12-01&append=sst",
     "expect": {"status": 200, "months": 24}, "why": "crosses 0 deg: concat branch, no truncation"},
    {"id": "J08", "kind": "json", "url": "/api/mhw?lon0=170&lon1=-170&lat0=-5&lat1=5&start=2020-01-01&end=2021-12-01&append=sst",
     "expect": {"status": 200, "months": 12}, "why": "crosses 180: lon_range computed as 340 -> truncated to 1 year"},
    {"id": "J09", "kind": "json", "url": "/api/mhw?lon0=130&lon1=140&lat0=10&lat1=20&start=2010-01-01&end=2019-12-01&append=sst_anomaly",
     "expect": {"status": 200, "months": 120}, "why": "10x10 deg, 10 years, no truncation"},
    {"id": "J10", "kind": "json", "url": "/api/mhw?lon0=135&lon1=150&lat0=15&lat1=30&start=2013-01-01",
     "expect": {"status": 200, "months": 12}, "why": "no end -> fixed clock; area 225 -> 1 year"},
    {"id": "J11", "kind": "json", "url": "/api/mhw?lon0=-180&lon1=180&lat0=-90&lat1=90&start=2026-10-01&end=2026-10-01",
     "expect": {"status": 200, "months": 1}, "why": "global, one month, concat branch, ~61 MB"},
    {"id": "J12", "kind": "json", "url": "/api/mhw?lon0=135&lon1=136&lat0=15&lat1=15&start=2020-01-01&end=2020-12-01",
     "expect": {"status": 200, "months": 12}, "why": "narrower than one cell: offset_lat branch"},
    # --- JSON: errors and empty results --------------------------------------
    {"id": "E01", "kind": "json", "url": "/api/mhw?lon0=135&lat0=15&append=foo",
     "expect": {"status": 400, "detail_prefix": "Invalid variable(s)"}, "why": "bad parameter"},
    {"id": "E02", "kind": "json", "url": "/api/mhw?lon0=100&lon1=101&lat0=40&lat1=41&start=2020-01-01&end=2020-03-01",
     "expect": {"status": 200, "all_level": 0}, "why": "land bbox + default level is NOT empty"},
    {"id": "E02b", "kind": "json", "url": "/api/mhw?lon0=100&lon1=101&lat0=40&lat1=41&start=2020-01-01&end=2020-03-01&append=sst",
     "expect": {"status": 400, "detail": "No data available after removing rows with NA values."},
     "why": "real empty result"},
    {"id": "E03", "kind": "json", "url": "/api/mhw?lon0=135&lat0=15&start=2030-01-01&end=2030-02-01",
     "expect": {"status": 400, "detail": "No data available after removing rows with NA values."},
     "why": "dates beyond data; nbytes==0 branch never fires (coords count)"},
    # --- CSV ---------------------------------------------------------------
    {"id": "C01", "kind": "csv", "url": f"/api/mhw/csv?lon0=135&lat0=15&{FULL}", "expect": {"status": 200}, "why": "as J01"},
    {"id": "C02", "kind": "csv", "url": f"/api/mhw/csv?lon0=135&lat0=15&{FULL}&append=sst,sst_anomaly,level,td",
     "expect": {"status": 200}, "why": "as J02"},
    {"id": "C03", "kind": "csv", "url": "/api/mhw/csv?lon0=-5&lon1=5&lat0=-5&lat1=5&start=2020-01-01&end=2021-12-01&append=sst",
     "expect": {"status": 200}, "why": "as J07"},
    {"id": "C04", "kind": "csv", "url": "/api/mhw/csv?lon0=100&lon1=101&lat0=40&lat1=41&start=2020-01-01&end=2020-03-01&append=sst",
     "expect": {"status": 400}, "why": "CSV route returns a JSON 400 body"},
    {"id": "C05", "kind": "csv", "url": "/api/mhw/csv?lon0=-180&lat0=0&start=2020-01-01&end=2020-12-01",
     "expect": {"status": 200}, "why": "same cell as C06, different filename (deg2str oriNeg)"},
    {"id": "C06", "kind": "csv", "url": "/api/mhw/csv?lon0=180&lat0=0&start=2020-01-01&end=2020-12-01",
     "expect": {"status": 200}, "why": "same cell as C05"},
    # --- PNG ---------------------------------------------------------------
    {"id": "P01", "kind": "png", "url": "/api/mhw/plot?bbox=(135,15,137,17)&period=1982-2011,2012-2021,2022&sstvar=sst&mode=month",
     "expect": {"status": 200}, "why": "month mode"},
    {"id": "P02", "kind": "png", "url": "/api/mhw/plot?bbox=(100,40)&period=1982-2011,2012-2021&sstvar=sst&mode=month",
     "expect": {"status": 200}, "why": "land point: empty lines, still 200"},
    {"id": "P03", "kind": "png", "url": "/api/mhw/plot?bbox=(-90,-10,-80,0),(-150,-5,-90,5)&period=2015-202306&sstvar=sst_anomaly",
     "expect": {"status": 200}, "why": "series mode, two bboxes"},
    {"id": "P04", "kind": "png", "url": "/api/mhw/plot?bbox=(20,55,21,56)&period=2015-202306&sstvar=sst_anomaly",
     "expect": {"status": 200}, "why": "18 interior NaN months -> line gaps"},
    {"id": "P05", "kind": "png", "url": "/api/mhw/plot?bbox=(135,15),(140,20)&period=all&palette=blue,green",
     "expect": {"status": 200}, "why": "period=all, palette, point bboxes"},
    {"id": "P06", "kind": "png", "url": "/api/mhw/plot?bbox=(135,15,137,17)&period=2025,2026&sstvar=sst_anomaly&mode=month",
     "expect": {"status": 200}, "why": "2026 == fixed clock year -> red line"},
]

BY_ID = {c["id"]: c for c in CASES}
