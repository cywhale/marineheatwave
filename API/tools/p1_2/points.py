"""Sample points and the comparison of API rows with fixture expectations
(P1-2 updater checks).  Pure Python (no numpy), so the loopback server check can
run in the tools venv; ``updater_check expect-points`` writes the expectations.
"""
from __future__ import annotations

SAMPLES = [(-14.875, 150.125), (0.125, 200.125), (35.125, 140.125), (-60.125, 300.125), (10.125, 125.125)]
FIELDS = ("sst", "sst_anomaly", "level", "td")
TOL = {"sst": 0.006}  # int16 * 0.01 quantisation of the NOAA fixture; everything else exact
APPEND = "sst,sst_anomaly,level,td"


def grid_index(lat: float, lon: float) -> tuple[int, int]:
    return int(round((lat + 89.875) / 0.25)), int(round((lon - 0.125) / 0.25))


def query_lon(lon: float) -> float:
    return lon if lon <= 180 else lon - 360


def url(kind: str, lat: float, lon: float, start: str, end: str) -> str:
    path = "/api/mhw" if kind == "json" else "/api/mhw/csv"
    return f"{path}?lon0={query_lon(lon)}&lat0={lat}&start={start}&end={end}&append={APPEND}"


def _num(v):
    if v is None or v == "" or (isinstance(v, str) and v.lower() in ("nan", "null", "none")):
        return None
    return float(v)


def row_month(r: dict) -> str:
    return str(r.get("date", ""))[:10]


def compare(rows: list[dict], want: dict[str, dict], default_tol: float = 1e-9) -> dict:
    """rows: API rows (JSON objects or CSV dicts); want: month "YYYY-MM-DD" ->
    {field: value or None}.

    The rows' months, in response order, must be exactly the wanted months in
    ascending order (Codex r22: no month outside the range, no duplicate, none
    missing, row count = month count; the API selects with a date slice, so its
    rows are chronological).  Each month's values must be within TOL
    (default_tol for the other fields: 1e-9 for JSON, 1e-5 for CSV text); a
    NaN / null is matched only by NaN / null."""
    want_months = sorted(want)
    got_months = [row_month(r) for r in rows]
    out = {"months": {}, "default_tol": default_tol, "n_rows": len(rows), "n_want": len(want_months),
           "dates_got": got_months, "dates_want": want_months, "date_list_match": got_months == want_months}
    out["match"] = out["date_list_match"]
    for month, w in want.items():
        hit = [r for r in rows if row_month(r) == month]
        m = {"n_rows": len(hit)}
        if len(hit) == 1:
            got = {k: _num(hit[0].get(k)) for k in FIELDS}
            m["got"], m["want"] = got, w
            m["match"] = all((got[k] is None and w[k] is None)
                             or (got[k] is not None and w[k] is not None and abs(got[k] - w[k]) <= TOL.get(k, default_tol))
                             for k in FIELDS)
        else:
            m["match"] = False
        out["months"][month] = m
        out["match"] = out["match"] and m["match"]
    return out
