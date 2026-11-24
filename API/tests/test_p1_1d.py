"""P1-1d targeted checks (spec 17 P1-1d): plots read one variable, no pyplot leak.

* disk level: series / month plots read only their sstvar's chunks;
* the plot input series (dates and order, NaN mask, values, dtype) equal the
  pre-P1-1d full-Dataset + dropna path (frozen copy in tests/reference);
* land (all NaN) still gives 200 with an empty line, Baltic (P04) keeps its
  interior NaN months so the line breaks instead of joining across them;
* current-year red line and palette branches are unchanged;
* month plots no longer register pyplot figure managers, on success or error,
  and leave no Figure alive.

No test reads a whole variable through .data / .values.
"""
import asyncio
import gc
from datetime import date

import matplotlib.axes
import matplotlib.figure
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pytest
from polars.testing import assert_frame_equal

from tests.reference import mhw_utils_fbb2926 as ref
from tests.test_p1_1a import COORDS, _vars_read, chunk_reads

# --------------------------------------------------------------- disk reads


@pytest.mark.parametrize("url,expected", [
    ("/api/mhw/plot?bbox=(135,15,137,17)&period=1982-2011,2012-2021,2022&sstvar=sst&mode=month", {"sst"}),
    ("/api/mhw/plot?bbox=(135,15,137,17)&period=2025,2026&sstvar=sst_anomaly&mode=month", {"sst_anomaly"}),
    ("/api/mhw/plot?bbox=(-90,-10,-80,0),(-150,-5,-90,5)&period=2015-202306&sstvar=sst_anomaly", {"sst_anomaly"}),
    ("/api/mhw/plot?bbox=(-5,-5,5,5)&period=2020-202112&sstvar=sst", {"sst"}),          # crosses 0
    ("/api/mhw/plot?bbox=(170,-5,-170,5)&period=2020-202112&sstvar=sst_anomaly", {"sst_anomaly"}),
    ("/api/mhw/plot?bbox=(-5,-5,5,5)&period=1982-2011&sstvar=sst&mode=month", {"sst"}),
])
def test_plot_reads_only_its_variable(client, url, expected):
    with chunk_reads() as keys:
        r = client.get(url)
    assert r.status_code == 200
    assert r.headers["content-type"] == "image/png"
    assert _vars_read(keys) - COORDS == expected


# ------------------------------------------- plot input series vs old path

D = date.fromisoformat
FULL_END = "2026-10-01"  # month_climatology's end under the fixed test clock
SERIES = [
    # (lon0, lat0, lon1, lat1, start, end, sstvar, mode)
    (135, 15, 137, 17, "1982-01-01", FULL_END, "sst", "month_mean"),               # P01
    (100, 40, None, None, "1982-01-01", FULL_END, "sst", "month_mean"),            # P02 land point
    (135, 15, 137, 17, "1982-01-01", FULL_END, "sst_anomaly", "month_mean"),       # P06 variable
    (-5, -5, 5, 5, "1982-01-01", FULL_END, "sst", "month_mean"),                   # month across 0
    (-90, -10, -80, 0, "2015-01-01", "2023-06-01", "sst_anomaly", "area_mean"),    # P03 bbox 1
    (-150, -5, -90, 5, "2015-01-01", "2023-06-01", "sst_anomaly", "area_mean"),    # P03 bbox 2
    (20, 55, 21, 56, "2015-01-01", "2023-06-01", "sst_anomaly", "area_mean"),      # P04 Baltic
    (135, 15, None, None, "1982-01-01", FULL_END, "sst", "area_mean"),             # P05 point
    (140, 20, None, None, "1982-01-01", FULL_END, "sst", "area_mean"),             # P05 second point
    (-5, -5, 5, 5, "2020-01-01", "2021-12-01", "sst", "area_mean"),                # crosses 0
    (170, -5, -170, 5, "2020-01-01", "2021-12-01", "sst_anomaly", "area_mean"),    # antimeridian (J08 behaviour)
    (100, 40, 101, 41, "2020-01-01", "2021-12-01", "sst", "area_mean"),            # land bbox
]


def _missing(s: pl.Series) -> list[bool]:
    """Missing-value mask: output_df's pl.from_pandas turns NaN into null."""
    return (s.is_null() | s.fill_null(0.0).is_nan()).to_list()


def _series(fn, args):
    lon0, lat0, lon1, lat1, start, end, var, mode = args
    m = mode if mode == "month_mean" else f"area_mean_{var}"
    _, df = asyncio.run(fn(lon0=lon0, lat0=lat0, lon1=lon1, lat1=lat1, start=D(start), end=D(end),
                           append=var, mode=m))
    return df


@pytest.fixture(scope="module")
def series_pairs(client):
    from src.mhw_utils import process_mhw_data

    return {a: (_series(process_mhw_data, a), _series(ref.process_mhw_data, a)) for a in SERIES}


@pytest.mark.parametrize("args", SERIES, ids=[f"{a[7]}_{a[6]}_{a[0]}_{a[1]}_{a[2]}_{a[3]}" for a in SERIES])
def test_plot_series_equal_old_full_dataset_path(series_pairs, args):
    new, old = series_pairs[args]
    var = args[6]
    assert new.columns == old.columns and new.schema == old.schema          # dtype
    assert new["date"].to_list() == old["date"].to_list()                   # dates and order
    assert _missing(new[var]) == _missing(old[var])                         # NaN mask
    assert_frame_equal(new, old, check_exact=True)                          # values


def test_land_is_all_nan_and_keeps_every_date(series_pairs):
    new, _ = series_pairs[SERIES[1]]
    assert new.height == 538 and all(_missing(new["sst"]))


def test_baltic_keeps_interior_nan_months(series_pairs):
    new, _ = series_pairs[SERIES[6]]
    v = _missing(new.sort("date")["sst_anomaly"])
    assert new.height == 102  # 2015-01 .. 2023-06, every month present
    interior = [i for i in range(1, len(v) - 1) if v[i] and not v[i - 1] or v[i] and not v[i + 1]]
    assert interior, "expected NaN months between valid months (sea ice)"
    # what dropna on the single variable would have done (the reason it was removed)
    assert new.drop_nulls().height < new.height


# --------------------------------------------------- what reaches ax.plot


@pytest.fixture()
def plot_calls(monkeypatch):
    calls = []
    real = matplotlib.axes.Axes.plot

    def spy(self, *args, **kw):
        calls.append({"x": list(args[0]), "y": np.asarray(args[1], dtype=float),
                      "color": kw.get("color"), "label": kw.get("label")})
        return real(self, *args, **kw)

    monkeypatch.setattr(matplotlib.axes.Axes, "plot", spy)
    return calls


def test_baltic_line_breaks_at_nan(client, plot_calls):
    r = client.get("/api/mhw/plot?bbox=(20,55,21,56)&period=2015-202306&sstvar=sst_anomaly")
    assert r.status_code == 200 and len(plot_calls) == 1
    y = plot_calls[0]["y"]
    assert len(y) == 102 and np.isnan(y).any() and not np.isnan(y).all()
    # NaN points sit between valid ones, so matplotlib draws separate segments
    nan_idx = np.flatnonzero(np.isnan(y))
    assert ((nan_idx > 0) & (nan_idx < len(y) - 1)).any()


def test_land_month_plot_is_200_with_empty_lines(client, plot_calls):
    r = client.get("/api/mhw/plot?bbox=(100,40)&period=1982-2011,2012-2021&sstvar=sst&mode=month")
    assert r.status_code == 200 and r.headers["content-type"] == "image/png"
    assert len(plot_calls) == 2 and all(np.isnan(c["y"]).all() for c in plot_calls)


def test_current_year_is_red_and_palette_branches(client, plot_calls):
    # default palette; 2026 is the fixed clock's current year -> red
    assert client.get("/api/mhw/plot?bbox=(135,15,137,17)&period=2025,2026&sstvar=sst_anomaly&mode=month").status_code == 200
    by_label = {c["label"]: c["color"] for c in plot_calls}
    assert by_label["2026"] == "red" and by_label["2025"] != "red"
    plot_calls.clear()
    # user palette shorter than periods: given colour first, cubehelix fills the rest; red still wins
    assert client.get("/api/mhw/plot?bbox=(135,15,137,17)&period=1982-2011,2012-2021,2026"
                      "&sstvar=sst&mode=month&palette=blue").status_code == 200
    by_label = {c["label"]: c["color"] for c in plot_calls}
    assert by_label["1982-2011"] == "blue" and by_label["2012-2021"] not in ("blue", "red")
    assert by_label["2026"] == "red"
    assert [len(c["y"]) for c in plot_calls] == [12, 12, 10]  # 2026 has Jan..Oct
    plot_calls.clear()
    # series palette (P05)
    assert client.get("/api/mhw/plot?bbox=(135,15),(140,20)&period=all&palette=blue,green").status_code == 200
    assert [c["color"] for c in plot_calls] == ["blue", "green"]


# ---------------------------------------------------- figure lifetime / leak


def _live_figures():
    gc.collect()
    return sum(isinstance(o, matplotlib.figure.Figure) for o in gc.get_objects())


def test_month_plot_registers_no_pyplot_figure(client):
    before_nums, before_live = plt.get_fignums(), _live_figures()
    for _ in range(3):
        r = client.get("/api/mhw/plot?bbox=(135,15,137,17)&period=1982-2011,2012-2021,2022&sstvar=sst&mode=month")
        assert r.status_code == 200
    assert plt.get_fignums() == before_nums
    assert _live_figures() == before_live  # nothing keeps the Figure after the PNG is written


def test_month_plot_error_registers_no_pyplot_figure(client):
    from src.mhw_plot import month_climatology

    before_nums, before_live = plt.get_fignums(), _live_figures()
    with pytest.raises(ValueError):  # "abc" fails after the figure was created
        asyncio.run(month_climatology(lon0=135, lat0=15, lon1=137, lat1=17, period="abc", sstvar="sst"))
    assert plt.get_fignums() == before_nums
    assert _live_figures() == before_live


def test_series_plot_registers_no_pyplot_figure(client):
    before_nums, before_live = plt.get_fignums(), _live_figures()
    assert client.get("/api/mhw/plot?bbox=(135,15,137,17)&period=2019-202006&sstvar=sst").status_code == 200
    assert plt.get_fignums() == before_nums
    assert _live_figures() == before_live


def test_mhw_plot_does_not_use_pyplot():
    import ast
    from pathlib import Path

    src = (Path(__file__).resolve().parents[1] / "src/mhw_plot.py").read_text()
    tree = ast.parse(src)
    mods = {a.name for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
    mods |= {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    assert "matplotlib.pyplot" not in mods
    assert sum(1 for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names if a.name == "numpy") == 1
