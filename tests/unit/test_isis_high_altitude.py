"""Unit tests for the WP0 Phase 5 ISIS-1 high-altitude comparison (scripts/isis_high_altitude.py)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def iha():
    spec = importlib.util.spec_from_file_location(
        "isis_high_altitude", REPO / "scripts" / "isis_high_altitude.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_bands_cover_900_to_1400(iha):
    assert iha.BANDS[0] == (900.0, 1000.0) and iha.BANDS[-1] == (1300.0, 1400.0)
    assert len(iha.BANDS) == 5


def test_band_rows_edges_and_realised_coverage(iha):
    rng = np.random.default_rng(3)
    n = 2000
    hs = rng.uniform(900.5, 1400.0, n)
    rz = np.where(rng.uniform(size=n) < 0.5, 120.0, 10.0)
    fn_obs = rng.lognormal(np.log(1.0), 0.3, n)
    # NeQuick 1 dB high everywhere, IRI-2001 exact, IRI-2001 corrected 10 % high:
    # the ensemble edge is IRI-2001's, so it realises 95 % coverage.
    models = {"epoch_NeQuick": fn_obs * 10 ** (1 / 20), "epoch_IRI2001": fn_obs.copy(),
              "epoch_IRI2001corr": fn_obs * 1.1}
    rows = iha.band_rows(fn_obs, models, hs, rz, np.full(n, 0.9))
    allrow = rows["all"]["900-1400"]
    assert allrow["n"] == n
    assert allrow["median_db"]["epoch_NeQuick"] == pytest.approx(1.0)
    assert allrow["ensemble_p05_mhz"] == pytest.approx(allrow["p05_mhz"]["observed"])
    assert allrow["realised_coverage_of_ensemble_edge"] == pytest.approx(0.95, abs=0.002)
    assert allrow["h0_fit_over_iri_law"][1] == pytest.approx(0.9)
    # Bands are half-open (lo, hi]; the five bands add up to the whole range.
    assert sum(rows["all"][f"{int(a)}-{int(b)}"]["n"] for a, b in iha.BANDS) == n
    assert set(rows) == {"solar_max", "intermediate", "solar_min", "all"}
    assert "900-1000" not in rows["intermediate"]          # no profiles there


def test_overstated_tail_lowers_realised_coverage(iha):
    rng = np.random.default_rng(4)
    n = 4000
    hs = np.full(n, 1200.0)
    rz = np.full(n, 150.0)
    fn_obs = rng.lognormal(0.0, 0.3, n)
    up = {k: fn_obs * 1.2 for k in ("epoch_NeQuick", "epoch_IRI2001", "epoch_IRI2001corr")}
    row = iha.band_rows(fn_obs, up, hs, rz, np.full(n, np.nan))["solar_max"]["1100-1200"]
    assert row["ensemble_over_observed_p05"] == pytest.approx(1.2, rel=1e-6)
    assert row["realised_coverage_of_ensemble_edge"] < 0.9
    assert "h0_fit_over_iri_law" not in row


def test_splits_separate_polar_from_low_latitude(iha):
    rng = np.random.default_rng(5)
    n = 3000
    lat = np.where(np.arange(n) < n // 2, 10.0, 75.0)
    fn_obs = rng.lognormal(0.0, 0.3, n)
    model = fn_obs.copy()
    model[lat >= 60] *= 2.0                         # the model twice too dense in the cap
    models = {k: model for k in ("epoch_NeQuick", "epoch_IRI2001", "epoch_IRI2001corr")}
    sp = iha.splits(fn_obs, models, np.full(n, 1200.0), np.full(n, 150.0), lat,
                    np.full(n, 12.0), np.full(n, 6.0))
    assert sp["solar_max_abs_lat_0_30"]["ensemble_over_observed_p05"] == pytest.approx(1.0)
    assert sp["solar_max_abs_lat_70_80"]["ensemble_over_observed_p05"] == pytest.approx(2.0)
    assert sp["solar_max_polar_1100_1200"]["realised_coverage_of_ensemble_edge"] < 0.5
    assert sp["n_solar_min_polar"] == 0
