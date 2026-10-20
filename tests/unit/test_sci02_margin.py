"""Unit tests for the SCI-02 margin analysis (WP0 Phase 3)."""

from __future__ import annotations

import numpy as np
import pytest

from tss.studies import sci02_margin as m


def test_parse_sample_name():
    info = m.parse_sample_name("O2_500km_max2000h_nequick.npz")
    assert info == dict(orbit_class="O2", altitude_km=500.0, solar="max2000h", topside="nequick")
    assert m.parse_sample_name("O1_900km_min_iri-2001-corrected.npz")["topside"] == \
        "iri-2001-corrected"
    assert m.parse_sample_name("ensemble.csv") is None


def test_margin_fractions_and_tails():
    # 1000 samples: hmF2 from 100 to 599 km, spacecraft at 500 km.
    hm = np.arange(100.0, 600.0, 0.5) * 1e3
    st = m.margin_stats(500e3, hm)
    assert st["n"] == 1000
    assert st["frac_margin_lt_0km"] == pytest.approx(0.199)    # hmF2 > 500 km
    assert st["frac_margin_lt_100km"] == pytest.approx(0.399)  # hmF2 > 400 km
    assert st["min_margin_km"] == pytest.approx(-99.5)
    assert st["frac_hmf2_below_envelope"] == pytest.approx(0.2)
    assert st["frac_hmf2_above_envelope"] == pytest.approx(0.199)
    assert st["hmf2_max_km"] == pytest.approx(599.5)


def test_margin_is_strictly_below():
    # A peak exactly at the spacecraft is margin 0, which is not "below 0".
    st = m.margin_stats(500e3, np.array([500e3, 400e3]))
    assert st["frac_margin_lt_0km"] == 0.0
    assert st["frac_margin_lt_100km"] == 0.5


def test_nan_hmf2_is_excluded_and_counted():
    st = m.margin_stats(500e3, np.array([300e3, np.nan, 550e3]))
    assert st["n"] == 2 and st["n_invalid"] == 1
    assert st["frac_margin_lt_0km"] == 0.5


def test_bins_are_half_open():
    lt = np.array([0.0, 2.99, 3.0, 23.9])
    hm = np.full(4, 300e3)
    rows = m.binned(lt, (0, 3, 6, 24), 500e3, hm, "lt_h")
    assert [r["n"] for r in rows] == [2, 1, 1]


def test_verdict_names_the_worst_cell():
    cells = [
        dict(run="a", orbit_class="O1", altitude_km=500.0, solar="max",
             **m.margin_stats(500e3, np.array([300e3, 510e3]))),
        dict(run="a", orbit_class="O1", altitude_km=900.0, solar="max",
             **m.margin_stats(900e3, np.array([300e3, 450e3]))),
        dict(run="a", orbit_class="O1", altitude_km=1200.0, solar="max",
             **m.margin_stats(1200e3, np.array([150e3, 1300e3]))),
    ]
    v = m.verdict(cells)
    assert v["n_cells"] == 2                       # 1200 km is outside the threshold range
    assert not v["above_peak"]
    assert v["min_fraction_above_peak"] == pytest.approx(0.5)
    assert not v["above_peak_fraction_met"]
    assert v["worst_below_peak"]["cell"] == "a:O1/500 km/max"
    assert v["highest_hmf2_p999"]["cell"] == "a:O1/500 km/max"


def test_cell_tables_read_samples(tmp_path):
    d = tmp_path / "samples"
    d.mkdir()
    n = 400
    rng = np.random.default_rng(1)
    hm = rng.uniform(250e3, 520e3, n).astype(np.float32)
    cols = np.array(["nmf2", "hmf2_km", "kp"])
    extra = np.stack([np.ones(n), hm / 1e3, rng.uniform(0, 9, n)], axis=1).astype(np.float32)
    for top in ("nequick", "iri-2001"):
        np.savez_compressed(d / f"O2_500km_max_{top}.npz", hmf2_m=hm,
                            local_time_hours=rng.uniform(0, 24, n).astype(np.float32),
                            extra=extra, extra_columns=cols, h_sat_m=500e3)
    cells, lt, kp = m.cell_tables(d, "t")
    assert len(cells) == 1 and len(lt) == 8 and len(kp) == 4
    assert sum(r["n"] for r in lt) == n and sum(r["n"] for r in kp) == n
    assert cells[0]["frac_margin_lt_0km"] == pytest.approx(np.mean(hm > 500e3))
    assert m.topside_independence(d) == 0.0


def test_mis01_fraction_reading():
    """MIS-01 as restated 2026-10-20: above the peak in at least 99.9 % of samples."""
    assert m.MIS01_ABOVE_PEAK_FRACTION == 0.999
    hm = np.full(10000, 300e3)
    hm[:7] = 520e3                                  # 0.07 %, the Phase 3 worst cell
    cells = [dict(run="a", orbit_class="O2", altitude_km=500.0, solar="max",
                  **m.margin_stats(500e3, hm))]
    v = m.verdict(cells)
    assert not v["above_peak"]
    assert v["above_peak_fraction_met"]
    assert v["min_fraction_above_peak"] == pytest.approx(0.9993)
    hm[:11] = 520e3                                 # 0.11 % fails
    cells[0].update(m.margin_stats(500e3, hm))
    assert not m.verdict(cells)["above_peak_fraction_met"]
