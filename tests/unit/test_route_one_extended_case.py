"""WP5 Phase 1, gate G1 decision 4 (2026-11-24): the 0.1 s fading case of the
route-one Monte Carlo is an extended reflector and is no longer run there."""

from __future__ import annotations

import numpy as np
import pytest

from tss.l2_instrument.fading import halo_extent_m
from tss.studies import route_one as r1


def test_compact_facet_cases_exclude_0p1_s():
    assert [tau for _, tau in r1.FADING_CASES] == [None, 0.5]
    assert r1.EXTENDED_FADING_CASES == (("tau_0p1", 0.1),)
    # The kept cases keep their 2026-10-06 codes, so their seeds are unchanged.
    assert r1.FADING_CASE_CODE["none"] == 0 and r1.FADING_CASE_CODE["tau_0p5"] == 2


def test_the_moved_case_is_wider_than_a_compact_facet():
    # D_length scenes switch to extended between 2 and 4 km; Fresnel 3.9 km.
    assert halo_extent_m(0.1, 5.0e6, 5.0e5) > 4.0e3
    assert halo_extent_m(0.5, 5.0e6, 5.0e5) < max(r1.FACETS_KM) * 1e3


def test_mc_table_has_only_the_compact_cases(tmp_path):
    class _Manifest:
        outputs: list = []

    result = r1.run(tmp_path, config={"route_one": {"n_trials": 8, "seed": 1}}, manifest=_Manifest())
    mc = np.loadtxt(tmp_path / "route_one_overflight_mc.csv", delimiter=",", skiprows=1)
    assert set(mc[:, 1].tolist()) == {0.0, 2.0}
    assert mc.shape == (2 * len(r1.FACETS_KM), 6)
    assert result.metrics["halo_km_tau_0p1_5mhz_500km"] == pytest.approx(13.8, abs=0.1)
    assert result.metrics["extended_fading_cases"] == "tau_0p1"
