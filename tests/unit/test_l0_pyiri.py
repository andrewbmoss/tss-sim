"""The l0/iri2020 tier: IRI-2020 through PyIRI.

Skipped where PyIRI is not installed.  The tests hold the binding, not PyIRI's
physics: that the tier returns what PyIRI returns, on the right axes and in SI,
that it refuses what it cannot do, and that the model's own peak parameters
reach the state.  How close PyIRI comes to Fortran IRI-2020 is recorded in
docs/pyiri-validation.md.
"""

from __future__ import annotations

import datetime as dt

import numpy as np
import pytest

pytest.importorskip("PyIRI")

import PyIRI  # noqa: E402
import PyIRI.sh_library as sh  # noqa: E402

from tss.config import load_campaign  # noqa: E402
from tss.constants import plasma_frequency  # noqa: E402
from tss.registry import build  # noqa: E402

EPOCH = "2024-03-01T12:00:00Z"
ALT = np.arange(1.0e5, 2.0e6 + 1.0, 1.0e4)


def _tier(**kw):
    return build("l0", "iri2020", **{"f107": 150.0, **kw})


def test_the_tier_needs_a_stated_f107():
    with pytest.raises(ValueError, match="F10.7"):
        build("l0", "iri2020")


@pytest.mark.parametrize("bad", [{"topside": "IRI-2001"}, {"foF2_coeff": "X"},
                                 {"hmF2_model": "X"}, {"f107": 20.0}])
def test_what_pyiri_cannot_do_is_refused(bad):
    with pytest.raises(ValueError):
        build("l0", "iri2020", **{"f107": 120.0, **bad})


def test_a_profile_is_a_plausible_ionosphere_in_si():
    state = _tier().state_at(42.6, -71.5, EPOCH, ALT)
    ne = np.asarray(state.ne)
    assert np.all(np.isfinite(ne)) and np.all(ne >= 0.0)
    peak = int(np.argmax(ne))
    assert 2.0e5 < ALT[peak] < 4.5e5
    assert 1.0e10 < ne[peak] < 5.0e12
    # The topside falls monotonically above the peak.
    assert np.all(np.diff(ne[peak:]) <= 0.0)
    # The model's own peak, not the grid's, reaches the state.
    extra = state.extra
    assert abs(extra["hmf2_km"] * 1e3 - ALT[peak]) <= 1.0e4
    assert extra["nmf2"] >= 0.99 * ne[peak]
    assert extra["fof2_hz"] == pytest.approx(float(plasma_frequency(extra["nmf2"])), rel=2e-3)
    assert "PyIRI" in state.source and "F10.7=150" in state.source


def test_the_tier_returns_what_pyiri_returns():
    tier = _tier()
    ne = tier.profiles([42.6], [-71.5], [dt.datetime(2024, 3, 1, 12)], ALT)[0]
    *_, edp = sh.IRI_density_1day(
        2024, 3, 1, np.array([12.0]), np.array([-71.5]), np.array([42.6]), ALT / 1e3,
        150.0, coeff_dir=PyIRI.coeff_dir, foF2_coeff="URSI", hmF2_model="SHU2015",
        coord="GEO", old_output=False,
    )
    assert np.allclose(ne, np.asarray(edp)[0][:, 0], rtol=1e-12, atol=0.0)


def test_samples_at_different_times_match_single_calls():
    tier = _tier()
    when = [dt.datetime(2024, 3, 1, 3), dt.datetime(2024, 7, 9, 15, 30)]
    both = tier.profiles([10.0, -35.0], [100.0, 20.0], when, ALT)
    one = tier.profiles([10.0], [100.0], [when[0]], ALT)[0]
    two = tier.profiles([-35.0], [20.0], [when[1]], ALT)[0]
    assert np.allclose(both[0], one) and np.allclose(both[1], two)


def test_longitude_convention_does_not_matter():
    tier = _tier()
    when = [dt.datetime(2024, 3, 1, 12)]
    assert np.allclose(tier.profiles([42.6], [-71.5], when, ALT),
                       tier.profiles([42.6], [288.5], when, ALT))


def test_solar_activity_raises_the_peak():
    low = _tier(f107=70.0).state_at(0.0, 0.0, EPOCH, ALT).extra["nmf2"]
    high = _tier(f107=200.0).state_at(0.0, 0.0, EPOCH, ALT).extra["nmf2"]
    assert high > 1.5 * low


def test_the_baseline_campaign_states_its_solar_level():
    l0 = load_campaign("baseline").config["tiers"]["l0"]
    assert l0["name"] == "iri2020" and l0["f107"] == 120.0
