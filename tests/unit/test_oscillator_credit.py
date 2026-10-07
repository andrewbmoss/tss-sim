"""The range-correlation credit is an integral over the noise spectrum.

De-chirping at a round-trip delay takes a lag difference of the oscillator
phase, and the coherent combination sees the change of that difference across
the coherent interval.  These tests hold the closed forms to a numerical
integral of the two filters, and record how far the single-offset shortcuts are
from it, since Requirements Rev G carried one of them.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import integrate

from tss.l2_instrument.oscillator import (
    FM_NOISE_TYPES,
    INS26_ALLAN_GOAL,
    INS26_GOAL_INTERVAL_S,
    Oscillator,
    range_correlation_credit,
    range_correlation_credit_db,
)

T_GOAL = 1.1539  # INS-15 goal interval, unrounded
TAU_MAX = 7.0e-3  # round trip at INS-04's 1050 km


def _numeric_credit(t, tau, exponent):
    """Ratio of uncredited to surviving phase jitter by direct integration."""
    grid = np.logspace(-6, 4.5, 4000)

    def surviving(f):
        return f**exponent * (tau * np.sinc(f * tau)) ** 2 * 4.0 * np.sin(np.pi * f * t) ** 2

    def allan(f):
        return f**exponent * 2.0 * np.sin(np.pi * f * t) ** 4 / (np.pi * f * t) ** 2

    s = integrate.quad(surviving, 1e-7, 5e4, limit=200000, points=grid)[0]
    a = integrate.quad(allan, 1e-7, 5e4, limit=200000, points=grid)[0]
    return float(np.sqrt(t**2 * a / s))


@pytest.mark.filterwarnings("ignore")
@pytest.mark.parametrize("noise", sorted(FM_NOISE_TYPES))
@pytest.mark.parametrize(("t", "tau"), [(T_GOAL, TAU_MAX), (0.5, 4.0e-3), (0.05, 7.0e-3)])
def test_closed_form_matches_the_integral(noise, t, tau):
    expected = _numeric_credit(t, tau, FM_NOISE_TYPES[noise])
    assert range_correlation_credit(t, tau, noise) == pytest.approx(expected, rel=2e-3)


def test_credit_at_the_ins26_goal_by_noise_type():
    assert range_correlation_credit(T_GOAL, TAU_MAX, "white_fm") == pytest.approx(9.08, abs=0.02)
    assert range_correlation_credit(T_GOAL, TAU_MAX, "flicker_fm") == pytest.approx(53.4, abs=0.2)
    assert range_correlation_credit(T_GOAL, TAU_MAX, "random_walk_fm") == pytest.approx(95.3, abs=0.2)


def test_single_offset_shortcuts_against_the_integral():
    flicker_db = 20.0 * np.log10(range_correlation_credit(T_GOAL, TAU_MAX, "flicker_fm"))
    half = float(range_correlation_credit_db(1.0 / (2.0 * T_GOAL), TAU_MAX))
    full = float(range_correlation_credit_db(1.0 / T_GOAL, TAU_MAX))
    # 1/(2T) is where the difference across T is most sensitive, and it
    # reproduces the flicker case.  1/T is its null and gives 6 dB less.
    assert half == pytest.approx(flicker_db, abs=0.2)
    assert full == pytest.approx(half - 6.02, abs=0.01)
    assert full == pytest.approx(28.4, abs=0.05)
    white_db = 20.0 * np.log10(range_correlation_credit(T_GOAL, TAU_MAX, "white_fm"))
    assert half - white_db > 15.0


def test_delay_equal_to_the_interval_gives_the_allan_convention():
    """With no separation the factor is 1/sqrt(2) for every type.

    The surviving variance is then that of the difference of two adjacent
    averages, which is twice the Allan variance by its definition.
    """
    for noise in FM_NOISE_TYPES:
        assert range_correlation_credit(1.0, 1.0, noise) == pytest.approx(np.sqrt(0.5))


def test_bad_arguments_are_refused():
    with pytest.raises(ValueError):
        range_correlation_credit(1.0, 7e-3, "pink")
    with pytest.raises(ValueError):
        range_correlation_credit(1e-3, 7e-3)
    with pytest.raises(ValueError):
        Oscillator(noise_type="pink")


def test_oscillator_applies_the_credit_of_its_noise_type():
    raw = Oscillator().phase_jitter_rad(3.0e6, T_GOAL)
    for noise in FM_NOISE_TYPES:
        osc = Oscillator(noise_type=noise)
        credited = osc.phase_jitter_credited_rad(3.0e6, T_GOAL, TAU_MAX)
        assert float(raw / credited) == pytest.approx(range_correlation_credit(T_GOAL, TAU_MAX, noise))
    # Not single-reference: nothing to correlate with, so no credit.
    lone = Oscillator(single_reference=False)
    assert float(lone.phase_jitter_credited_rad(3.0e6, T_GOAL, TAU_MAX)) == pytest.approx(float(raw))


def test_score_reports_the_credit_and_scores_without_it():
    score = Oscillator().score_ins26(
        3.0e6, INS26_GOAL_INTERVAL_S, TAU_MAX, threshold=INS26_ALLAN_GOAL
    )
    assert score["scored_on"] == "uncredited"
    assert score["verdict"] == "pass"  # the 2e-9 part against the 3.2e-9 goal
    assert score["noise_type"] == "flicker_fm"
    assert score["credit_factor"] == pytest.approx(53.2, abs=0.5)
    assert set(score["credit_factor_by_noise_type"]) == set(FM_NOISE_TYPES)
    assert Oscillator.from_config({"noise_type": "white_fm"}).noise_type == "white_fm"
