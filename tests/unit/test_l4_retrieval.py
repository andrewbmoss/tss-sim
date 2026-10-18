"""Unit tests for L4, retrieval.

The lamination inversion is tested two ways and the two answer different
questions.

Against its own forward operator, the inversion must be exact: a profile that
really is piecewise linear in plasma frequency against height is recovered to
machine precision, and any error there is arithmetic.

Against a closed-form exponential topside, it must be accurate rather than
exact, because such a profile is not piecewise linear and what is being
measured is the lamina approximation itself.  That case has a known answer:
``R'(f) = 2 H arccosh(f / f_0)`` and ``h(f) = h_s - 2 H ln(f / f_0)``.  The
stated tolerance is 500 m of height error with 80 laminae, which is one
twentieth of an INS-05 threshold range cell, and the error is checked to fall
as the grid refines so that a coincidental pass is distinguishable from a
correct one.
"""

from __future__ import annotations

import numpy as np
import pytest

from tss.constants import density_from_plasma_frequency
from tss.interfaces import DensityProfile, Ionogram, PlatformState
from tss.l4_retrieval import (
    LaminationRetrieval,
    depletion,
    error_budget,
    lamina_thicknesses_m,
    lamination_matrix,
    trace_points,
    virtual_range_from_laminae_m,
    with_uncertainty,
)
from tss.registry import available

ALT_S = 5.0e5          # spacecraft altitude, m
F0 = 5.0e5             # plasma frequency at the spacecraft, Hz
SCALE_H = 8.0e4        # topside scale height of the analytic case, m


def _exponential_case(n: int = 80) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Frequencies, group ranges and true heights for the analytic topside.

    ``f_N(h) = F0 exp((ALT_S - h) / 2H)``, for which the group range integral
    has the closed form ``2 H arccosh(f / F0)``.
    """
    freq = F0 * np.exp(np.linspace(0.02, 2.0, n))
    virtual = 2.0 * SCALE_H * np.arccosh(freq / F0)
    heights = ALT_S - 2.0 * SCALE_H * np.log(freq / F0)
    return freq, virtual, heights


def _ionogram(freq: np.ndarray, virtual: np.ndarray, **flags) -> Ionogram:
    base = {
        "range_resolution_m": 1.8e3,
        "unambiguous_range_m": 2.6e6,
        "ionogram_period_s": 3.3,
        "mode_discrimination": False,
    }
    base.update(flags)
    return Ionogram(
        freq_hz=freq,
        virtual_range_m=virtual[:, None],
        snr_db=np.full((freq.size, 1), 25.0),
        modes=("unresolved",),
        flags=base,
    )


def _platform() -> PlatformState:
    return PlatformState(
        alt_m=ALT_S,
        lat_deg=0.0,
        lon_deg=0.0,
        velocity_ned_mps=np.array([7600.0, 0.0, 0.0]),
        epoch="2026-01-01T00:00:00Z",
    )


# ---------------------------------------------------------------------------
# The lamination algebra
# ---------------------------------------------------------------------------


def test_lamination_matrix_is_lower_triangular():
    """A lamina deeper than the reflection level contributes nothing.

    The wave never reaches it, so the coefficient is zero.  If the matrix were
    not triangular the forward substitution would be wrong and the inversion
    would need a full solve.
    """
    freq = np.array([1.0, 2.0, 4.0, 8.0]) * 1.0e6
    g = lamination_matrix(freq, 0.5e6)
    assert np.allclose(np.triu(g, 1), 0.0)
    assert np.all(np.diag(g) > 0.0)


def test_lamination_matrix_refuses_an_evanescent_frequency():
    """A sounding below the local plasma frequency cannot have echoed.

    Accepting it would take the arcsine of a number above one.  Raising is
    right here: this is a programming or configuration error, not a physical
    state the tier should report as a masked result.
    """
    with pytest.raises(ValueError, match="evanescent"):
        lamination_matrix(np.array([4.0e5, 1.0e6]), 5.0e5)


def test_lamination_inverts_its_own_forward_operator_exactly():
    """``R' = G t`` and the forward substitution are exact inverses."""
    freq = np.array([1.0, 1.4, 2.0, 3.0, 4.5, 6.0, 9.0]) * 1.0e6
    thickness = np.array([2.0e4, 3.0e4, 5.0e4, 4.0e4, 6.0e4, 2.5e4, 3.5e4])
    virtual = virtual_range_from_laminae_m(thickness, freq, 0.8e6)
    recovered = lamina_thicknesses_m(virtual, freq, 0.8e6)
    assert np.allclose(recovered, thickness, rtol=0.0, atol=1.0e-6)


def test_lamination_recovers_the_analytic_exponential_topside():
    """Heights within 500 m of the closed-form answer with 80 laminae.

    500 m is one twentieth of an INS-05 threshold range cell, so the lamina
    approximation is well inside the measurement error it will be used with.
    """
    freq, virtual, truth = _exponential_case(80)
    heights = ALT_S - np.cumsum(lamina_thicknesses_m(virtual, freq, F0))
    assert np.max(np.abs(heights - truth)) < 500.0


def test_lamination_error_falls_as_the_grid_refines():
    """The approximation converges, which distinguishes correct from lucky.

    A single tolerance can be met by accident.  An error that falls by more
    than a factor of four as the lamina count is doubled twice is the lamina
    approximation behaving as an approximation should.
    """
    errors = []
    for n in (20, 40, 80):
        freq, virtual, truth = _exponential_case(n)
        heights = ALT_S - np.cumsum(lamina_thicknesses_m(virtual, freq, F0))
        errors.append(float(np.max(np.abs(heights - truth))))
    assert errors[0] > errors[1] > errors[2]
    assert errors[0] / errors[2] > 4.0


# ---------------------------------------------------------------------------
# The registered tier
# ---------------------------------------------------------------------------


def test_lamination_tier_is_registered_and_inverts():
    assert "lamination" in available("l4")["l4"]
    freq, virtual, truth = _exponential_case(80)
    profile = LaminationRetrieval(start_plasma_hz=F0).invert(
        _ionogram(freq, virtual), _platform()
    )
    assert isinstance(profile, DensityProfile)
    assert profile.alt_m.size == freq.size + 1          # plus the spacecraft point
    assert np.all(np.diff(profile.alt_m) > 0.0)         # ascending altitude
    recovered = np.sort(profile.alt_m)[::-1][1:]
    assert np.max(np.abs(recovered - truth)) < 500.0


def test_lamination_tier_records_which_starting_value_it_used():
    """The three starting-value routes give different error budgets.

    ``docs/architecture.md`` requires the retrieval to say which one it used,
    because a profile that does not say is not usable as evidence.
    """
    freq, virtual, _ = _exponential_case(40)
    given = LaminationRetrieval(start_plasma_hz=F0).invert(
        _ionogram(freq, virtual), _platform()
    )
    assert given.flags["starting_value"] == "in_situ"

    fallback = LaminationRetrieval().invert(_ionogram(freq, virtual), _platform())
    assert fallback.flags["starting_value"] == "lowest_trace"
    assert fallback.flags["dropped_below_start"] == 1


def test_lamination_tier_flags_a_non_monotonic_result_rather_than_clamping():
    """A non-positive lamina is the only symptom the algebra gives.

    A depletion breaks the monotonicity the method assumes.  Clamping the
    thickness to zero would remove the one piece of evidence that the
    assumption failed, so the value is returned unclamped and flagged.
    """
    freq = np.array([1.0, 1.5, 2.0, 3.0]) * 1.0e6
    virtual = np.array([5.0e4, 1.2e5, 6.0e4, 2.0e5])   # group range falls: not monotonic
    profile = LaminationRetrieval(start_plasma_hz=0.8e6).invert(
        _ionogram(freq, virtual), _platform()
    )
    assert profile.flags["n_nonpositive_laminae"] >= 1
    assert profile.flags["heights_monotonic"] is False
    assert "monotonicity_warning" in profile.flags


def test_lamination_tier_returns_an_empty_profile_for_an_empty_trace():
    """No detection is a result, not an exception.  Contract rule 3."""
    freq = np.array([1.0, 2.0]) * 1.0e6
    empty = Ionogram(
        freq_hz=freq,
        virtual_range_m=np.full((2, 1), np.nan),
        snr_db=np.full((2, 1), np.nan),
        modes=("unresolved",),
    )
    profile = LaminationRetrieval().invert(empty, _platform())
    assert profile.alt_m.size == 0
    assert profile.flags["no_trace"] is True


def test_lamination_tier_states_that_the_field_was_ignored():
    """The no-field group index is the largest known bias and it is flagged."""
    freq, virtual, _ = _exponential_case(20)
    profile = LaminationRetrieval(start_plasma_hz=F0).invert(
        _ionogram(freq, virtual), _platform()
    )
    assert profile.flags["magnetoionic"] is False
    assert profile.flags["vertical_incidence_assumed"] is True


def test_lamination_tier_does_not_mutate_its_input():
    """Contract rule 1."""
    freq, virtual, _ = _exponential_case(20)
    ionogram = _ionogram(freq, virtual)
    before = ionogram.virtual_range_m.copy()
    LaminationRetrieval(start_plasma_hz=F0).invert(ionogram, _platform())
    assert np.array_equal(ionogram.virtual_range_m, before)


def test_lamination_tier_describe_round_trips_through_json():
    import json

    record = LaminationRetrieval().describe()
    assert json.loads(json.dumps(dict(record)))["name"] == "lamination"


def test_trace_points_drops_undetected_frequencies_rather_than_interpolating():
    """An interpolated trace point would enter the inversion as a measurement."""
    freq = np.array([1.0, 2.0, 3.0, 4.0]) * 1.0e6
    virtual = np.array([1.0e5, np.nan, 2.0e5, 3.0e5])
    f, r, _ = trace_points(_ionogram(freq, virtual))
    assert f.size == 3
    assert np.all(np.isfinite(r))


def test_density_from_frequency_round_trips():
    """``Ne = (f / 8.98)^2`` against its inverse."""
    from tss.constants import plasma_frequency

    freq = np.array([5.0e5, 3.0e6, 1.0e7])
    assert np.allclose(plasma_frequency(density_from_plasma_frequency(freq)), freq)


# ---------------------------------------------------------------------------
# Stubs
# ---------------------------------------------------------------------------


def test_magnetoionic_lamination_is_no_longer_a_stub():
    """The field-corrected inversion runs.  Its own tests are in test_magnetoionic.py."""
    profile = LaminationRetrieval(
        magnetoionic=True, geomagnetic_model="dipole"
    ).invert(_ionogram(*_exponential_case(10)[:2]), _platform())
    assert profile.flags["magnetoionic"] is True
    assert profile.flags["field_model"]["model"] == "centred dipole"


@pytest.mark.parametrize(
    "call",
    [
        lambda: error_budget.range_covariance_m2(None),
        lambda: error_budget.propagate_covariance(np.zeros((2, 2)), np.ones(2), 1.0),
        lambda: error_budget.bias_budget_m(None, None),
        lambda: with_uncertainty(None, None, None),
    ],
)
def test_error_budget_is_a_stub_naming_wp6(call):
    """Until this lands, ``ne_sigma`` stays ``None`` and no accuracy is quoted."""
    with pytest.raises(NotImplementedError, match="WP6"):
        call()


@pytest.mark.parametrize(
    "call",
    [
        lambda: depletion.detect_depletion_crossing(None),
        lambda: depletion.segment_trace(None, np.zeros(1)),
        lambda: depletion.invert_across_depletion(None, None),
    ],
)
def test_depletion_inversion_is_a_stub_naming_wp10(call):
    with pytest.raises(NotImplementedError, match="WP10"):
        call()
