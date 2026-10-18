"""The magnetic field in the forward model, the inversion and the scoring.

The ordinary wave's group index near reflection is larger than the field-free
one by the cosecant of the angle between the field and the vertical (Budden
1961, eq. 12.12).  A field-free inversion of a real topside trace is therefore
wrong by tens of kilometres at mid and high magnetic latitude.  These tests
hold the group index to its definition, the magnetoionic lamination to a known
profile, and the scoring to the rule that a profile is scored against SCI-05
and SCI-07 only when the field was treated the same way on both sides.
"""

from __future__ import annotations

import numpy as np
import pytest

from tss.constants import density_from_plasma_frequency
from tss.interfaces import DensityProfile, Ionogram, IonosphereState, PlatformState
from tss.l1_propagation import magnetoionic as mi
from tss.l1_propagation.analytic import X_GROUP_FLAG, AnalyticPropagation
from tss.l1_propagation.analytic import group_range_m as field_free_group_range_m
from tss.l4_retrieval.lamination import (
    LaminationRetrieval,
    MagnetoionicLaminationRetrieval,
    constant_field,
    geomagnetic_field_model,
    lamina_thicknesses_m,
    lamina_thicknesses_magnetoionic_m,
    lamination_matrix,
    magnetoionic_lamination_row,
    virtual_range_from_laminae_magnetoionic_m,
)
from tss.l5_assessment.requirements_check import check_requirements
from tss.l5_assessment.scoring import profile_metrics, profile_scoring_basis

ALT_S = 8.0e5
F0 = 1.0e6
SCALE_H = 1.5e5
F_H = 0.9e6


# ---------------------------------------------------------------------------
# The index and the group index
# ---------------------------------------------------------------------------


def test_index_zeros_are_the_reflection_levels():
    assert float(mi.refractive_index_squared(1.0, 0.4, 0.6, "O")) == 0.0
    assert float(mi.refractive_index_squared(0.6, 0.4, 0.6, "X")) == pytest.approx(0.0, abs=1e-12)
    # No field: both waves are the isotropic one.
    assert float(mi.refractive_index_squared(0.3, 0.0, 0.6, "O")) == pytest.approx(0.7)
    # Field across the wave normal: the ordinary wave does not see it.
    assert float(mi.refractive_index_squared(0.3, 0.5, np.pi / 2, "O")) == pytest.approx(0.7)


def test_group_index_is_the_frequency_derivative_of_f_mu():
    rng = np.random.default_rng(0)
    checked = 0
    for mode in ("O", "X"):
        for _ in range(1500):
            f = rng.uniform(0.4, 10.0)
            f_h = rng.uniform(0.3, 1.5)
            theta = rng.uniform(0.02, 1.55)
            f_n = rng.uniform(0.05, 1.0) * f
            if not float(mi.refractive_index_squared((f_n / f) ** 2, f_h / f, theta, mode)) > 1e-3:
                continue

            def f_mu(ff, f_n=f_n, f_h=f_h, theta=theta, mode=mode):
                return ff * float(mi.refractive_index((f_n / ff) ** 2, f_h / ff, theta, mode))

            numeric = (f_mu(f * (1 + 1e-6)) - f_mu(f * (1 - 1e-6))) / (2e-6 * f)
            if not np.isfinite(numeric):
                continue
            exact = float(mi.group_index((f_n / f) ** 2, f_h / f, theta, mode))
            assert exact == pytest.approx(numeric, rel=2e-5)
            checked += 1
    assert checked > 1500


def test_group_index_limits():
    # No field: 1 / mu.
    assert float(mi.group_index(0.5, 0.0, 0.3, "O")) == pytest.approx(1.0 / np.sqrt(0.5))
    # Near reflection: cosec(theta) / sqrt(1 - X), Budden's equation 12.12.
    for theta in np.radians([10.0, 30.0, 70.0]):
        eps = 1.0e-6
        assert float(mi.group_index(1.0 - eps, 0.4, theta, "O")) * np.sqrt(eps) == pytest.approx(
            1.0 / np.sin(theta), rel=1e-3
        )
    # Evanescent: no group index.
    assert np.isnan(float(mi.group_index(1.2, 0.4, 0.6, "O")))


def test_index_with_collisions_reduces_to_the_collisionless_one():
    n = mi.appleton_hartree(5.0e6, 3.0e6, 1.0e6, 0.5, 0.0, "O")
    assert float((n**2).real) == pytest.approx(
        float(mi.refractive_index_squared(0.36, 0.2, 0.5, "O"))
    )
    lossy = mi.appleton_hartree(5.0e6, 3.0e6, 1.0e6, 0.5, 1.0e3, "X")
    assert lossy.real > 0.0 and lossy.imag < 0.0


# ---------------------------------------------------------------------------
# The forward model
# ---------------------------------------------------------------------------


def _state(theta_deg: float | None, *, f_h_hz: float = F_H) -> IonosphereState:
    """Exponential topside with a uniform field at ``theta_deg`` from vertical."""
    alt = np.linspace(2.0e5, 1.0e6, 3201)
    f_n = F0 * np.exp((ALT_S - alt) / (2.0 * SCALE_H))
    b = None
    if theta_deg is not None:
        magnitude = f_h_hz / 2.799249e10
        t = np.radians(theta_deg)
        b = np.vstack([np.full(alt.size, magnitude * np.sin(t)), np.zeros(alt.size),
                       np.full(alt.size, magnitude * np.cos(t))])
    return IonosphereState(
        alt_m=alt, ne=density_from_plasma_frequency(f_n), nu=np.full(alt.size, 100.0), b_field_t=b
    )


def _platform(lat: float = 0.0) -> PlatformState:
    return PlatformState(
        alt_m=ALT_S, lat_deg=lat, lon_deg=0.0,
        velocity_ned_mps=np.array([7600.0, 0.0, 0.0]), epoch="2026-01-01T00:00:00Z",
    )


def test_field_at_reads_gyrofrequency_and_angle():
    f_h, theta = mi.field_at(_state(30.0), np.array([5.0e5]))
    assert float(f_h[0]) == pytest.approx(F_H, rel=1e-4)
    assert float(np.degrees(theta[0])) == pytest.approx(30.0, abs=1e-6)
    with pytest.raises(ValueError):
        mi.field_at(_state(None), np.array([5.0e5]), default_b_tesla=3.0e-5)


def test_group_range_equals_the_field_free_one_when_the_field_is_horizontal():
    state = _state(90.0)
    h_refl = float(mi.reflection_altitude_m(state, 2.0e6, ALT_S, "O"))
    with_field = mi.group_range_m(state, 2.0e6, h_refl, ALT_S, "O")
    assert with_field == pytest.approx(field_free_group_range_m(state, 2.0e6, h_refl, ALT_S), rel=1e-9)


def test_ordinary_group_range_grows_as_the_field_steepens():
    h_refl = float(mi.reflection_altitude_m(_state(90.0), 1.5e6, ALT_S, "O"))
    ranges = [
        mi.group_range_m(_state(t), 1.5e6, h_refl, ALT_S, "O") for t in (90.0, 50.0, 30.0, 10.0)
    ]
    assert all(b > a for a, b in zip(ranges, ranges[1:], strict=False))
    assert ranges[0] == pytest.approx(2.0 * SCALE_H * np.arccosh(1.5), rel=1e-3)
    assert ranges[-1] - ranges[0] > 4.0e4


def test_analytic_tier_uses_the_field_when_it_has_one():
    freqs = np.array([1.5e6, 3.0e6])
    rng = np.random.default_rng(0)

    tier = AnalyticPropagation()
    with_field = tier.echoes(_state(30.0), _platform(), freqs, rng)
    assert tier.group_index_used == "magnetoionic"
    assert "group_index=magnetoionic" in with_field.source
    assert X_GROUP_FLAG not in set(with_field.flags)

    plain = AnalyticPropagation(group_index="field_free")
    without = plain.echoes(_state(30.0), _platform(), freqs, rng)
    assert plain.group_index_used == "field_free"
    assert X_GROUP_FLAG in set(without.flags)
    o = np.asarray(with_field.mode) == "O"
    assert np.all(with_field.group_delay_s[o] > without.group_delay_s[o])

    no_field = AnalyticPropagation()
    no_field.echoes(_state(None), _platform(), freqs, rng)
    assert no_field.group_index_used == "field_free"
    with pytest.raises(ValueError):
        AnalyticPropagation(group_index="magnetoionic").echoes(_state(None), _platform(), freqs, rng)


# ---------------------------------------------------------------------------
# The inversion
# ---------------------------------------------------------------------------


def _trace(theta_deg: float, freq: np.ndarray) -> np.ndarray:
    """Ordinary-mode virtual range of the exponential topside, by quadrature."""
    out = []
    for f_k in freq:
        n = 8000
        p0 = np.arcsin(F0 / f_k)
        phi = p0 + (np.pi / 2 - p0) * (np.arange(n) + 0.5) / n
        mu_g = mi.group_index(np.sin(phi) ** 2, F_H / f_k, np.radians(theta_deg), "O")
        out.append(np.sum(mu_g * 2.0 * SCALE_H / np.tan(phi)) * (np.pi / 2 - p0) / n)
    return np.array(out)


def _grid() -> tuple[np.ndarray, np.ndarray]:
    freq = F0 * 1.05 ** np.arange(1, 80)
    freq = freq[freq <= F0 * np.exp((ALT_S - 3.0e5) / (2.0 * SCALE_H))]
    return freq, 2.0 * SCALE_H * np.log(freq / F0)


def test_coefficients_reduce_to_the_arcsine_form_without_a_field():
    freq, _ = _grid()
    k = 12
    edges = np.concatenate(([F0], freq))
    row = magnetoionic_lamination_row(
        freq[k], edges[: k + 1], edges[1 : k + 2], np.zeros(k + 1), np.full(k + 1, 0.5)
    )
    assert row == pytest.approx(lamination_matrix(freq, F0)[k, : k + 1], rel=1e-9)


@pytest.mark.parametrize(
    ("theta_deg", "field_free_error_km"), [(10.0, 32.0), (30.0, 24.0), (50.0, 13.0)]
)
def test_magnetoionic_lamination_recovers_what_the_field_free_one_misses(
    theta_deg, field_free_error_km
):
    freq, depth = _grid()
    trace = _trace(theta_deg, freq)
    field = constant_field(F_H, np.radians(theta_deg))

    plain = np.cumsum(lamina_thicknesses_m(trace, freq, F0))
    corrected = np.cumsum(lamina_thicknesses_magnetoionic_m(trace, freq, F0, ALT_S, field))

    assert np.max(np.abs(plain - depth)) / 1e3 == pytest.approx(field_free_error_km, abs=1.5)
    assert np.max(np.abs(corrected - depth)) < 3.0e3


def test_magnetoionic_round_trip_is_exact():
    freq, _ = _grid()
    field = constant_field(F_H, np.radians(25.0))
    trace = _trace(25.0, freq)
    thickness = lamina_thicknesses_magnetoionic_m(trace, freq, F0, ALT_S, field)
    back = virtual_range_from_laminae_magnetoionic_m(thickness, freq, F0, ALT_S, field)
    assert back == pytest.approx(trace, abs=1e-6)


def test_field_model_follows_the_dipole():
    field, provenance = geomagnetic_field_model(60.0, 0.0, "2026-01-01T00:00:00Z", 8.0e5, prefer="dipole")
    f_h, theta = field(np.array([3.0e5, 8.0e5]))
    assert provenance["model"] == "centred dipole"
    assert f_h[0] > f_h[1] > 0.0
    assert np.degrees(theta[0]) < 30.0


def _ionogram(freq: np.ndarray, virtual: np.ndarray, modes=("unresolved",)) -> Ionogram:
    return Ionogram(
        freq_hz=freq,
        virtual_range_m=virtual[:, None],
        snr_db=np.full((freq.size, 1), 25.0),
        modes=modes,
        flags={},
    )


def test_tier_records_the_field_it_used_and_refuses_an_x_column():
    freq, depth = _grid()
    theta = 30.0
    tier = MagnetoionicLaminationRetrieval(
        start_plasma_hz=F0, field=constant_field(F_H, np.radians(theta))
    )
    profile = tier.invert(_ionogram(freq, _trace(theta, freq)), _platform())
    assert profile.flags["magnetoionic"] is True
    assert profile.flags["field_model"]["field_angle_from_vertical_deg"] == pytest.approx(theta)
    retrieved_depth = ALT_S - profile.alt_m[::-1][1:]
    assert np.max(np.abs(retrieved_depth - depth)) < 3.0e3
    with pytest.raises(ValueError, match="ordinary"):
        tier.invert(_ionogram(freq, _trace(theta, freq), modes=("X",)), _platform())


# ---------------------------------------------------------------------------
# The scoring
# ---------------------------------------------------------------------------


def _profile(flags: dict) -> DensityProfile:
    alt = np.linspace(4.0e5, 7.5e5, 20)
    f_n = F0 * np.exp((ALT_S - alt) / (2.0 * SCALE_H))
    return DensityProfile(alt_m=alt, ne=density_from_plasma_frequency(f_n), flags=flags)


@pytest.mark.parametrize(
    ("flags", "scorable"),
    [
        ({"magnetoionic": True, "forward_group_index": "magnetoionic"}, True),
        ({"magnetoionic": False, "forward_group_index": "magnetoionic"}, False),
        ({"magnetoionic": False, "forward_group_index": "field_free"}, False),
        ({"magnetoionic": True, "forward_group_index": "field_free"}, False),
        ({"magnetoionic": True}, True),
        ({}, True),
    ],
)
def test_profile_is_scored_only_when_the_field_is_treated_consistently(flags, scorable):
    profile = _profile(flags)
    assert profile_scoring_basis(profile)[0] is scorable

    metrics = profile_metrics(profile, _state(30.0))
    # The errors are computed either way, so a reader can still see them.
    assert metrics["height_error_rms_m"] is not None
    records = check_requirements(metrics)
    for rid in ("SCI-05", "SCI-07"):
        if scorable:
            assert records[rid]["verdict"] != "not_evaluated"
        else:
            assert records[rid]["verdict"] == "not_evaluated"
            assert "field" in records[rid]["reason"]


def test_field_free_tier_is_never_scored_against_the_profile_requirements():
    freq, _ = _grid()
    profile = LaminationRetrieval(start_plasma_hz=F0).invert(
        _ionogram(freq, _trace(30.0, freq)), _platform()
    )
    assert profile.flags["magnetoionic"] is False
    assert profile_scoring_basis(profile)[0] is False


# ---------------------------------------------------------------------------
# The study that holds the numbers the Technical Volume quotes
# ---------------------------------------------------------------------------


def test_field_inversion_study_is_registered_and_its_inversion_table_holds():
    from tss.studies import available_studies, field_inversion

    assert available_studies()["field_inversion"]["status"] == "implemented"
    table = field_inversion.inversion_only()
    assert table[:, 0] == pytest.approx(field_inversion.ANGLES_DEG)
    # The field-free error falls as the field tips toward the horizontal.
    assert np.all(np.diff(table[:, 1]) < 0.0)
    assert table[1, 1] / 1e3 == pytest.approx(32.3, abs=1.0)
    assert np.max(table[:, 2]) < 3.0e3
