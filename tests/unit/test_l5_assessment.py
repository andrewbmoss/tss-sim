"""Unit tests for L5, assessment.

The tests that matter most here are not about arithmetic.  They are about the
discipline of Section 3.1 of the proposal: a requirement whose metric was not
computed must be reported ``not_evaluated`` and must never be reported as met.
A scoring routine that quietly dropped unmeasured requirements would print a
short table of passes that reads as a clean sweep, and nothing downstream would
notice.  Several tests here exist only to make that impossible.

The four senses are tested at, above and below both threshold and goal.  The
verdict boundaries are inclusive on both ``min`` and ``max``, so a value
exactly at the threshold is a pass, and that is checked explicitly rather than
left to be discovered.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from tss.constants import density_from_plasma_frequency
from tss.interfaces import DensityProfile, Ionogram, IonosphereState, PlatformState
from tss.l5_assessment import (
    METRIC_FOR_REQUIREMENT,
    StandardAssessment,
    check_requirements,
    compute_metrics,
    feasibility,
    profile_metrics,
    render_json,
    render_text,
    summarise_verdicts,
)
from tss.registry import available
from tss.requirements import Requirement, load_requirements

ALT_S = 5.0e5
F0 = 5.0e5
SCALE_H = 8.0e4


def _truth() -> IonosphereState:
    alt = np.linspace(2.0e5, ALT_S, 400)
    f_n = F0 * np.exp((ALT_S - alt) / (2.0 * SCALE_H))
    return IonosphereState(
        alt_m=alt,
        ne=density_from_plasma_frequency(f_n),
        lat_deg=0.0,
        lon_deg=0.0,
        epoch="2026-01-01T00:00:00Z",
        source="analytic exponential topside",
    )


def _retrieved(offset_m: float = 0.0, scale: float = 1.0) -> DensityProfile:
    alt = np.linspace(2.5e5, ALT_S, 60)
    f_n = F0 * np.exp((ALT_S - (alt - offset_m)) / (2.0 * SCALE_H))
    return DensityProfile(
        alt_m=alt,
        ne=density_from_plasma_frequency(f_n) * scale,
        flags={"starting_value": "truth"},
    )


def _ionogram(**flags) -> Ionogram:
    # 100 points keeps df / min(f, 3 MHz) at 0.024, inside the INS-02 goal.
    freq = F0 * np.exp(np.linspace(0.02, 2.0, 100))
    base = {
        "range_resolution_m": 1.8e3,
        "unambiguous_range_m": 2.6e6,
        "ionogram_period_s": 3.3,
        "noise_bandwidth_hz": 67.2,
        "mode_discrimination": False,
    }
    base.update(flags)
    return Ionogram(
        freq_hz=freq,
        virtual_range_m=(2.0 * SCALE_H * np.arccosh(freq / F0))[:, None],
        snr_db=np.full((freq.size, 1), 25.0),
        modes=("unresolved",),
        flags=base,
    )


def _platform() -> PlatformState:
    return PlatformState(
        alt_m=ALT_S,
        lat_deg=0.0,
        lon_deg=0.0,
        velocity_ecef_mps=np.array([7600.0, 0.0, 0.0]),
        epoch="2026-01-01T00:00:00Z",
    )


# ---------------------------------------------------------------------------
# The four senses, at, above and below threshold and goal
# ---------------------------------------------------------------------------


def _req(sense, threshold, goal):
    return Requirement(
        id="TST-01", title="test", sense=sense, threshold=threshold, goal=goal,
        source="test fixture",
    )


@pytest.mark.parametrize(
    "achieved,expected",
    [
        (25.0, "goal"),        # above goal
        (20.0, "goal"),        # exactly at goal
        (15.0, "threshold"),   # between
        (10.0, "threshold"),   # exactly at threshold
        (5.0, "fail"),         # below threshold
    ],
)
def test_min_sense_verdicts(achieved, expected):
    """``min`` means at or above.  The boundaries are inclusive."""
    assert _req("min", 10.0, 20.0).evaluate(achieved)["verdict"] == expected


@pytest.mark.parametrize(
    "achieved,expected",
    [
        (1.0, "goal"),         # below goal
        (5.0e3, "goal"),       # exactly at goal
        (8.0e3, "threshold"),  # between
        (1.0e4, "threshold"),  # exactly at threshold
        (2.0e4, "fail"),       # above threshold
    ],
)
def test_max_sense_verdicts(achieved, expected):
    """``max`` means at or below.  The boundaries are inclusive."""
    assert _req("max", 1.0e4, 5.0e3).evaluate(achieved)["verdict"] == expected


@pytest.mark.parametrize(
    "achieved,expected",
    [
        (5.0, "goal"),         # inside the goal band
        (4.0, "goal"),         # on the goal lower edge
        (6.0, "goal"),         # on the goal upper edge
        (8.0, "threshold"),    # inside the threshold band only
        (0.0, "threshold"),    # on the threshold lower edge
        (12.0, "fail"),        # outside both
        (-1.0, "fail"),        # outside both, the other way
    ],
)
def test_range_sense_verdicts(achieved, expected):
    """``range`` means inside ``[lower, upper]``, edges included.

    Scored with a goal band nested inside the threshold band, which is what the
    sense is designed for: outside the threshold band is a failure, inside it
    is at least a threshold pass, and inside the narrower goal band is a goal.
    """
    req = _req("range", [0.0, 10.0], [4.0, 6.0])
    assert req.evaluate(achieved)["verdict"] == expected


@pytest.mark.parametrize(
    "achieved,expected",
    [
        (3.3e5, "goal"),        # on the threshold edge, and so inside the goal band
        (1.0e6, "goal"),        # far worse, and still scored goal
        (3.0e5, "fail"),        # inside the goal band but below the threshold band
    ],
)
def test_range_sense_is_the_wrong_tool_for_band_coverage(achieved, expected):
    """The INS-01 case, checked so that the reason it is not scored is on record.

    INS-01's goal band *contains* its threshold band, because a wider band is
    better.  The ``range`` sense checks the threshold band first, so any value
    inside it is also inside the goal band and scores ``goal``: a sounder whose
    lowest usable frequency was 1 MHz would score ``goal`` on the sounding
    band.  The sense is behaving correctly for what it means; it is the wrong
    tool for coverage.  That is why
    :mod:`tss.l5_assessment.requirements_check` leaves INS-01 unevaluated and
    records the achieved band edges in ``Assessment.errors`` instead.
    """
    req = _req("range", [3.3e5, 2.0e7], [2.6e5, 2.0e7])
    assert req.evaluate(achieved)["verdict"] == expected


@pytest.mark.parametrize("achieved,expected", [(True, "goal"), (False, "fail")])
def test_bool_sense_verdicts(achieved, expected):
    """``bool`` is demonstrated or not.  False is a failure, not a silence."""
    assert _req("bool", True, None).evaluate(achieved)["verdict"] == expected


@pytest.mark.parametrize("sense", ["min", "max", "range", "bool"])
def test_every_sense_reports_not_evaluated_for_a_missing_value(sense):
    """``None`` is ``not_evaluated`` in all four senses, never a pass."""
    threshold = [1.0, 2.0] if sense == "range" else 1.0
    record = _req(sense, threshold, None).evaluate(None)
    assert record["verdict"] == "not_evaluated"
    assert record["margin"] is None


# ---------------------------------------------------------------------------
# An absent metric is never a silent pass
# ---------------------------------------------------------------------------


def test_empty_metrics_evaluate_nothing_and_pass_nothing():
    """The check that protects the whole discipline.

    Run with no metrics at all, every requirement must be ``not_evaluated``,
    none may be met, and every row must carry a reason.  A reason-less
    ``not_evaluated`` is nearly as bad as a silent pass, because it gives the
    reader nothing to act on.
    """
    records = check_requirements({})
    assert records, "the requirement set should not be empty"
    assert all(r["verdict"] == "not_evaluated" for r in records.values())
    assert all(r.get("reason") for r in records.values())


def test_a_metric_that_was_not_computed_names_itself_in_the_reason():
    """INS-09 with no Doppler estimate says which metric is missing."""
    records = check_requirements(compute_metrics(
        _retrieved(), _truth(), _ionogram(), platform=_platform()
    ))
    ins09 = records["INS-09"]
    assert ins09["verdict"] == "not_evaluated"
    assert ins09["metric"] == "residual_range_bias_cells"
    assert "residual_range_bias_cells" in ins09["reason"]


def test_a_requirement_owned_by_another_layer_says_so():
    """INS-23 is an L2 property and a WP9 bench measurement, and says both."""
    records = check_requirements({})
    assert records["INS-23"]["metric"] is None
    assert "L2" in records["INS-23"]["reason"]


def test_every_mapped_metric_name_is_produced_by_compute_metrics():
    """The mapping table cannot name a metric the scorer never emits.

    A row that named a metric nothing produces would report ``not_evaluated``
    for ever and nobody would find out why.
    """
    metrics = compute_metrics(_retrieved(), _truth(), _ionogram(), platform=_platform())
    missing = sorted(set(METRIC_FOR_REQUIREMENT.values()) - set(metrics))
    assert not missing, f"mapped but never produced: {missing}"


def test_every_mapped_requirement_id_exists_in_the_requirement_set():
    """The mapping table cannot name a requirement that is not in the YAML."""
    known = set(load_requirements())
    assert set(METRIC_FOR_REQUIREMENT) <= known
    from tss.l5_assessment import NOT_SCORED_HERE

    assert set(NOT_SCORED_HERE) <= known


# ---------------------------------------------------------------------------
# The metrics themselves
# ---------------------------------------------------------------------------


def test_a_perfect_retrieval_scores_zero_error():
    """Truth against itself gives zero RMS and zero bias in both forms."""
    truth = _truth()
    exact = DensityProfile(alt_m=truth.alt_m.copy(), ne=truth.ne.copy())
    metrics = profile_metrics(exact, truth)
    assert metrics["ne_rms_fractional"] == pytest.approx(0.0, abs=1.0e-12)
    assert metrics["ne_bias_fractional"] == pytest.approx(0.0, abs=1.0e-12)
    assert metrics["height_error_rms_m"] == pytest.approx(0.0, abs=1.0e-6)


def test_a_uniformly_scaled_retrieval_shows_the_scale_as_a_bias():
    """A 10 per cent density overestimate reads as a 0.1 fractional bias.

    Built on the truth grid so that no interpolation error enters and the
    expected number is exactly 0.1.  On a different grid the linear
    interpolation of an exponential profile contributes a few parts in a
    thousand of its own, which is real and is not what this test is about.
    """
    truth = _truth()
    scaled = DensityProfile(alt_m=truth.alt_m.copy(), ne=truth.ne * 1.1)
    metrics = profile_metrics(scaled, truth)
    assert metrics["ne_bias_fractional"] == pytest.approx(0.1, rel=1.0e-9)
    assert metrics["ne_rms_fractional"] == pytest.approx(0.1, rel=1.0e-9)


def test_a_displaced_retrieval_shows_the_displacement_as_a_height_error():
    """A profile shifted 2 km up reads as a 2 km height bias at matched density.

    This is the metric a topside inversion should be judged on.  The same
    displacement read as a density error would be whatever the local gradient
    made of it, which is a different number in a steep part of the profile than
    in a flat one.
    """
    metrics = profile_metrics(_retrieved(offset_m=2.0e3), _truth())
    assert metrics["height_error_bias_m"] == pytest.approx(2.0e3, rel=5.0e-3)


def test_along_track_sampling_is_period_times_horizontal_speed():
    metrics = compute_metrics(_retrieved(), _truth(), _ionogram(), platform=_platform())
    assert metrics["along_track_sampling_m"] == pytest.approx(3.3 * 7600.0, rel=1.0e-9)


def test_along_track_sampling_is_none_without_a_platform():
    """Missing input gives ``None``, which becomes ``not_evaluated``, not zero."""
    metrics = compute_metrics(_retrieved(), _truth(), _ionogram(), platform=None)
    assert metrics["along_track_sampling_m"] is None
    assert check_requirements(metrics)["SCI-09"]["verdict"] == "not_evaluated"


def test_frequency_step_metric_expresses_both_branches_of_ins02():
    """``df / min(f, 3 MHz)`` is below 0.05 for a compliant grid and above for not.

    The crossover is ``150 kHz / 0.05``.  Below it the fractional rule binds,
    above it the absolute cap does, and dividing by ``min(f, crossover)`` turns
    the two branches into one number that can be scored against 0.05.
    """
    fine = _ionogram()
    assert check_requirements(
        compute_metrics(_retrieved(), _truth(), fine)
    )["INS-02"]["verdict"] in ("goal", "threshold")

    coarse = Ionogram(
        freq_hz=np.array([1.0e6, 1.5e6, 2.0e6]),
        virtual_range_m=np.full((3, 1), 3.0e5),
        snr_db=np.full((3, 1), 20.0),
        modes=("unresolved",),
        flags={"range_resolution_m": 1.8e3},
    )
    records = check_requirements(compute_metrics(_retrieved(), _truth(), coarse))
    assert records["INS-02"]["verdict"] == "fail"


# ---------------------------------------------------------------------------
# The tier and the report
# ---------------------------------------------------------------------------


def test_standard_tier_is_registered_and_returns_an_assessment():
    assert "standard" in available("l5")["l5"]
    tier = StandardAssessment(platform=_platform())
    assessment = tier.score(_retrieved(), _truth(), _ionogram(), "leo_500km_single")
    assert assessment.configuration == "leo_500km_single"
    assert len(assessment.requirements) == len(load_requirements())
    assert assessment.flags["summary"]["fail"] >= 1        # INS-08, scalar receive
    assert "ne_rms_fractional" in assessment.errors


def test_standard_tier_refuses_an_unnamed_configuration():
    """A result without its configuration is not a result."""
    with pytest.raises(ValueError, match="configuration"):
        StandardAssessment().score(_retrieved(), _truth(), _ionogram(), "")


def test_standard_tier_carries_the_flags_that_change_how_numbers_read():
    """Which starting value, and whether the modes were labelled."""
    assessment = StandardAssessment().score(
        _retrieved(), _truth(), _ionogram(), "leo_500km_single"
    )
    assert assessment.flags["l4.starting_value"] == "truth"
    assert assessment.flags["l3.mode_discrimination"] is False


def test_report_shows_every_requirement_including_the_unevaluated_ones():
    """A table that showed only what was scored would read as a clean sweep."""
    assessment = StandardAssessment(platform=_platform()).score(
        _retrieved(), _truth(), _ionogram(), "leo_500km_single"
    )
    text = render_text(assessment)
    for rid in load_requirements():
        assert rid in text
    assert "NOT EVALUATED" in text
    assert "reason:" in text


def test_report_json_round_trips():
    assessment = StandardAssessment().score(
        _retrieved(), _truth(), _ionogram(), "leo_500km_single"
    )
    payload = json.loads(render_json(assessment))
    assert payload["configuration"] == "leo_500km_single"
    assert set(payload["summary"]) == {"goal", "threshold", "fail", "not_evaluated"}


def test_summarise_keeps_zero_counts():
    """A summary line has the same shape every time it is printed."""
    assert summarise_verdicts({}) == {
        "goal": 0, "threshold": 0, "fail": 0, "not_evaluated": 0
    }


# ---------------------------------------------------------------------------
# Stubs
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "call",
    [
        lambda: feasibility.feasible_mask([]),
        lambda: feasibility.binding_requirement([]),
        lambda: feasibility.feasibility_surface([], {}),
    ],
)
def test_feasibility_is_a_stub_naming_wp7(call):
    with pytest.raises(NotImplementedError, match="WP7"):
        call()


def test_profile_metrics_score_the_topside_only():
    """Retrieved points below the truth peak are counted, not scored.

    A topside sounder has no echo from below the F2 peak, so a retrieved point
    there is an artefact.  SW-01 scopes the density and height errors to the
    topside, and a single wild point below the peak must not dominate them.
    """
    alt = np.linspace(1.0e5, ALT_S, 400)
    hm = 3.0e5
    ne = 1.0e12 * np.exp(1.0 - (alt - hm) / SCALE_H - np.exp(-(alt - hm) / SCALE_H))
    truth = IonosphereState(alt_m=alt, ne=ne, lat_deg=0.0, lon_deg=0.0,
                            epoch="2026-01-01T00:00:00Z", source="Chapman layer")
    top = np.linspace(3.2e5, ALT_S, 20)
    good = DensityProfile(alt_m=top, ne=np.interp(top, alt, ne), flags={})
    with_artefact = DensityProfile(
        alt_m=np.concatenate(([1.8e5], top)),
        ne=np.concatenate(([1.0e14], np.interp(top, alt, ne))),
        flags={},
    )
    clean = profile_metrics(good, truth)
    dirty = profile_metrics(with_artefact, truth)
    assert dirty["ne_rms_fractional"] == pytest.approx(clean["ne_rms_fractional"], abs=1e-9)
    assert dirty["n_points_below_peak"] == 1
    assert dirty["overlap_min_alt_m"] >= hm - 2.0e3


def test_band_edges_are_scored_from_the_sweep_not_the_detected_echoes():
    """INS-01a/b ask what the instrument can sweep, not where echoes appeared.

    The lowest detected echo is the plasma frequency at the spacecraft and the
    highest is fxF2, both set by the ionosphere.  Scoring INS-01a against the
    detected edge would fail a sounder that swept down to 0.33 MHz simply
    because the local plasma frequency that day was 3 MHz.
    """
    ionogram = _ionogram(sweep_min_frequency_hz=3.3e5, sweep_max_frequency_hz=2.0e7)
    metrics = compute_metrics(_retrieved(), _truth(), ionogram, platform=_platform())
    # The detected edges of the fixture sit well inside the band.
    assert metrics["min_sounding_frequency_hz"] > 3.3e5
    assert metrics["max_sounding_frequency_hz"] < 2.0e7
    records = check_requirements(metrics)
    assert records["INS-01a"]["metric"] == "sweep_min_frequency_hz"
    assert records["INS-01a"]["verdict"] == "threshold"
    assert records["INS-01b"]["verdict"] in ("threshold", "goal")

    narrow = _ionogram(sweep_min_frequency_hz=9.0e5, sweep_max_frequency_hz=2.0e7)
    records = check_requirements(
        compute_metrics(_retrieved(), _truth(), narrow, platform=_platform())
    )
    assert records["INS-01a"]["verdict"] == "fail"
