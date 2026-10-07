"""Layer: L5
WP: WP0, WP7
Claims: C1, C2, C3, C4, C5
Status: implemented
Requirements: INS-01a, INS-01b, INS-02, INS-04, INS-05, INS-07, INS-08, INS-09, INS-10, INS-11, INS-14, INS-15, INS-18, INS-23, INS-24, INS-26, SCI-05, SCI-07, SCI-09, SCI-23

Map computed metrics onto requirement IDs and score them.

The mapping is a table in this module and not a set of conditionals scattered
through the assessment tier, because the mapping is itself a claim: it asserts
that a particular measured quantity is what a particular requirement means.
That assertion should be readable in one place, reviewable, and changeable
without touching the scoring arithmetic.

The rule that matters most here is the second one.  A requirement whose metric
was not computed is reported ``not_evaluated``, with the metric name and the
reason attached.  It is never reported as met.  The failure this prevents is
specific: a campaign that does not exercise a capability produces no metric for
it, and a scoring routine that skipped such requirements would print a table of
passes that looked like a clean sweep.  :meth:`tss.requirements.RequirementSet.score`
already reports the full set with absent entries as ``not_evaluated``; this
module adds the reason, so that a reader can tell a capability that was not
attempted from one that was attempted and not measured.

The mapping table
-----------------

Each row is (requirement ID, metric name, note).  A row with a metric of
``None`` is a requirement this layer deliberately does not score, and the note
says who owns it or why it cannot be scored as written.

==========  ================================  ==========================================
Requirement Metric                             Note
==========  ================================  ==========================================
INS-01a     sweep_min_frequency_hz             Lowest frequency swept (sense max).
INS-01b     sweep_max_frequency_hz             Highest frequency swept (sense min).
INS-02      frequency_step_fraction_max        df / min(f, 3 MHz), so the fractional
                                               rule and the 150 kHz cap are one number.
INS-04      unambiguous_range_m                From the sample rate and the sweep rate.
INS-05      virtual_range_resolution_by_       Window ENBW over the coherent block,
            frequency (falls back to           per column, scored per band (Rev G,
            virtual_range_resolution_m)        amended 2026-09-30).
INS-07      ionogram_period_s                  The sweep period actually transmitted.
INS-08      mode_discrimination                Boolean.  False from a scalar receiver,
                                               and False is a failure, not a silence.
INS-09      residual_range_bias_cells          Residual after the WP5 separation, in
                                               INS-05 cells.  None until WP5 exists.
INS-10      snr_db_min                         Worst per-cell SNR over the detected
                                               trace.
INS-11      none                               Receiver dynamic range, owned by L2.
INS-14      none                               Transmit power, owned by L2 and by the
                                               configuration.
INS-15      coherent_interval_s                Only a tier that integrates coherently
                                               reports one.  See "INS-15" below.
INS-18      none                               Vector channel phase match, owned by L2.
INS-23      none                               Transmit to receive isolation, owned by
                                               L2 and measured on the bench in WP9.
INS-24      none                               Residual leakage floor, same owners.
INS-26      none                               Oscillator stability, owned by L2.
SCI-05      ne_rms_fractional                  Density error against truth (SW-01).
                                               Scored only when the field is treated
                                               consistently; see "SCI-05 and SCI-07".
SCI-07      height_error_rms_m                 Height error against truth (SW-01).
                                               Same condition.
SCI-09      along_track_sampling_m             Ionogram period times horizontal speed.
SCI-23      drift_resolution_mps               Drift accuracy from a second look (Rev G).
                                               None until WP5 exists; rename pending.
==========  ================================  ==========================================

INS-01 (split into INS-01a and INS-01b in Rev G)
-------------------------------------------------
Up to Rev F, INS-01 was the sounding band, carried with ``sense: range`` and a threshold of
[0.33, 20] MHz against a goal of [0.26, 20] MHz.  The ``range`` sense in
:mod:`tss.requirements` compares one achieved value against those bounds, and
the goal band contains the threshold band, so any achieved value inside the
threshold band is also inside the goal band and would be scored ``goal``.  A
sounder whose lowest usable frequency was 5 MHz would be scored ``goal`` on the
band requirement, which is not what the requirement means.

Up to Rev F, INS-01 was therefore left ``not_evaluated``.  Rev G splits it
into two one-sided requirements: INS-01a, the lowest sounding frequency with
``sense: max``, and INS-01b, the highest, with ``sense: min``.  Both are
scored from the band the transmitted sweep covered, ``sweep_min_frequency_hz``
and ``sweep_max_frequency_hz``, because INS-01 is a requirement on what the
instrument can sweep.  The first Rev G mapping scored them from the lowest and
highest *detected* frequencies instead.  Those are set by the ionosphere (the
plasma frequency at the spacecraft and fxF2), so INS-01a failed in any
ionosphere denser than 0.33 MHz at the spacecraft whatever the instrument
could do.  The detected edges are still carried in ``Assessment.errors`` as
diagnostics.

SCI-05 and SCI-07
-----------------
Both are scored from the L4 profile against its L0 truth, and both carry a
precondition: the echoes must have been made with the magnetoionic group index
and the inversion must have corrected for it.  The field changes the group
path of the ordinary wave by enough to move a topside profile by 14 to 34 km
at mid and high magnetic latitude, against the 20 km of SCI-07, so a field-free
inversion fails for a reason that is known and is not the instrument's, and a
field-free inversion of field-free echoes passes for a reason that says nothing
about a real ionosphere.  In either case the verdict is ``not_evaluated`` with
the reason attached, and the error metrics are still carried in
``Assessment.errors``.  :func:`tss.l5_assessment.scoring.profile_scoring_basis`
is the rule.

INS-15
------
INS-15 asks for a coherent processing interval of 0.05 s at threshold, from
the SCI-15 accuracy, and 1.15 s at goal, from the SCI-15 resolution, at 3 MHz
and at twice the plasma frequency at the spacecraft (Rev G, amended
2026-10-06 for the refractive index at the spacecraft).
It is mapped to ``coherent_interval_s``, which only a tier that performs
motion-compensated coherent integration reports.  The baseline tier's
per-transform block duration is carried separately as ``block_duration_s`` and
is deliberately not scored against INS-15.  Scoring it would report a failure
of a capability the baseline does not attempt, and INS-15 exists to constrain
the C3 processing of WP5, not the conventional chain of WP4.  The distinction
between "did not meet" and "did not attempt" is exactly what ``not_evaluated``
is for.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from tss.requirements import RequirementSet, load_requirements

__all__ = [
    "MAPPING_TABLE",
    "FIELD_CONSISTENT_ONLY",
    "METRIC_FOR_REQUIREMENT",
    "NOT_SCORED_HERE",
    "CARRY_IN_ERRORS",
    "achieved_from_metrics",
    "check_requirements",
]

#: The table above, as data.  ``(requirement_id, metric_name_or_None, note)``.
#: Everything else in this module is derived from it, so the table and the
#: behaviour cannot drift apart.
MAPPING_TABLE: tuple[tuple[str, str | None, str], ...] = (
    (
        "INS-01a",
        "sweep_min_frequency_hz",
        "lowest frequency the transmitted sweep covered, scored with sense max.  "
        "Not the lowest detected echo, which the ionosphere sets.",
    ),
    (
        "INS-01b",
        "sweep_max_frequency_hz",
        "highest frequency the transmitted sweep covered, scored with sense min.  "
        "Not the highest detected echo, which is fxF2.",
    ),
    (
        "INS-02",
        "frequency_step_fraction_max",
        "df / min(f, 150 kHz / 0.05), which expresses the fractional rule and the "
        "absolute cap as one fraction scored against 0.05.",
    ),
    ("INS-04", "unambiguous_range_m", "c fs / (4 K_f), from the record and the sweep."),
    (
        "INS-05",
        "virtual_range_resolution_by_frequency",
        "window equivalent noise bandwidth over the coherent block duration, not the "
        "Fourier bin spacing, which zero padding would otherwise flatter.  Scored per "
        "band: <=15 km below 0.5 MHz and <=10 km at and above at threshold.  A tier that "
        "reports only one resolution is scored on it against the strictest band.",
    ),
    ("INS-07", "ionogram_period_s", "the sweep period the instrument actually transmitted."),
    (
        "INS-08",
        "mode_discrimination",
        "boolean.  A scalar receive chain reports False and is scored 'fail', which "
        "is correct: the baseline configuration does not meet INS-08 and the vector "
        "receive of WP11 is what would.",
    ),
    (
        "INS-09",
        "residual_range_bias_cells",
        "residual apparent-range bias after the WP5 separation, in INS-05 range "
        "cells, threshold and goal both 1.0.  A tier that has not separated Doppler "
        "reports nothing and is scored 'not_evaluated'.",
    ),
    (
        "INS-10",
        "snr_db_min",
        "worst per-range-cell signal to noise ratio over the detected trace.  The "
        "cell noise bandwidth is the one from tss.l3_signal.transform, about 37 Hz "
        "for a 1 km cell at K_f = 5.6 MHz/s, which is what makes it closable at "
        "sub-watt transmit power.",
    ),
    (
        "INS-11",
        None,
        "instantaneous receiver dynamic range is an L2 property; it is measured by "
        "the receiver model and the bench work of WP9, not by the processing.",
    ),
    (
        "INS-14",
        None,
        "transmit power is set by the configuration and reported by L2.  L5 does not "
        "re-derive a configured value as though it had measured it.",
    ),
    (
        "INS-15",
        "coherent_interval_s",
        "only reported by a tier that integrates coherently.  The baseline tier's "
        "per-block duration is carried as block_duration_s and is not scored here; "
        "see the module docstring.",
    ),
    ("INS-18", None, "vector receive channel phase match is an L2 property, WP11."),
    (
        "INS-23",
        None,
        "transmit to receive isolation is an L2 property and ultimately a WP9 bench "
        "measurement.  docs/verification.md states that C4 is a hardware claim with "
        "a modelling component.",
    ),
    ("INS-24", None, "residual leakage floor, same owners as INS-23."),
    ("INS-26", None, "oscillator Allan deviation is an L2 property, WP3."),
    (
        "SCI-05",
        "ne_rms_fractional",
        "rms fractional density error of the L4 profile against its L0 truth, the "
        "simulated form of SCI-05's absolute accuracy (SW-01).",
    ),
    (
        "SCI-07",
        "height_error_rms_m",
        "rms height error of the L4 profile points against their L0 truth (SW-01).",
    ),
    (
        "SCI-09",
        "along_track_sampling_m",
        "ionogram period times the platform's horizontal speed.  Sampling, not "
        "resolution; see tss.l5_assessment.scoring.",
    ),
    (
        "SCI-23",
        "drift_resolution_mps",
        "achieved drift accuracy from a second look (Rev G).  None until WP5 "
        "exists, so 'not_evaluated' rather than a failure.",
    ),
)

#: Requirement ID to metric name, for the rows that have a metric.
METRIC_FOR_REQUIREMENT: dict[str, str] = {
    rid: metric for rid, metric, _ in MAPPING_TABLE if metric is not None
}

#: Requirements scored from the retrieved profile, which need the field
#: treated consistently by the echoes and the inversion.
FIELD_CONSISTENT_ONLY: tuple[str, ...] = ("SCI-05", "SCI-07")

#: A metric to use when the mapped one was not computed.  INS-05 is mapped to
#: the per-column resolution; a tier that reports only the scalar is scored on
#: that against the strictest band.
FALLBACK_METRIC: dict[str, str] = {"INS-05": "virtual_range_resolution_m"}

#: Requirement ID to the reason it is not scored here, for the rows that do not.
NOT_SCORED_HERE: dict[str, str] = {
    rid: note for rid, metric, note in MAPPING_TABLE if metric is None
}

#: Metrics that are reported in ``Assessment.errors`` whether or not any
#: requirement is scored against them.  These are the numbers a reader wants
#: even when the requirement schema cannot express the comparison.
CARRY_IN_ERRORS: tuple[str, ...] = (
    "ne_rms_fractional",
    "ne_bias_fractional",
    "ne_rms_absolute",
    "ne_bias_absolute",
    "height_error_rms_m",
    "height_error_bias_m",
    "inversion_magnetoionic",
    "forward_group_index",
    "n_points_compared",
    "n_height_points_compared",
    "min_sounding_frequency_hz",
    "max_sounding_frequency_hz",
    "sweep_min_frequency_hz",
    "sweep_max_frequency_hz",
    "virtual_range_resolution_m",
    "unambiguous_range_m",
    "noise_bandwidth_hz",
    "snr_db_min",
    "snr_db_median",
    "along_track_sampling_m",
    "block_duration_s",
    "n_detected",
)


def achieved_from_metrics(metrics: Mapping[str, Any]) -> dict[str, Any]:
    """Requirement ID to achieved value, for the metrics that were computed.

    A metric that is absent or ``None`` produces no entry, which is what makes
    :meth:`tss.requirements.RequirementSet.score` report the requirement as
    ``not_evaluated``.  A metric that is ``False`` does produce an entry, because
    False is a measurement: the capability was looked for and was not there.
    """
    achieved: dict[str, Any] = {}
    for rid, metric in METRIC_FOR_REQUIREMENT.items():
        if rid in FIELD_CONSISTENT_ONLY and metrics.get("profile_scorable") is False:
            continue
        value = metrics.get(metric)
        if value is None and rid in FALLBACK_METRIC:
            value = metrics.get(FALLBACK_METRIC[rid])
        if value is None:
            continue
        achieved[rid] = value
    return achieved


def check_requirements(
    metrics: Mapping[str, Any],
    requirements: RequirementSet | str | None = None,
) -> dict[str, dict[str, Any]]:
    """Score a metric mapping against the requirement set.

    Parameters
    ----------
    metrics
        Output of :func:`tss.l5_assessment.scoring.compute_metrics`, or any
        mapping with the same keys.
    requirements
        A loaded :class:`~tss.requirements.RequirementSet`, a name to load one
        from, or ``None`` for the default set.

    Returns
    -------
    dict
        One record per requirement in the set, as
        :meth:`tss.requirements.RequirementSet.score` produces, with two extra
        keys added: ``metric``, naming the metric the requirement was scored
        from or ``None``, and ``reason``, present only on a ``not_evaluated``
        verdict and saying why.

    Notes
    -----
    Every requirement in the set appears in the result, including those no
    metric is mapped to.  That is the point: the table of verdicts is the full
    requirement set every time, so a campaign's coverage is legible from its
    own output rather than from what is missing from it.
    """
    if requirements is None:
        req_set = load_requirements()
    elif isinstance(requirements, RequirementSet):
        req_set = requirements
    else:
        req_set = load_requirements(requirements)

    scored = req_set.score(achieved_from_metrics(metrics))

    for rid, record in scored.items():
        metric = METRIC_FOR_REQUIREMENT.get(rid)
        if metric is not None and metrics.get(metric) is None and rid in FALLBACK_METRIC:
            if metrics.get(FALLBACK_METRIC[rid]) is not None:
                metric = FALLBACK_METRIC[rid]
        record["metric"] = metric
        if record.get("verdict") != "not_evaluated":
            continue
        if rid in FIELD_CONSISTENT_ONLY and metrics.get("profile_scorable") is False:
            record["reason"] = str(
                metrics.get("profile_not_scorable_reason")
                or "the magnetic field was not treated consistently"
            )
        elif metric is None:
            record["reason"] = NOT_SCORED_HERE.get(
                rid,
                "no metric is mapped to this requirement in "
                "tss.l5_assessment.requirements_check.MAPPING_TABLE",
            )
        elif metric not in metrics:
            record["reason"] = (
                f"metric {metric!r} is not in the metric set this run produced"
            )
        else:
            record["reason"] = (
                f"metric {metric!r} was not computed by this run; the tier that "
                "would produce it did not run or did not report it"
            )
    return scored
