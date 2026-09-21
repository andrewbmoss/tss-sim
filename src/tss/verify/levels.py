"""Layer: support
WP: WP8
Claims: C1, C2, C3, C4, C5
Status: partial
Requirements: -

The V1 to V4 verification structure, and the checks registered at each level.

``docs/verification.md`` sets out what each level answers and what it cannot
reach.  This module is that document made executable: a registry of named
checks, a runner that collects their outcomes, and a record type that carries a
verdict, the quantity compared, the observed discrepancy and the tolerance.

    V1  Internal consistency, including closed-form limiting cases.
    V2  Cross-fidelity: two independent implementations of the same physics.
    V3  Literature replication.
    V4  Archived real measurements.

Verdicts
--------
``pass``
    The check ran and agreed within tolerance.
``fail``
    The check ran and did not agree, or raised an unexpected exception.  A fail
    is never suppressed and never resolved by widening a tolerance without a
    recorded justification.
``not_run``
    The check is registered and its physics is not written.  The reason names
    the work package.  This is the verdict that matters most in the present
    state of the codebase, because most of what V2 to V4 would compare does not
    exist yet, and a harness that omitted those checks would report a short
    clean table.
``skipped``
    The check could run in principle but something external is absent: MATLAB,
    PHaRLAP, or a data set that has not been fetched.

``docs/verification.md`` lists four verdicts, ``pass``, ``fail``, ``skipped``
and ``partial``, and does not list ``not_run``.  The distinction being drawn
here is between a check that cannot run on this machine, which is ``skipped``
and might pass on another machine, and a check whose subject has not been
written, which is ``not_run`` and will not pass anywhere.  Collapsing the two
would let the codebase's real coverage hide behind a missing MATLAB install.
The document needs the extra row.

How a check is written
----------------------
A check is a zero-argument function decorated with :func:`check`.  It returns a
mapping which may carry ``ok``, ``observed``, ``tolerance`` and ``detail``; an
absent ``ok`` means it passed.  To report that it could not run, it raises
:class:`NotRun` with a reason, or lets a :class:`NotImplementedError` from a
stub propagate, which is turned into ``not_run`` carrying the stub's own
message.  That last route is deliberate: it means a check written against an
unwritten algorithm reports the work package automatically, and it starts
passing or failing honestly the moment the algorithm lands, with no change to
the check.

Any other exception is a ``fail``.  A check that crashes is a check that found
something.
"""

from __future__ import annotations

import traceback
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any

__all__ = [
    "LEVELS",
    "VERDICTS",
    "CheckResult",
    "NotRun",
    "Skipped",
    "check",
    "registered",
    "run_level",
    "run_all",
]

#: The four levels, in the order a report shows them.
LEVELS: tuple[str, ...] = ("V1", "V2", "V3", "V4")

#: The verdicts a check can carry.  See the module docstring on ``not_run``.
VERDICTS: tuple[str, ...] = ("pass", "fail", "partial", "skipped", "not_run")


class NotRun(Exception):
    """Raised by a check whose subject is not written, carrying the reason.

    The reason names the work package.  A ``not_run`` check will not pass on
    any machine until that work package lands, which is what distinguishes it
    from :class:`Skipped`.
    """

    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = str(reason)


class Skipped(NotRun):
    """Raised by a check that could run but for an absent external dependency.

    MATLAB, PHaRLAP, or a data set that has not been fetched.  A ``skipped``
    check might pass on another machine; a ``not_run`` one will not.  Keeping
    the two apart is what stops the codebase's real coverage hiding behind a
    missing MATLAB install.
    """


@dataclass(frozen=True)
class CheckResult:
    """One verification check and what it found."""

    level: str
    name: str
    verdict: str
    quantity: str = ""
    observed: float | None = None
    tolerance: float | None = None
    detail: str = ""
    reason: str = ""
    reference: str = ""

    def as_dict(self) -> dict[str, Any]:
        """JSON-serialisable record, for the verification report."""
        return {
            "level": self.level,
            "name": self.name,
            "verdict": self.verdict,
            "quantity": self.quantity,
            "observed": self.observed,
            "tolerance": self.tolerance,
            "detail": self.detail,
            "reason": self.reason,
            "reference": self.reference,
        }


@dataclass
class _Registration:
    level: str
    name: str
    fn: Callable[[], Mapping[str, Any] | None]
    quantity: str = ""
    reference: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)


_CHECKS: dict[str, list[_Registration]] = {level: [] for level in LEVELS}


def check(
    level: str,
    name: str,
    *,
    quantity: str = "",
    reference: str = "",
    tags: tuple[str, ...] = (),
) -> Callable[[Callable[[], Mapping[str, Any] | None]], Callable[[], Mapping[str, Any] | None]]:
    """Register a zero-argument verification check at a level.

    ``quantity`` names what is being compared, which is what a report row shows
    next to the number.  ``reference`` names the source a V3 or V4 check is
    compared against.
    """
    if level not in _CHECKS:
        raise ValueError(f"unknown verification level {level!r}; expected {LEVELS}")

    def _decorate(fn: Callable[[], Mapping[str, Any] | None]):
        if any(reg.name == name for reg in _CHECKS[level]):
            raise ValueError(f"{level} check {name!r} is already registered")
        _CHECKS[level].append(
            _Registration(
                level=level,
                name=name,
                fn=fn,
                quantity=quantity,
                reference=reference,
                tags=tuple(tags),
            )
        )
        return fn

    return _decorate


def registered(level: str | None = None) -> dict[str, list[str]]:
    """Names of the registered checks, for ``tss verify`` and for tests."""
    if level is not None:
        return {level: [reg.name for reg in _CHECKS[level]]}
    return {lev: [reg.name for reg in regs] for lev, regs in _CHECKS.items()}


def _run_one(reg: _Registration) -> CheckResult:
    """Run one check and turn whatever happens into a record."""
    try:
        outcome = reg.fn() or {}
    except Skipped as exc:
        return CheckResult(
            level=reg.level,
            name=reg.name,
            verdict="skipped",
            quantity=reg.quantity,
            reason=exc.reason,
            reference=reg.reference,
        )
    except NotRun as exc:
        return CheckResult(
            level=reg.level,
            name=reg.name,
            verdict="not_run",
            quantity=reg.quantity,
            reason=exc.reason,
            reference=reg.reference,
        )
    except NotImplementedError as exc:
        # A check written against a stub reports the stub's own message, which
        # names the work package.  It starts passing or failing honestly the
        # moment the algorithm lands, with no change to the check.
        return CheckResult(
            level=reg.level,
            name=reg.name,
            verdict="not_run",
            quantity=reg.quantity,
            reason=str(exc),
            reference=reg.reference,
        )
    except Exception:  # noqa: BLE001 - a crashing check is a finding
        return CheckResult(
            level=reg.level,
            name=reg.name,
            verdict="fail",
            quantity=reg.quantity,
            detail=traceback.format_exc(limit=3).strip().splitlines()[-1],
            reason="the check raised an unexpected exception",
            reference=reg.reference,
        )

    verdict = str(outcome.get("verdict") or ("pass" if outcome.get("ok", True) else "fail"))
    if verdict not in VERDICTS:
        raise ValueError(f"{reg.level}/{reg.name} returned unknown verdict {verdict!r}")
    return CheckResult(
        level=reg.level,
        name=reg.name,
        verdict=verdict,
        quantity=str(outcome.get("quantity", reg.quantity)),
        observed=_maybe_float(outcome.get("observed")),
        tolerance=_maybe_float(outcome.get("tolerance")),
        detail=str(outcome.get("detail", "")),
        reason=str(outcome.get("reason", "")),
        reference=str(outcome.get("reference", reg.reference)),
    )


def _maybe_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def run_level(level: str) -> list[CheckResult]:
    """Run every check registered at one level, in registration order."""
    if level not in _CHECKS:
        raise ValueError(f"unknown verification level {level!r}; expected {LEVELS}")
    return [_run_one(reg) for reg in _CHECKS[level]]


def run_all() -> dict[str, list[CheckResult]]:
    """Run every level and return the results keyed by level."""
    return {level: run_level(level) for level in LEVELS}


# ===========================================================================
# V1  Internal consistency
# ===========================================================================
#
# These are the closed-form limiting cases docs/verification.md asks for.  Each
# one is a place a sign error or a factor of two would actually be caught; a
# check that only confirmed an array shape would catch almost nothing.


@check("V1", "beat_range_round_trip", quantity="relative range error, round trip")
def _v1_beat_range_round_trip() -> Mapping[str, Any]:
    """Range to beat frequency and back is exact.

    ``f_b = 2 K_f R / c`` and ``R = c f_b / (2 K_f)`` are exact inverses, so the
    round trip should be limited only by floating point.  This is the mapping
    every apparent range in the package goes through, so an error here is an
    error everywhere.
    """
    import numpy as np

    from tss.constants import beat_frequency
    from tss.l3_signal.base import range_from_beat_hz

    ranges = np.array([1.0e3, 5.0e4, 3.0e5, 1.0e6, 2.5e6])
    for k_f in (1.0e6, 2.8e6, 5.6e6, 1.12e7):
        back = range_from_beat_hz(beat_frequency(ranges, k_f), k_f)
        error = float(np.max(np.abs(back - ranges) / ranges))
        if error > 1.0e-12:
            return {"ok": False, "observed": error, "tolerance": 1.0e-12,
                    "detail": f"worst relative error at K_f = {k_f:.3g} Hz/s"}
    return {"observed": 0.0, "tolerance": 1.0e-12,
            "detail": "exact over four sweep rates and five ranges"}


@check("V1", "range_cell_noise_bandwidth", quantity="relative error against 1/T")
def _v1_range_cell_noise_bandwidth() -> Mapping[str, Any]:
    """The per-cell noise bandwidth equals the window ENBW over the block time.

    This is the identity the sub-watt sensitivity argument of C1 rests on: the
    noise admitted by one range cell is set by the coherent block duration and
    not by the swept bandwidth.  It is checked as a relation rather than
    asserted in prose.
    """

    from tss.constants import range_cell_noise_bandwidth
    from tss.l3_signal.transform import window_gains, window_taper

    worst = 0.0
    for k_f in (1.0e6, 5.6e6):
        for n, fs in ((4096, 2.0e5), (8192, 5.0e5)):
            duration = n / fs
            for window in ("rect", "hann", "blackman"):
                _, enbw = window_gains(window_taper(window, n))
                from tss.constants import C_LIGHT

                resolution = enbw * C_LIGHT / (2.0 * k_f * duration)
                bandwidth = float(range_cell_noise_bandwidth(resolution, k_f))
                expected = enbw / duration
                worst = max(worst, abs(bandwidth - expected) / expected)
    ok = worst < 1.0e-12
    return {"ok": ok, "observed": worst, "tolerance": 1.0e-12,
            "detail": "B_cell = enbw / T over three windows, two block lengths, two sweep rates"}


@check("V1", "sub_watt_cell_bandwidth_number", quantity="cell noise bandwidth, Hz")
def _v1_sub_watt_cell_bandwidth_number() -> Mapping[str, Any]:
    """A 1 km cell at K_f = 5.6 MHz/s admits about 37 Hz.

    The number quoted throughout the package and in the proposal.  It is
    checked because a quoted number that nothing verifies drifts.
    """
    from tss.constants import range_cell_noise_bandwidth

    observed = float(range_cell_noise_bandwidth(1.0e3, 5.6e6))
    ok = abs(observed - 37.36) < 0.05
    return {"ok": ok, "observed": observed, "tolerance": 0.05,
            "detail": "2 K_f dR / c for dR = 1 km, K_f = 5.6 MHz/s"}


@check("V1", "transform_tone_lands_in_expected_bin", quantity="bin index error")
def _v1_transform_tone_bin() -> Mapping[str, Any]:
    """A synthetic echo at a known range peaks in the predicted bin.

    A stationary point target with no noise and no window, de-chirped and
    transformed, must peak at the bin
    :func:`tss.l3_signal.transform.range_bin_index` predicts.  This couples the
    de-chirp sign convention, the orientation and the bin-to-range mapping in
    one test, which is where those three most often disagree.
    """
    import numpy as np

    from tss.constants import C_LIGHT
    from tss.interfaces import ReceivedSignal, SweepSchedule
    from tss.l3_signal.base import sweep_phase_rad
    from tss.l3_signal.dechirp import dechirp
    from tss.l3_signal.transform import range_bin_index, transform

    sweep = SweepSchedule(f_start_hz=1.0e6, f_stop_hz=5.0e6, k_f_hz_per_s=5.6e6,
                          sweep_period_s=1.0, direction=1)
    fs = 2.0e5
    n = 4096
    t = np.arange(n) / fs
    worst = 0.0
    for range_m in (1.0e5, 4.0e5, 9.0e5):
        tau = 2.0 * range_m / C_LIGHT
        samples = np.exp(1j * sweep_phase_rad(sweep, t - tau))[None, :]
        signal = ReceivedSignal(samples=samples, fs_hz=fs, sweep=sweep)
        spectrum = transform(dechirp(signal), window="rect", n_fft=n)
        peak = int(np.argmax(spectrum.power_db[0]))
        expected = float(range_bin_index(range_m, fs, n, 5.6e6))
        worst = max(worst, abs(peak - expected))
    ok = worst <= 0.5
    return {"ok": ok, "observed": worst, "tolerance": 0.5,
            "detail": "three ranges, rectangular window, no noise"}


@check("V1", "doppler_bias_is_antisymmetric_in_sweep_direction",
       quantity="half-sum residual, m")
def _v1_doppler_bias_antisymmetry() -> Mapping[str, Any]:
    """The up-sweep and down-sweep range biases are equal and opposite.

    This is the property the whole of C3 rests on, so it is checked against the
    forward model of :mod:`tss.l1_propagation.doppler` directly.  If the
    half-sum of the two apparent ranges is not the true range, the estimator of
    WP5 has nothing to recover.
    """

    from tss.l1_propagation.doppler import apparent_range_bias_m

    worst = 0.0
    for velocity in (-500.0, -50.0, 1.0, 250.0, 1000.0):
        for f_hz in (1.0e6, 1.0e7, 2.0e7):
            up = float(apparent_range_bias_m(velocity, f_hz, 5.6e6, direction=1))
            down = float(apparent_range_bias_m(velocity, f_hz, 5.6e6, direction=-1))
            worst = max(worst, abs(up + down))
    ok = worst < 1.0e-9
    return {"ok": ok, "observed": worst, "tolerance": 1.0e-9,
            "detail": "dd(up) + dd(down) = 0 over five velocities and three frequencies"}


@check("V1", "range_correlation_small_argument_limit",
       quantity="relative error against 4 (pi f_m tau)^2")
def _v1_range_correlation_limit() -> Mapping[str, Any]:
    """Range-correlation suppression goes as the square of the offset frequency.

    ``4 sin^2(pi f_m tau)`` tends to ``4 (pi f_m tau)^2`` for small
    ``f_m tau``, which is the statement that close-in oscillator noise is
    rejected hardest.  That statement is what relaxes the waveform purity
    requirement in Chapter 5 of the Technical Volume, so it is checked.
    """
    import numpy as np

    from tss.l3_signal.dechirp import range_correlation_suppression

    tau = 2.0e-3
    f_m = np.array([0.1, 1.0, 10.0])
    exact = range_correlation_suppression(f_m, tau)
    approx = 4.0 * (np.pi * f_m * tau) ** 2
    worst = float(np.max(np.abs(exact - approx) / approx))
    # The leading neglected term is (pi f_m tau)^2 / 3, which at the worst
    # point of this set, 10 Hz against a 2 ms delay, is 1.3e-3.  The tolerance
    # is set just above it rather than at a round number, so that the check
    # would notice an error of the same size as the term it is ignoring.
    ok = worst < 2.0e-3
    return {"ok": ok, "observed": worst, "tolerance": 2.0e-3,
            "detail": "tau = 2 ms, offsets 0.1 to 10 Hz; neglected term is "
                      "(pi f_m tau)^2 / 3 = 1.3e-3 at the worst point"}


@check("V1", "lamination_forward_inverse_round_trip", quantity="max thickness error, m")
def _v1_lamination_round_trip() -> Mapping[str, Any]:
    """The lamination inversion exactly inverts its own forward operator.

    ``R' = G t`` and the forward substitution that solves it are exact
    inverses, so a profile that really is piecewise linear in plasma frequency
    against height is recovered to machine precision.  Any error in this check
    is arithmetic, not physics.
    """
    import numpy as np

    from tss.l4_retrieval.lamination import (
        lamina_thicknesses_m,
        virtual_range_from_laminae_m,
    )

    freq = np.array([1.0, 1.4, 2.0, 3.0, 4.5, 6.0, 9.0]) * 1.0e6
    thickness = np.array([2.0e4, 3.0e4, 5.0e4, 4.0e4, 6.0e4, 2.5e4, 3.5e4])
    virtual = virtual_range_from_laminae_m(thickness, freq, 0.8e6)
    recovered = lamina_thicknesses_m(virtual, freq, 0.8e6)
    worst = float(np.max(np.abs(recovered - thickness)))
    ok = worst < 1.0e-6
    return {"ok": ok, "observed": worst, "tolerance": 1.0e-6,
            "detail": "seven laminae, exact forward operator"}


@check("V1", "lamination_recovers_exponential_topside", quantity="max height error, m")
def _v1_lamination_exponential() -> Mapping[str, Any]:
    """The inversion recovers a closed-form exponential topside.

    For ``f_N(h) = f_0 exp((h_s - h) / 2H)`` the group range has the closed
    form ``R'(f) = 2 H arccosh(f / f_0)`` and the reflection height is
    ``h_s - 2H ln(f / f_0)``.  This is a genuine limiting case: the profile is
    not piecewise linear in plasma frequency, so what is being tested is the
    lamina approximation and not the algebra.
    """
    import numpy as np

    from tss.l4_retrieval.lamination import lamina_thicknesses_m

    scale_h = 8.0e4
    alt_s = 5.0e5
    f0 = 5.0e5
    freq = f0 * np.exp(np.linspace(0.02, 2.0, 80))
    virtual = 2.0 * scale_h * np.arccosh(freq / f0)
    heights = alt_s - np.cumsum(lamina_thicknesses_m(virtual, freq, f0))
    truth = alt_s - 2.0 * scale_h * np.log(freq / f0)
    worst = float(np.max(np.abs(heights - truth)))
    ok = worst < 500.0
    return {"ok": ok, "observed": worst, "tolerance": 500.0,
            "detail": "80 laminae, H = 80 km, tolerance is one twentieth of an INS-05 cell"}


@check("V1", "assessment_reports_not_evaluated", quantity="unevaluated requirements")
def _v1_assessment_not_evaluated() -> Mapping[str, Any]:
    """An unmeasured requirement is reported ``not_evaluated``, never met.

    The check that protects the discipline of Section 3.1: a campaign that did
    not exercise a capability must not produce a verdict that reads as though
    it did.  It is run with an empty metric set, where every requirement should
    be unevaluated and none should be met.
    """
    from tss.l5_assessment.requirements_check import check_requirements

    records = check_requirements({})
    met = [rid for rid, r in records.items() if r["verdict"] in ("goal", "threshold")]
    unevaluated = [rid for rid, r in records.items() if r["verdict"] == "not_evaluated"]
    missing_reason = [rid for rid in unevaluated if not records[rid].get("reason")]
    ok = not met and len(unevaluated) == len(records) and not missing_reason
    return {"ok": ok, "observed": float(len(unevaluated)), "tolerance": float(len(records)),
            "detail": f"{len(met)} requirements scored as met from an empty metric set "
                      f"(must be 0); {len(missing_reason)} unevaluated rows carry no reason"}


@check("V1", "requirement_set_validates", quantity="problems found")
def _v1_requirement_set_validates() -> Mapping[str, Any]:
    """Every requirement carries a threshold, a sense and a traced source.

    :meth:`tss.requirements.RequirementSet.validate` refuses a requirement with
    no source, or with ``source: derived`` and no derivation named.  An
    untraceable number cannot be used to score a result, which is the first
    governing principle of Section 3.1.
    """
    from tss.requirements import load_requirements

    problems = load_requirements().validate()
    return {"ok": not problems, "observed": float(len(problems)), "tolerance": 0.0,
            "detail": "; ".join(problems[:3]) if problems else "requirement set is well formed"}


@check("V1", "l2_l3_snr_consistency", quantity="SNR discrepancy, dB")
def _v1_l2_l3_snr_consistency() -> Mapping[str, Any]:
    """The SNR L3 measures matches what the L2 link budget predicts.

    Not run.  This is the cross-layer check ``docs/verification.md`` calls one
    of the more productive V1 tests, because it couples a closed-form
    calculation to a sample-level simulation.  It needs the sample-level
    receiver of WP3 and a stated processing loss budget, and until both exist
    there is nothing to compare.
    """
    raise NotRun(
        "WP3 and WP4: the closed-form link budget and the sample-level receiver "
        "cannot yet be compared.  It needs the L2 baseline receiver to produce "
        "samples with a known noise density and a stated processing loss budget "
        "for the window and the detector.  See docs/verification.md section V1, "
        "'Conservation and consistency between layers'."
    )


@check("V1", "c3_estimator_recovers_truth_velocity", quantity="velocity error, m/s")
def _v1_c3_estimator() -> Mapping[str, Any]:
    """The up and down sweep estimator recovers the velocity that was put in.

    Not run until WP5 exists.  This is the definitive V1 test of C3: the
    forward model of :mod:`tss.l1_propagation.doppler` is complete, so the
    truth velocity is known, and the estimator either returns it or does not.
    The check is registered now so that it starts reporting the moment the
    estimator lands.
    """
    import numpy as np

    from tss.l3_signal.novel.sweep_pairing import range_rate_from_pair_mps

    # Raises NotImplementedError, which the harness turns into not_run with the
    # stub's own message naming WP5.
    range_rate_from_pair_mps(np.array([1.0]), np.array([1.0]), np.array([1.0e7]), 5.6e6)
    return {"ok": False, "detail": "unreachable while the estimator is a stub"}


@check("V1", "leakage_cancellation_meets_ins24", quantity="residual margin, dB")
def _v1_leakage_cancellation() -> Mapping[str, Any]:
    """The de-chirped leakage residual sits below the weakest admitted echo.

    Not run until WP5 exists.  Note that even when it does, this check bounds
    what processing can remove and not what the instrument achieves: the floor
    is set by front-end compression, quantisation and intermodulation, which
    are WP9 bench measurements.  C4 is a hardware claim with a modelling
    component.
    """
    from tss.l3_signal.novel.leakage_cancel import residual_floor_db

    residual_floor_db(None)  # raises NotImplementedError naming WP5
    return {"ok": False, "detail": "unreachable while the canceller is a stub"}


# ===========================================================================
# V2  Cross-fidelity
# ===========================================================================
#
# Two independent implementations of the same physics.  None of these can run
# yet, and the reasons differ: one needs a tier that is a stub, one needs
# MATLAB, one needs a model that is not written.  They are registered
# separately rather than as one "V2 blocked" row so that each unblocks on its
# own.


@check("V2", "analytic_vs_raytrace_group_delay", quantity="group delay difference, s",
       reference="PHaRLAP 4.7.4 against tss.l1_propagation.analytic")
def _v2_analytic_vs_raytrace() -> Mapping[str, Any]:
    """Analytic and ray-traced group delay agree at vertical incidence.

    Not run.  At vertical incidence in a horizontally uniform ionosphere the
    two must give the same group path, and any disagreement is a bug in one of
    them rather than physics.  It needs the PHaRLAP tier of WP2.
    """
    raise NotRun(
        "WP2: the ray-traced L1 tier is a stub and PHaRLAP runs MATLAB-side through "
        "the file bridge.  See docs/matlab-python-bridge.md and "
        "tss.l1_propagation.raytrace."
    )


@check("V2", "python_matlab_transcription", quantity="relative difference",
       reference="matlab/tss_verify_transcription.m")
def _v2_transcription() -> Mapping[str, Any]:
    """Ported functions match their MATLAB peers to machine precision.

    Not run.  ``docs/verification.md`` calls this the level that most often
    finds real bugs, because two implementations written at different times by
    the same author make different mistakes.  It needs MATLAB on the path.
    """
    raise Skipped(
        "MATLAB is not available in this environment, so the transcription check "
        "cannot run.  It is skipped rather than absent; see "
        "docs/matlab-python-bridge.md for the procedure and the tolerance standard."
    )


@check("V2", "iri2020_vs_nequick2", quantity="fractional density difference",
       reference="IRI-2020 against NeQuick 2")
def _v2_iri_vs_nequick() -> Mapping[str, Any]:
    """Two independent topside formulations over the same epochs.

    Not run.  This is a sensitivity test and not a correctness test: a result
    that changes materially between the two depends on the topside model rather
    than on the instrument and has to be reported that way.  It needs the
    NeQuick 2 tier, which is the outstanding half of WP1.
    """
    from tss.registry import available

    if "nequick2" not in available("l0")["l0"]:
        raise NotRun(
            "WP1: no l0 tier named 'nequick2' is registered, so there is no second "
            "topside formulation to compare against.  See docs/architecture.md "
            "section 'L0 Ionosphere state'."
        )
    raise NotRun("WP1: the comparison study is not written.")


@check("V2", "linkbudget_vs_sample_receiver", quantity="SNR difference, dB")
def _v2_linkbudget_vs_receiver() -> Mapping[str, Any]:
    """The closed-form link budget against the sample-level receiver.

    Not run.  The same comparison as the V1 cross-layer check, seen from the
    other side of the fidelity ladder.  It needs the WP3 sample-level receiver.
    """
    raise NotRun(
        "WP3: the sample-level receiver model does not yet produce a record with a "
        "stated noise density that the closed-form link budget can be compared "
        "against."
    )


# ===========================================================================
# V3  Literature replication
# ===========================================================================
#
# Each replication needs two things: digitised reference values under
# data/reference/, and a study that reproduces the published setup.  Neither
# exists for any of the five, so each reports not_run naming which piece is
# missing.  Several of these compare against digitised figures rather than
# tabulated numbers, so agreement to a few per cent is the practical ceiling
# and the tolerance should be set accordingly when they are written.


def _reference_missing(stem: str, what: str) -> NotRun:
    """Build the not_run for a replication whose reference data is absent."""
    from tss.config import repo_root

    path = repo_root() / "data" / "reference" / f"{stem}.json"
    return NotRun(
        f"WP8: no digitised reference values at {path}.  {what}  "
        "See docs/verification.md section V3 and docs/data-sources.md."
    )


@check("V3", "barry_1971_chirp_range_resolution", quantity="range resolution, m",
       reference="Barry 1971, HF chirp sounder")
def _v3_barry_1971() -> Mapping[str, Any]:
    """Chirp sounder technique and the range resolution relations of INS-05."""
    raise _reference_missing(
        "barry1971",
        "The replication would reproduce the de-chirp range resolution against "
        "sweep rate and integration time, which is the relation INS-05 rests on.",
    )


@check("V3", "poole_1985_ionogram_formation", quantity="virtual height, m",
       reference="Poole 1985, chirp ionogram formation")
def _v3_poole_1985() -> Mapping[str, Any]:
    """Chirp ionogram formation and the treatment of the sweep against the trace."""
    raise _reference_missing(
        "poole1985",
        "The replication would form an ionogram from a published sweep and trace "
        "configuration and compare the virtual heights.",
    )


@check("V3", "zawdie_2017_topside_ionogram", quantity="virtual height, m",
       reference="Zawdie 2017, modelled topside ionograms")
def _v3_zawdie_2017() -> Mapping[str, Any]:
    """Modelled topside ionograms, the closest published analogue to L1 to L3."""
    raise _reference_missing(
        "zawdie2017",
        "The replication would run L0 to L3 over the published ionosphere and "
        "compare the modelled topside ionogram.",
    )


@check("V3", "franklin_maclean_1969_sky_noise", quantity="noise density, dBm/Hz",
       reference="Franklin and MacLean 1969, in-orbit sky noise")
def _v3_franklin_maclean_1969() -> Mapping[str, Any]:
    """The external noise reference the link budget uses.

    This is the reference that runs out below about 0.6 MHz.  Below that the
    link budget extrapolates, which is the origin of the band-edge caveat on C1
    and the reason WP12 exists.  A replication that agreed above 0.6 MHz and
    said nothing below it would be reported that way, as a partial.
    """
    raise _reference_missing(
        "franklin_maclean1969",
        "The replication would compare the link budget's external noise model "
        "against the measured values above 0.6 MHz and record that below 0.6 MHz "
        "there is nothing to compare against.",
    )


@check("V3", "chartier_2022_propagation", quantity="group path difference, m",
       reference="Chartier 2022, PHaRLAP propagation comparison")
def _v3_chartier_2022() -> Mapping[str, Any]:
    """The PHaRLAP-based propagation comparison.

    The PHaRLAP installation used here comes from the author's Chartier
    verification work, so the ray-tracing configuration can be checked against
    a known-good setup rather than assembled from the manual.  That makes this
    the replication most likely to run first once WP2 lands.
    """
    raise NotRun(
        "WP2 and WP8: this replication needs the ray-traced L1 tier, which is a "
        "stub, and MATLAB with PHaRLAP on the path, which is not available here.  "
        "See tss.l1_propagation.raytrace and docs/verification.md section V3."
    )


# ===========================================================================
# V4  Archived real measurements
# ===========================================================================
#
# What V4 can and cannot reach is the most important limitation in
# docs/verification.md and it is repeated in the check docstrings rather than
# left in the document, because a reader of a verification report should not
# have to go and find it.


@check("V4", "isis1_profile_comparison", quantity="fractional density difference",
       reference="NASA SPDF hand-scaled ISIS-1 topside profiles")
def _v4_isis1_profiles() -> Mapping[str, Any]:
    """Retrieved profiles against hand-scaled ISIS-1 profiles.

    Not run without the cached archive.  This is the strongest evidence
    available for the ionospheric half of the chain: it exercises L0, L1, L4
    and the profile-accuracy part of L5 against measurements made by a real
    topside sounder in a real ionosphere.
    """
    from tss.config import repo_root

    cache = repo_root() / "data" / "external"
    candidates = sorted(cache.glob("**/isis1*")) if cache.exists() else []
    if not candidates:
        raise Skipped(
            f"no cached ISIS-1 archive under {cache}; run 'tss fetch spdf_isis' "
            "first.  The test-suite portion runs against tests/fixtures/ and never "
            "touches the network."
        )
    raise NotRun(
        "WP8: the archive comparison study (tss.studies.isis_validation) is not "
        "written, so the cached archive cannot yet be driven through the pipeline."
    )


@check("V4", "fmcw_half_not_reachable_from_archive", quantity="")
def _v4_fmcw_not_reachable() -> Mapping[str, Any]:
    """The FMCW half of the chain cannot be validated against the archive.

    Registered as a check so that the limitation appears in every verification
    report rather than only in ``docs/verification.md``.  Alouette and ISIS were
    pulse sounders: they transmitted a short pulse and switched the receiver on
    afterwards.  So the archive contains no range-Doppler coupling, which is
    C3; no continuous-wave transmit-to-receive leakage, which is C4; and no
    sweep nonlinearity, coherent processing interval or oscillator stability
    term.  Those are validated by V1 to V3 and, for the numbers C4 turns on, on
    the bench in WP9.
    """
    raise NotRun(
        "By construction.  Alouette and ISIS were pulse sounders, so the archive "
        "contains no range-Doppler coupling (C3), no continuous-wave leakage (C4) "
        "and no FMCW-specific waveform terms.  This is a permanent limitation of "
        "V4 and not a gap that will close; see docs/verification.md section V4, "
        "'What it cannot reach'."
    )
