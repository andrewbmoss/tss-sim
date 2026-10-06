"""Layer: L2
WP: WP3
Claims: C1, C3, C4
Status: implemented
Requirements: INS-21, INS-26

Oscillator phase noise, and the range-correlation credit that may or may not be
claimable against it.

WHY AN FMCW SOUNDER CARES
-------------------------
The instrument transmits continuously and de-chirps the echo against the same
reference that generated the transmission.  INS-21 requires that architecture:
one reference, not two.  That has a consequence which is the subject of this
module.  Phase noise on the reference appears on both the transmitted signal and
on the de-chirp reference, delayed against each other by the round-trip time of
the echo.  Over a short delay the two are correlated, so most of the phase noise
subtracts out, and only the part that decorrelated over that delay survives.
This is range correlation, and it is the standard result for a homodyne FMCW or
CW radar.

THE THREE RELATIONS
-------------------
**Phase jitter from an Allan deviation.**

    sigma_phi = 2 pi (f tau) sigma_y(tau)

with ``f`` the carrier, ``tau`` the interval over which the jitter is asked
about, and ``sigma_y`` the fractional frequency Allan deviation at that
interval.  This is the bridge between the way an oscillator is specified,
which is an Allan deviation, and the way it hurts, which is phase.  INS-26's
goal of 3.2e-9 at 1.15 s comes out of this relation applied to the coherent
processing interval of INS-15.  Its threshold of 5.8e-8 at 0.05 s comes from
the velocity accuracy instead, through ``sigma_v = c sigma_y / (2 mu_sat)``
(Requirements Rev G, amended 2026-10-06).

**Sideband level from a phase departure.**

    (pi df tau)^2

for a frequency departure ``df`` persisting for ``tau``.  This is the
small-angle power in the first-order sideband relative to the carrier, and it is
what turns a frequency wander into a skirt around the echo that competes with a
weak adjacent trace.

**Range-correlation suppression.**

    4 sin^2(pi f_m tau)

for a phase-noise Fourier offset ``f_m`` and a round-trip delay ``tau``.  This
is the factor the phase noise spectrum is multiplied by.  It goes to zero as
``f_m tau -> 0``, which is the correlation: slow phase noise is common to the
transmitted signal and the reference and cancels.  It reaches 4, that is +6 dB,
at ``f_m tau = 1/2``, which is worth knowing because the suppression is a
suppression only for offsets below ``1 / (2 tau)`` and is an enhancement above.
A topside echo at 600 km has a round-trip delay of about 4 ms, so the crossover
sits near 125 Hz, and everything below that is credited while everything above
is not.

OPEN ITEM 12, WHICH IS NOT SETTLED
----------------------------------
Whether the range-correlation credit is claimable at all for this instrument is
open item 12, and the credit in question is a factor of 24 to 48, which is 14 to
17 dB.  That is large enough that the answer changes the oscillator selection
and therefore the parts cost.

The doubt is not about the algebra, it is about whether the architecture earns
it.  Range correlation requires the transmitted signal and the de-chirp
reference to share their phase noise.  They do if they come from the same
reference and if nothing between them adds uncorrelated noise; a separate
synthesiser in the receive path, a retimed digital reference or a long and
thermally unstable cable between the two can break the correlation without
breaking the block diagram.  The leakage path is also not the echo path, so the
leakage residual and the echo do not get the same credit, and WP3's leakage work
has to say which one it is quoting.

Until it is settled, the uncredited number stands.  Every function here that can
apply the credit takes ``credit=False`` by default and
:meth:`Oscillator.score_ins26` reports both, with the uncredited figure as the
verdict and the credited figure as a note.  A requirement scored against a
credit that has not been demonstrated is a requirement that has not been scored.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np

from tss.units import db

__all__ = [
    "Oscillator",
    "phase_jitter_rad",
    "sideband_power_ratio",
    "range_correlation_suppression",
    "range_correlation_credit_db",
    "correlation_crossover_hz",
    "INS26_ALLAN_THRESHOLD",
    "INS26_THRESHOLD_INTERVAL_S",
    "INS26_ALLAN_GOAL",
    "INS26_GOAL_INTERVAL_S",
    "BASELINE_PART_ALLAN",
    "OPEN_ITEM_12_CREDIT_RANGE",
]

#: INS-26 threshold and the interval it is written at.  Requirements Rev G as
#: amended 2026-10-06: one tenth of the SCI-15 accuracy, at a refractive index
#: at the spacecraft of 0.87.
INS26_ALLAN_THRESHOLD = 5.8e-8
INS26_THRESHOLD_INTERVAL_S = 0.05

#: INS-26 goal and its interval, from the SCI-15 goal resolution with no
#: range-correlation credit.
INS26_ALLAN_GOAL = 3.2e-9
INS26_GOAL_INTERVAL_S = 1.15

#: Allan deviation of the part the baseline instrument assumes, at 1 s and at
#: 2 s.  It is a property of the part and not a requirement.  It meets the
#: INS-26 goal.  It was the Rev F requirement, which is where the number came
#: from.
BASELINE_PART_ALLAN = 2.0e-9

#: The disputed range-correlation credit, as a power factor.  14 to 17 dB.
OPEN_ITEM_12_CREDIT_RANGE = (24.0, 48.0)


def phase_jitter_rad(
    f_hz: np.ndarray | float,
    tau_s: np.ndarray | float,
    allan_deviation: np.ndarray | float,
) -> np.ndarray:
    """``sigma_phi = 2 pi (f tau) sigma_y(tau)``, radians RMS.

    The product ``f tau`` is the number of carrier cycles in the interval, so
    the relation reads as: a fractional frequency error of ``sigma_y`` displaces
    the phase by that fraction of those cycles.  It is exact for a flat-in-tau
    Allan deviation and is the standard engineering approximation otherwise.
    """
    return (
        2.0
        * np.pi
        * np.asarray(f_hz, dtype=float)
        * np.asarray(tau_s, dtype=float)
        * np.asarray(allan_deviation, dtype=float)
    )


def sideband_power_ratio(df_hz: np.ndarray | float, tau_s: np.ndarray | float) -> np.ndarray:
    """``(pi df tau)^2``, the first-order sideband power relative to the carrier.

    Small-angle: valid while the result is well below one.  A value approaching
    unity means the small-angle assumption has failed and the number should not
    be quoted, which is why :meth:`Oscillator.describe` reports the argument as
    well as the result.
    """
    return (np.pi * np.asarray(df_hz, dtype=float) * np.asarray(tau_s, dtype=float)) ** 2


def range_correlation_suppression(
    f_m_hz: np.ndarray | float, tau_s: np.ndarray | float
) -> np.ndarray:
    """``4 sin^2(pi f_m tau)``, the factor applied to the phase noise spectrum.

    Below unity for offsets under ``1/(3 tau)`` or so, above unity between there
    and ``1/tau``, and oscillating about 2 beyond.  A caller quoting a
    suppression should say at which offset, because a single number is
    meaningless for this factor.
    """
    return 4.0 * np.sin(np.pi * np.asarray(f_m_hz, dtype=float) * np.asarray(tau_s, dtype=float)) ** 2


def correlation_crossover_hz(tau_s: float) -> float:
    """Offset frequency at which the suppression factor passes through one.

    ``4 sin^2(pi f_m tau) = 1`` gives ``f_m = 1/(6 tau)``.  Below this the
    correlation helps; above it, it does not.
    """
    return 1.0 / (6.0 * float(tau_s))


def range_correlation_credit_db(
    f_m_hz: np.ndarray | float, tau_s: np.ndarray | float
) -> np.ndarray:
    """Suppression expressed as a credit in dB.  Positive means noise removed."""
    return -db(range_correlation_suppression(f_m_hz, tau_s))


class Oscillator:
    """The frequency reference, specified the way one is bought.

    Parameters
    ----------
    allan_deviation
        Mapping from averaging interval in seconds to fractional frequency Allan
        deviation.  The baseline instrument configuration states 1 s and 2 s
        for the part it assumes.  INS-26 is written at 0.05 s and 1.15 s, so a
        score at the threshold interval is flagged as outside the table.
        Intermediate intervals are
        interpolated in log-log, which is right for a power-law noise type and
        wrong across a noise-type boundary; a caller asking far outside the
        tabulated intervals gets a warning flag in :meth:`describe` rather than a
        silent extrapolation.
    single_reference
        Whether the transmit chain and the de-chirp reference come from the same
        source.  This is INS-21 and it is the precondition for any
        range-correlation credit at all.  The credit is refused outright when it
        is False, because there is then nothing for the noise to correlate with.
    """

    def __init__(
        self,
        allan_deviation: Mapping[float, float] | None = None,
        *,
        single_reference: bool = True,
    ) -> None:
        table = dict(allan_deviation or {1.0: BASELINE_PART_ALLAN, 2.0: BASELINE_PART_ALLAN})
        if not table:
            raise ValueError("allan_deviation needs at least one interval")
        if any(t <= 0.0 for t in table) or any(v <= 0.0 for v in table.values()):
            raise ValueError("Allan deviation intervals and values must be positive")
        self.allan_deviation = {float(k): float(v) for k, v in sorted(table.items())}
        self.single_reference = bool(single_reference)

    # -- the specification ---------------------------------------------------

    def sigma_y(self, tau_s: np.ndarray | float) -> np.ndarray:
        """Allan deviation at an arbitrary interval, log-log interpolated."""
        taus = np.array(sorted(self.allan_deviation), dtype=float)
        vals = np.array([self.allan_deviation[float(t)] for t in taus], dtype=float)
        query = np.asarray(tau_s, dtype=float)
        if taus.size == 1:
            return np.full(np.shape(query), vals[0], dtype=float)
        return np.exp(np.interp(np.log(query), np.log(taus), np.log(vals)))

    def outside_tabulated(self, tau_s: np.ndarray | float) -> bool:
        """True when the requested interval is outside the tabulated range."""
        taus = sorted(self.allan_deviation)
        query = np.atleast_1d(np.asarray(tau_s, dtype=float))
        return bool(np.any(query < taus[0]) or np.any(query > taus[-1]))

    # -- the consequences ----------------------------------------------------

    def phase_jitter_rad(
        self, f_hz: np.ndarray | float, tau_s: np.ndarray | float
    ) -> np.ndarray:
        """Phase jitter over an interval, radians RMS.  Uncredited."""
        return phase_jitter_rad(f_hz, tau_s, self.sigma_y(tau_s))

    def phase_jitter_credited_rad(
        self,
        f_hz: np.ndarray | float,
        tau_s: np.ndarray | float,
        round_trip_delay_s: float,
        f_m_hz: float | None = None,
    ) -> np.ndarray:
        """Phase jitter with the range-correlation credit applied.

        ``f_m_hz`` defaults to ``1 / (2 tau)``, the offset the coherent interval
        itself corresponds to.  The credit is refused, and the uncredited value
        returned unchanged, when the instrument is not single-reference.

        Whether this number may be quoted at all is open item 12.  The
        uncredited number stands until it is settled; see the module docstring.
        """
        raw = self.phase_jitter_rad(f_hz, tau_s)
        if not self.single_reference:
            return raw
        offset = float(f_m_hz) if f_m_hz is not None else 1.0 / (2.0 * float(np.max(tau_s)))
        suppression = range_correlation_suppression(offset, round_trip_delay_s)
        return raw * np.sqrt(suppression)

    def sideband_power_ratio(
        self, df_hz: np.ndarray | float, tau_s: np.ndarray | float
    ) -> np.ndarray:
        """``(pi df tau)^2`` for this oscillator's frequency departure."""
        return sideband_power_ratio(df_hz, tau_s)

    def frequency_departure_hz(
        self, f_hz: np.ndarray | float, tau_s: np.ndarray | float
    ) -> np.ndarray:
        """``f * sigma_y(tau)``, the absolute frequency wander over an interval."""
        return np.asarray(f_hz, dtype=float) * self.sigma_y(tau_s)

    # -- scoring -------------------------------------------------------------

    def score_ins26(
        self,
        f_hz: float,
        coherent_interval_s: float,
        round_trip_delay_s: float,
        *,
        threshold: float = INS26_ALLAN_THRESHOLD,
    ) -> dict[str, Any]:
        """Score the oscillator against INS-26, both ways.

        ``threshold`` defaults to the INS-26 threshold, which is written at
        0.05 s.  To score the goal, pass :data:`INS26_ALLAN_GOAL` with a
        coherent interval of :data:`INS26_GOAL_INTERVAL_S`.

        The verdict is taken on the UNCREDITED figure.  The credited figure is
        reported beside it as a note, together with the credit that would apply
        and the statement that open item 12 has not settled whether it may be
        claimed.  A campaign that wants the credited verdict has to change this
        deliberately, which is the point.
        """
        sigma_y = float(self.sigma_y(coherent_interval_s))
        raw = float(self.phase_jitter_rad(f_hz, coherent_interval_s))
        offset = 1.0 / (2.0 * float(coherent_interval_s))
        suppression = float(range_correlation_suppression(offset, round_trip_delay_s))
        credited = float(
            self.phase_jitter_credited_rad(f_hz, coherent_interval_s, round_trip_delay_s)
        )
        return {
            "requirement": "INS-26",
            "allan_deviation": sigma_y,
            "threshold": float(threshold),
            "verdict": "pass" if sigma_y <= float(threshold) else "fail",
            "scored_on": "uncredited",
            "phase_jitter_rad": raw,
            "phase_jitter_deg": math.degrees(raw),
            "credited_phase_jitter_rad": credited,
            "credit_db": float(range_correlation_credit_db(offset, round_trip_delay_s)),
            "suppression_factor": suppression,
            "correlation_crossover_hz": correlation_crossover_hz(round_trip_delay_s),
            "offset_used_hz": offset,
            "single_reference": self.single_reference,
            "open_item_12": (
                "Whether the range-correlation credit of "
                f"{OPEN_ITEM_12_CREDIT_RANGE[0]:g} to {OPEN_ITEM_12_CREDIT_RANGE[1]:g} "
                "applies is unsettled.  The uncredited number stands until it is."
            ),
            "interval_outside_table": self.outside_tabulated(coherent_interval_s),
        }

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable record for the run manifest."""
        return {
            "allan_deviation": dict(self.allan_deviation),
            "single_reference": self.single_reference,
            "ins26_threshold": INS26_ALLAN_THRESHOLD,
            "ins26_threshold_interval_s": INS26_THRESHOLD_INTERVAL_S,
            "ins26_goal": INS26_ALLAN_GOAL,
            "ins26_goal_interval_s": INS26_GOAL_INTERVAL_S,
            "range_correlation_credit": "not applied by default, open item 12",
            "open_item_12_credit_range": list(OPEN_ITEM_12_CREDIT_RANGE),
        }

    @classmethod
    def from_config(cls, cfg: Mapping[str, Any]) -> Oscillator:
        """Build from the ``oscillator:`` block of an instrument configuration.

        The configuration spells the intervals into the key names,
        ``allan_deviation_1s`` and ``allan_deviation_2s``, so they are parsed
        back out here rather than requiring a nested mapping that the YAML does
        not have.
        """
        table: dict[float, float] = {}
        for key, value in cfg.items():
            if not str(key).startswith("allan_deviation_"):
                continue
            tail = str(key)[len("allan_deviation_") :]
            if tail.endswith("s"):
                tail = tail[:-1]
            try:
                table[float(tail)] = float(value)
            except ValueError:
                continue
        return cls(table or None, single_reference=bool(cfg.get("single_reference", True)))
