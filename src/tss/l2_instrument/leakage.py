"""Layer: L2
WP: WP3
Claims: C4
Status: stub
Requirements: INS-23, INS-24

Transmit-to-receive leakage and its de-chirped residual.  Not implemented.
WP3, and it is on the critical path.

WHY THIS IS THE HARDEST PART OF THE INSTRUMENT
----------------------------------------------
The sounder transmits continuously while it receives.  The transmitter is
therefore the strongest signal at the receiver input at all times, by a very
large margin: the echo arrives after free-space spreading over twice the range
and after ionospheric absorption, while the leakage takes a path of a metre or
two.  INS-23 sets the isolation threshold at 40 dB, which is not enough on its
own, and the rest has to come from the de-chirp.

The de-chirp is what makes it survivable.  Leakage at a delay of a few
nanoseconds de-chirps to a beat frequency of a few tens of millihertz, which is
essentially at zero, while an echo at 600 km de-chirps to tens of kilohertz.  So
the leakage and the echo separate in the transform even though they do not
separate in time.  INS-24 is a requirement on what is left after that
separation, which is the residual, and the residual is what this module has to
produce.

WHAT MAKES THE RESIDUAL NON-TRIVIAL
-----------------------------------
If the leakage were a clean scaled copy of the transmitted sweep at a fixed
delay, its de-chirped residual would be a single tone at near zero frequency and
it would be removed by a high-pass filter.  Three things stop that.

1.  **The isolation is not flat with frequency.**  A coupling path through a
    deployed structure has resonances, and the sweep crosses them, so the
    leakage amplitude is modulated as the sweep proceeds.  Amplitude modulation
    of a tone at zero beat frequency puts sidebands at the modulation rate,
    which is not at zero, and those land in the range cells where echoes live.
2.  **The coupling delay is not fixed.**  A deployed antenna moves.  The
    coupling delay of the baseline configuration is stated as 3 ns and carries
    open item 13, because the deployed geometry is not fixed yet.  A slowly
    varying delay is a slowly varying phase on the residual, and the rate of
    variation sets how well an adaptive canceller can track it.
3.  **The leakage is strong enough to matter through the receiver's
    nonlinearity.**  It is the signal that decides whether the front end clips,
    and a clipped leakage produces intermodulation that is no longer at zero
    beat frequency at all.  That is the coupling between this module and
    :mod:`tss.l2_instrument.receiver` and it is why the two cannot be developed
    independently.

WHAT HAS TO BE DECIDED FIRST
----------------------------
The frequency-dependent isolation profile is the input this module does not
have.  A flat 40 dB is the requirement threshold, not a measurement, and the
whole question is what the departures from flat do.  Until WP3's antenna and
structure work produces a coupling model, or a measurement from the prototype
gives one, any residual this module computed would be a residual for an
isolation profile that was invented to produce it.

The companion at L3 is ``tss.l3_signal.leakage_suppress``, which is also not
started.  The two score different things: INS-23 is the isolation this module
models, INS-24 is the residual after both this module and that one have had
their turn.

Raises
------
NotImplementedError
    On every call.  WP3, ``docs/work-packages.md`` row WP3 and
    ``docs/architecture.md``, the L2 tier table row ``leakage``.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from tss.interfaces import ReceivedSignal, SweepSchedule

__all__ = [
    "LeakageModel",
    "isolation_profile_db",
    "leakage_signal",
    "dechirped_residual",
    "score_ins23_ins24",
    "BASELINE_ISOLATION_DB",
    "BASELINE_COUPLING_DELAY_S",
]

#: INS-23 threshold, dB.  A requirement, not a measurement.
BASELINE_ISOLATION_DB = 40.0

#: Coupling delay of the baseline configuration, seconds.  Open item 13: the
#: deployed geometry is not fixed, so this is a placeholder with the right order
#: of magnitude and not a number any result may rest on.
BASELINE_COUPLING_DELAY_S = 3.0e-9

_MESSAGE = (
    "WP3: transmit-to-receive leakage is not implemented.  It is on the critical "
    "path for C4.  The blocking input is a frequency-dependent isolation profile; "
    "a flat 40 dB is the INS-23 threshold rather than a coupling model, and the "
    "whole question is what the departures from flat do to the de-chirped "
    "residual.  See docs/work-packages.md row WP3."
)


class LeakageModel:
    """The coupling path from transmitter to receiver.  Not implemented.

    Parameters
    ----------
    isolation_db
        Isolation at the reference frequency, positive decibels.  Default
        :data:`BASELINE_ISOLATION_DB`, the INS-23 threshold.
    coupling_delay_s
        One-way delay of the coupling path.  Default
        :data:`BASELINE_COUPLING_DELAY_S`, which carries open item 13.
    isolation_ripple_db
        Peak-to-peak departure of the isolation from flat across the band.  This
        is the parameter the whole module turns on and it has no defensible
        default, so it has none.
    ripple_period_hz
        Frequency period of that ripple, which for a structural resonance is set
        by the path length.  Together with ``isolation_ripple_db`` it sets where
        the residual sidebands land.
    delay_jitter_s
        Variation of the coupling delay with time, from a moving deployed
        structure.  Sets how fast an adaptive canceller has to track.

    Raises
    ------
    NotImplementedError
        Always.  WP3.
    """

    def __init__(
        self,
        isolation_db: float = BASELINE_ISOLATION_DB,
        coupling_delay_s: float = BASELINE_COUPLING_DELAY_S,
        isolation_ripple_db: float | None = None,
        ripple_period_hz: float | None = None,
        delay_jitter_s: float = 0.0,
    ) -> None:
        del (
            isolation_db,
            coupling_delay_s,
            isolation_ripple_db,
            ripple_period_hz,
            delay_jitter_s,
        )
        raise NotImplementedError(_MESSAGE)


def isolation_profile_db(
    f_hz: np.ndarray,
    model: LeakageModel,
) -> np.ndarray:
    """Isolation against sounding frequency, positive decibels.  Not implemented.

    The intended return is an array parallel to ``f_hz``.  This is the quantity
    INS-23 is scored on, and the departures from flat in it are what produce the
    residual INS-24 is scored on, so the two requirements share this one array.

    Raises
    ------
    NotImplementedError
        Always.  WP3.
    """
    del f_hz, model
    raise NotImplementedError(_MESSAGE)


def leakage_signal(
    sweep: SweepSchedule,
    fs_hz: float,
    model: LeakageModel,
    rng: np.random.Generator,
    *,
    n_samples: int | None = None,
    transmit_power_w: float = 1.0,
) -> np.ndarray:
    """The leakage at the receiver input, as complex samples.  Not implemented.

    Intended to be added to the echo samples before the receiver is applied, not
    after, because whether the front end clips is decided by the sum and the
    clipping is the coupling to :mod:`tss.l2_instrument.receiver`.

    Randomness comes from ``rng`` and from nowhere else, per contract rule 4.

    Raises
    ------
    NotImplementedError
        Always.  WP3.
    """
    del sweep, fs_hz, model, rng, n_samples, transmit_power_w
    raise NotImplementedError(_MESSAGE)


def dechirped_residual(
    signal: ReceivedSignal,
    model: LeakageModel,
    *,
    suppression_db: float = 0.0,
) -> np.ndarray:
    """What is left of the leakage after the de-chirp.  Not implemented.

    ``suppression_db`` is the credit taken for whatever
    ``tss.l3_signal.leakage_suppress`` achieved, so that this function can be
    used to score INS-24 both with and without the L3 canceller and the two can
    be reported separately.  Neither exists yet.

    Raises
    ------
    NotImplementedError
        Always.  WP3.
    """
    del signal, model, suppression_db
    raise NotImplementedError(_MESSAGE)


def score_ins23_ins24(
    f_hz: np.ndarray,
    model: LeakageModel,
    signal: ReceivedSignal | None = None,
) -> Mapping[str, Any]:
    """Score the isolation and the residual.  Not implemented.

    INS-23 scores the isolation profile, INS-24 the de-chirped residual relative
    to the weakest echo that has to be detected beside it.  They are scored
    together because a design can trade one against the other and reporting
    either alone hides the trade.

    Raises
    ------
    NotImplementedError
        Always.  WP3.
    """
    del f_hz, model, signal
    raise NotImplementedError(_MESSAGE)
