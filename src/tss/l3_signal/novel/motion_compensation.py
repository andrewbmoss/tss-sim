"""Layer: L3
WP: WP5
Claims: C3
Status: stub
Requirements: INS-15, SCI-23

Motion-compensated coherent integration over the known ephemeris.  Not
implemented.  WP5.

What it is for
--------------
Pairing an up-sweep with a down-sweep, which is
:mod:`tss.l3_signal.novel.sweep_pairing`, removes the first-order range bias.
It does not by itself deliver a velocity resolution.  Resolving one velocity
from another needs a coherent interval long enough for their Doppler shifts to
separate, and that is what this module is for: integrating several sweeps
coherently, after compensating for the phase history the platform's own motion
imposes, so that the residual Doppler is the ionosphere's and not the
spacecraft's.

The interval the requirement demands
-------------------------------------
Two returns whose line-of-sight velocities differ by ``delta_v`` have Doppler
shifts differing by ``2 delta_v f / c``.  Resolving them needs a coherent
interval of at least the reciprocal of that difference:

    T >= c / (2 * delta_v * f)

The consequence worth stating is that the product is frequency independent:

    f * T = c / (2 * delta_v)

so the demand is on the velocity resolution alone and not on where in the band
the sounding sits.  For the 50 m/s threshold of SCI-23 the product is about
3.0e6, which is 3.0 s at 1 MHz, 0.30 s at 10 MHz and 0.15 s at 20 MHz.  For the
25 m/s goal it is twice that.  The relation is implemented in
:func:`tss.l1_propagation.doppler.coherent_interval_for_velocity_resolution_s`
and INS-15 takes its 1 s threshold and 2 s goal from it, evaluated at the low
end of the band where the demand is hardest.  That is the free-space relation.
Inside the plasma the interval lengthens by ``1 / mu_sat``, the reciprocal of
the refractive index at the platform, which the function takes as a keyword.

Why that interval is the difficulty rather than the specification
-------------------------------------------------------------------
One to two seconds is a long time for this platform.  At 7.6 km/s the
spacecraft moves 7.6 to 15 km along track within a single velocity estimate.
That is comparable to the along-track sampling SCI-09 asks for, 50 km at
threshold and 25 km at goal, so the aperture over which the estimate is formed
is a substantial fraction of the cell the estimate is supposed to describe.

The geometry therefore changes materially inside one coherent interval, and the
platform cannot be treated as quasi-stationary.  Three consequences follow, and
they are the reasons this is research rather than a standard integration:

* The line-of-sight direction to the reflection point rotates during the
  interval, so the projection that produces the Doppler shift is not constant
  and the phase history is not a straight line in time.
* The reflection point itself moves along track with the spacecraft, so the
  signal being integrated is not a single scatterer observed repeatedly but a
  succession of nearby reflection points.  Whether they stay coherent over the
  interval is an open question about the ionosphere, not about the processor,
  and it is the same question the aperture synthesis of
  :mod:`tss.l3_signal.novel.aperture_synthesis` asks from the other side.
* The range itself changes during the interval by more than a range cell in
  some geometries, so the integration has to follow the echo across cells
  rather than integrating one cell.

The compensation is nonetheless tractable in principle, because the platform
motion is known.  The ephemeris gives the position and velocity to far better
accuracy than the phase needs, so the phase history can be predicted and
removed rather than estimated.  What is not known is the reflection point's
motion, which is the signal, and the ionosphere's own decorrelation time, which
bounds how much of the interval is usable.

Oscillator stability enters here
---------------------------------
Coherent integration over 1 to 2 s requires the reference to stay coherent over
1 to 2 s.  INS-26 takes its Allan deviation threshold from
``sigma_phi = 2 pi (f T) sigma_y(T)`` evaluated at this interval, so the
coherent interval chosen here sets an oscillator requirement directly.  Whether
the range-correlation credit described in :mod:`tss.l3_signal.dechirp` applies
to it is open item 12 in the proposal and it would relax INS-26 by more than an
order of magnitude if it does.  It is not assumed here.

Open questions
--------------
1.  Over what interval the topside echo stays coherent.  Nothing in the
    archived pulse-sounder data bears on it, because a pulse sounder never
    integrated coherently across sweeps.
2.  Whether the along-track motion during the interval degrades the velocity
    estimate or improves the along-track resolution.  Both are true at once and
    which dominates is the trade this module and
    :mod:`tss.l3_signal.novel.aperture_synthesis` have to settle between them.
3.  How the compensation interacts with the up and down sweep pairing.  The
    pairing wants the two sweeps close together in time; the integration wants
    many sweeps.  Those pull in opposite directions and the resolution of that
    conflict is the design of the sweep schedule, which is why
    :class:`~tss.interfaces.SweepSchedule` carries ``direction`` and a campaign
    supplies a sequence of schedules rather than one.

Raises
------
NotImplementedError
    On every call.  WP5.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from tss.interfaces import PlatformState, ReceivedSignal, SweepSchedule

__all__ = [
    "required_coherent_interval_s",
    "phase_history_rad",
    "compensate_motion",
    "integrate_coherently",
]

_MESSAGE = (
    "WP5: motion-compensated coherent integration is not implemented.  The "
    "required interval follows from T >= c / (2 delta_v f), which is implemented "
    "in tss.l1_propagation.doppler.coherent_interval_for_velocity_resolution_s, "
    "but the compensation itself is unstarted: over the 1 to 2 s interval INS-15 "
    "demands, the spacecraft moves 7.6 to 15 km and the geometry changes "
    "materially within one velocity estimate, so the platform cannot be treated "
    "as quasi-stationary.  See docs/work-packages.md row WP5.  Claim C3."
)


def required_coherent_interval_s(
    delta_v_mps: float | np.ndarray,
    f_hz: float | np.ndarray,
) -> np.ndarray:
    """Coherent interval needed for a velocity resolution.  Not implemented here.

    The relation ``T >= c / (2 delta_v f)`` is settled arithmetic and it is
    implemented once, in
    :func:`tss.l1_propagation.doppler.coherent_interval_for_velocity_resolution_s`.
    This function exists only so that a reader looking for it in the module
    that depends on it is sent to the implementation rather than finding a
    second copy that can drift from the first.

    Raises
    ------
    NotImplementedError
        Always, with a pointer to the implementation.  WP5.
    """
    del delta_v_mps, f_hz
    raise NotImplementedError(
        "WP5: use tss.l1_propagation.doppler.coherent_interval_for_velocity_resolution_s, "
        "which implements T >= c / (2 delta_v f).  It is not duplicated here so that "
        "the two cannot drift apart.  " + _MESSAGE
    )


def phase_history_rad(
    platform: Sequence[PlatformState],
    f_hz: float,
    u_los: np.ndarray,
) -> np.ndarray:
    """Predicted phase against time from the platform motion.  Not implemented.

    Intended to integrate ``2 pi f_d(t)`` over the coherent interval using the
    ephemeris, with the line of sight recomputed at each step rather than held
    fixed, because holding it fixed is exactly the quasi-stationary assumption
    that does not hold over 1 to 2 s.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del platform, f_hz, u_los
    raise NotImplementedError(_MESSAGE)


def compensate_motion(
    signal: ReceivedSignal,
    platform: Sequence[PlatformState],
    sweep: SweepSchedule,
) -> ReceivedSignal:
    """Remove the platform's own phase history from a record.  Not implemented.

    Intended to return a new :class:`~tss.interfaces.ReceivedSignal` whose
    residual Doppler is the ionosphere's rather than the spacecraft's, with the
    range migration across cells followed rather than ignored.  The input is
    not modified.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del signal, platform, sweep
    raise NotImplementedError(_MESSAGE)


def integrate_coherently(
    signals: Sequence[ReceivedSignal],
    platform: Sequence[PlatformState],
    *,
    interval_s: float = 1.0,
) -> ReceivedSignal:
    """Coherently integrate several sweeps after compensation.  Not implemented.

    ``interval_s`` defaults to the INS-15 threshold of 1 s.  How much of the
    interval is usable depends on the coherence time of the topside echo, which
    is open question 1 in the module docstring and has no measured answer.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del signals, platform, interval_s
    raise NotImplementedError(_MESSAGE)
