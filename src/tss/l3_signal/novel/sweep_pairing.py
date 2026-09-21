"""Layer: L3
WP: WP5
Claims: C3
Status: stub
Requirements: INS-05, INS-09, INS-15

Up-sweep and down-sweep pairing: unbiased range from the half-sum, radial
velocity from the half-difference.  Not implemented.  This is the principal
contribution of the dissertation and it is the module the concept turns on.

The geometry
------------
The platform is in a low Earth orbit at about 7.6 km/s and it sounds downward.
For a vertical sounding the line of sight ``u_los`` points from the spacecraft
toward the reflection point, which is local down, and the closing velocity is
the projection of the platform velocity minus the reflector velocity onto that
line:

    v_close = dot(v_platform - v_reflector, u_los)

For a circular orbit and an exactly vertical line of sight the orbital velocity
is transverse and contributes nothing.  What contributes is everything that
makes the real case not that: orbit eccentricity, the oblateness of the Earth,
attitude pointing error, the vertical plasma drift of the reflecting layer, and
the fact that the reflection point is not directly below the spacecraft when
the ionosphere has a horizontal gradient.  The residual line-of-sight velocity
is small compared with the orbital speed and it is not small compared with what
INS-09 allows.

The reflector's own motion enters through the same projection as the
platform's.  A vertical plasma drift and a radial platform motion are not
separable from one line of sight.  That is a property of the measurement and
not of the estimator, and any velocity this module returns is the sum of the
two.  The simulator can separate them because it holds the truth; the
instrument cannot.

Why the half-sum cancels the bias to first order
-------------------------------------------------
In a linear sweep of rate ``K_f`` the de-chirped beat frequency for a one-way
group range ``R`` is ``f_R = 2 K_f R / c``.  Motion adds the two-way Doppler
shift ``f_d = 2 v_close f / c`` to the received instantaneous frequency, and the
processor, having no way to tell the two contributions apart within one sweep,
converts the sum back to range.  Writing it with the range rate
``Rdot = -v_close``, the apparent range on a sweep of direction ``d`` is

    R_app(d) = R + d * Rdot * f / K_f

where ``dd = Rdot * f / K_f`` is :func:`tss.constants.doppler_range_bias` and
``d`` is +1 for an up-sweep and -1 for a down-sweep.  The derivation and the
sign are in :mod:`tss.l1_propagation.doppler`; that module is the forward model
this one has to invert, and it is deliberately complete so that the inversion
is testable against a known answer.

Pair an up-sweep with a down-sweep taken at the same sounding frequency and the
two apparent ranges are

    R_up   = R + dd
    R_down = R - dd

so the half-sum and the half-difference separate them exactly:

    R  = (R_up + R_down) / 2
    dd = (R_up - R_down) / 2      and therefore     Rdot = dd * K_f / f

The cancellation is exact only because ``dd`` enters the two sweeps with equal
magnitude and opposite sign, and that in turn holds only while ``Rdot``, ``f``
and ``K_f`` are the same for both members of the pair.  Nothing else about the
ionosphere, the absorption or the antenna has to be known.  That is the appeal
of the method and it is why it works with one transmit and one receive pair,
without a second antenna, a second receiver or an external velocity reference.

The residual, which is the actual research
-------------------------------------------
The three quantities that must match between the members of a pair do not
match, and the second-order terms are what INS-09 has to be scored against.

1.  **The velocity changes between the sweeps.**  The two sweeps are separated
    by at least one sweep period ``T_sweep``.  If the line-of-sight
    acceleration is ``a``, the range rates differ by ``a T_sweep`` and the
    half-sum retains

        delta_R ~ (a * T_sweep / 2) * f / K_f

    which is :func:`tss.l1_propagation.doppler.residual_bias_second_order_m`.
    That expression is an order-of-magnitude bound, not an error model: it
    assumes the acceleration is constant across the pair and ignores the
    variation of the line of sight during each sweep.  Making it an error model
    is part of WP5.

2.  **The sounding frequency is not the same at the same instant of an up-sweep
    and a down-sweep.**  A given frequency is visited early in one sweep and
    late in the other, so the pair is separated in time by an amount that
    depends on where in the band the frequency sits: near the middle of the
    band the separation is about one sweep period, and at the band edges it
    approaches two.  The residual is therefore frequency dependent and largest
    at the edges, which is the opposite of the first-order bias, which is
    largest at the top of the band.  The two effects do not have the same
    shape in frequency and cannot be lumped.

3.  **The range itself changes between the sweeps.**  The reflection point
    moves, both because the ionosphere evolves and because the spacecraft has
    flown several kilometres along track.  The half-sum then returns the mean
    of two different ranges rather than one range measured twice, and the
    difference is an along-track average over a distance comparable to the
    along-track sampling of SCI-09.

What INS-09 requires
--------------------
INS-09 is scored as the residual apparent-range bias remaining after this
separation, expressed in range cells of INS-05, with a threshold and a goal of
1.0.  INS-05's threshold cell is 10 km and its goal cell is 5 km, so the
requirement is that the residual stay under one cell of whichever is being
claimed.  At the baseline ``K_f`` of 5.6 MHz/s and the top of the band at
20 MHz, one metre per second of unremoved line-of-sight velocity is 3.6 m of
range bias, so a 10 km cell is exhausted by about 2.8 km/s of unremoved
velocity and a 5 km cell by about 1.4 km/s.  The first-order term is therefore
enormous compared with the requirement and the whole margin lives in how well
the pairing cancels it.  That is why the residual, and not the bias, is what
gets scored, and it is why a bound that can be evaluated from an ephemeris is
worth having before the estimator exists.

The companion conditions, plus or minus 1 km/s unambiguous and 100 m/s
resolution, are INS-09a and INS-09b in the full Rev E set and are not yet in
``config/requirements/requirements.yaml``.  They constrain the pairing too: the
unambiguous velocity follows from how large a half-difference can be
interpreted without the two traces being mis-associated.

The association problem
-----------------------
Everything above assumes the two members of a pair are known to be the same
echo.  They are not known to be.  The up-sweep trace and the down-sweep trace
are displaced in opposite directions by an amount that is unknown before the
estimate is made, so associating them is circular: the association needs the
velocity and the velocity needs the association.  On a clean single trace the
circularity is broken by continuity, since the two traces have the same shape
and differ by a constant offset at each frequency.  Near the cusp, with O and X
superposed and the trace spread, it is not obviously breakable at all.  This is
the part that is genuinely unsolved and it is the reason the module is a stub
rather than ten lines of arithmetic.

Raises
------
NotImplementedError
    On every call.  WP5.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from tss.interfaces import Ionogram, PlatformState, SweepSchedule

__all__ = [
    "SweepPair",
    "pair_sweeps",
    "unbiased_range_m",
    "range_rate_from_pair_mps",
    "associate_traces",
    "residual_bias_cells",
]

_MESSAGE = (
    "WP5: up and down sweep pairing is not implemented.  This is claim C3, the "
    "principal contribution, and it is unstarted research rather than unfinished "
    "code.  The half-sum arithmetic is stated in the module docstring of "
    "tss.l3_signal.novel.sweep_pairing; what is missing is the association of the "
    "two traces, the second-order residual error model, and the demonstration "
    "that the residual stays inside one INS-05 range cell as INS-09 requires.  "
    "See docs/work-packages.md row WP5 and tss.l1_propagation.doppler for the "
    "forward model this has to invert."
)


class SweepPair:
    """One up-sweep and one down-sweep, associated frequency by frequency.

    Intended to carry the two ionograms, the association between their trace
    points, the recovered range and range rate, and the residual bias estimate.
    It is declared here so that the shape of the answer is fixed before the
    estimator exists, which is the same reason
    ``config/campaigns/c3_doppler.yaml`` names a tier that does not work yet.

    Constructing one raises, because a :class:`SweepPair` that could be built
    without the estimator would be an empty container that downstream code
    could read zeros out of.
    """

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        del args, kwargs
        raise NotImplementedError(_MESSAGE)


def associate_traces(
    up: Ionogram,
    down: Ionogram,
    *,
    max_offset_m: float | None = None,
) -> np.ndarray:
    """Match up-sweep trace points to down-sweep trace points.  Not implemented.

    Intended to return an index array giving, for each frequency column of
    ``up``, the column of ``down`` that holds the same echo, or -1 where no
    association could be made.  The circularity described in the module
    docstring is the whole difficulty: the displacement between the two traces
    is the quantity being measured.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del up, down, max_offset_m
    raise NotImplementedError(_MESSAGE)


def unbiased_range_m(
    range_up_m: np.ndarray,
    range_down_m: np.ndarray,
) -> np.ndarray:
    """Range from the half-sum of a paired up and down sweep.  Not implemented.

    The arithmetic is ``(R_up + R_down) / 2``.  It is not implemented as a
    one-line function because a function that takes two range arrays and
    averages them would be used on unassociated traces, would return a number,
    and the number would look like a measurement.  The association is the part
    that does not exist, so nothing downstream of it is offered.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del range_up_m, range_down_m
    raise NotImplementedError(_MESSAGE)


def range_rate_from_pair_mps(
    range_up_m: np.ndarray,
    range_down_m: np.ndarray,
    freq_hz: np.ndarray,
    k_f_hz_per_s: float,
) -> np.ndarray:
    """Range rate from the half-difference.  Not implemented.

    The arithmetic is ``Rdot = (R_up - R_down) K_f / (2 f)``, and the closing
    velocity is its negative.  Remember that what comes back is the sum of the
    platform's radial motion and the reflecting layer's vertical drift; the two
    are not separable from one line of sight, which is stated in
    :mod:`tss.l1_propagation.doppler` as a property of the measurement.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del range_up_m, range_down_m, freq_hz, k_f_hz_per_s
    raise NotImplementedError(_MESSAGE)


def residual_bias_cells(
    pair: SweepPair,
    range_resolution_m: float,
) -> np.ndarray:
    """Residual apparent-range bias in INS-05 range cells.  Not implemented.

    This is the quantity INS-09 is scored against, with a threshold and a goal
    of 1.0 cell.  It needs the second-order error model described in the module
    docstring and it cannot be estimated from the pair alone, because the
    residual is by construction the part of the bias the pairing did not
    remove.  The intended route is to evaluate the bound from the ephemeris and
    to check it against the simulator, where the truth velocity is known.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del pair, range_resolution_m
    raise NotImplementedError(_MESSAGE)


def pair_sweeps(
    up: Ionogram,
    down: Ionogram,
    *,
    platform: PlatformState | None = None,
    sweep_up: SweepSchedule | None = None,
    sweep_down: SweepSchedule | None = None,
    settings: Mapping[str, Any] | None = None,
) -> Ionogram:
    """Separate range from Doppler using an up and down sweep pair.  Not implemented.

    Intended to return a new :class:`~tss.interfaces.Ionogram` whose
    ``virtual_range_m`` is the unbiased half-sum range, whose ``doppler_mps``
    is populated from the half-difference, and whose
    ``residual_range_bias_m`` carries the part of the bias the estimator knows
    it has not removed, so that L4 and L5 can propagate it instead of silently
    absorbing it into the retrieved profile.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del up, down, platform, sweep_up, sweep_down, settings
    raise NotImplementedError(_MESSAGE)
