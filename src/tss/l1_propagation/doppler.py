"""Layer: L1
WP: WP2, WP5
Claims: C3
Status: implemented
Requirements: INS-09, SCI-23

The forward model C3 has to invert: platform-motion Doppler and the
apparent-range bias it puts on an FMCW measurement.

This module is deliberately not a stub.  The generation of the bias is
geometry and arithmetic and both are settled; what is unsettled is the
inversion, which lives in :mod:`tss.l3_signal.novel.sweep_pairing`.  Keeping the
forward model honest and complete here is what makes the inversion testable at
all, because the test of C3 is that the estimator recovers the range and the
velocity this module put in.

Geometry and sign conventions
-----------------------------
``u_los`` is the unit vector pointing **from the platform toward the reflection
point**.  For a vertical topside sounding that is local down.  The closing
velocity is

    v_close = dot(v_platform - v_reflector, u_los)

so a positive closing velocity means the two-way path is shortening.  The range
rate is its negative, ``Rdot = -v_close``.  The two-way Doppler shift is

    f_d = 2 * mu_sat * v_close * f / c

positive for a closing geometry, which is the usual radar convention.
``mu_sat`` is the refractive index at the platform, explained next.

The refractive index at the platform
------------------------------------
A topside sounder is inside the plasma it sounds, so the wave it launches has
the local wavelength ``c / (f mu_sat)`` and not the free-space one.  The phase
of the echo is ``-(2 pi f / c) * 2 P`` with ``P`` the phase path, which for a
vertical sounding is the integral of the refractive index from the reflection
level to the platform (Budden 1961, eq. 10.14).  Moving the platform end of
that integral by ``dz`` changes ``P`` by ``mu_sat * dz``, because the integrand
at that end is ``mu_sat``.  The Doppler shift ``-(2 f / c) dP/dt`` is therefore

    f_d = 2 * f * mu_sat * v_close / c

and every relation in this module that converts a velocity to a frequency or a
range bias carries the same factor.  Only relative motion enters, so a rigid
vertical drift of the whole profile carries the factor as well.  A profile that
changes shape, rather than moving rigidly, changes ``P`` through the integrand
and is not described by one velocity.

``mu_sat`` is 1 well above the local plasma frequency and falls to zero at it:
0.87 at twice the plasma frequency at the platform, 0.75 at 1.5 times, 0.42 at
1.1 times.  :func:`refractive_index_at_platform` returns the field-free value
``sqrt(1 - (f_N_sat / f)^2)``, which is the ordinary-mode index when the field
is neglected and is the same approximation the analytic tier makes for group
delay.  Every function takes ``mu_sat`` as a keyword with a default of 1, which
is the free-space relation, so a caller that has no ionosphere state says so by
omission and a caller that has one must pass it.

The apparent-range bias
-----------------------
In a linear FMCW sweep of rate ``K_f`` the de-chirped beat frequency for a
two-way range ``R`` is ``f_R = 2 K_f R / c``.  Motion changes it.  The echo
arrives with instantaneous frequency

    f_rx(t) = f_tx(t) - K_f * tau + f_d

so on an up-sweep the beat is ``f_R - f_d`` and not ``f_R + f_d``: a closing
geometry makes the up-sweep beat frequency smaller.  That is the standard FMCW
result and it is easy to write down backwards.  The processor, which has no way
to tell the two terms apart in a single sweep, converts the difference back to
range:

    R_apparent(up) = c (f_R - f_d) / (2 K_f) = R - mu_sat * v_close * f / K_f

On a down-sweep ``K_f`` changes sign, the beat becomes ``f_R + f_d``, and the
displacement is the other way.  Writing it with the range rate
``Rdot = -v_close`` instead of the closing velocity removes the minus sign:

    R_apparent = R + direction * mu_sat * Rdot * f / K_f

and ``Rdot * f / K_f`` is :func:`tss.constants.doppler_range_bias`, the
free-space value.  Three
properties of that expression drive the whole of C3.

* It does not depend on range.  Every echo in the sweep is displaced by the
  same amount, so it cannot be found by looking at the shape of a single trace.
* It grows with sounding frequency, linearly once ``mu_sat`` has reached 1, so
  it is worst at the top of the band.  At 20 MHz with a 5.6 MHz/s sweep, 1 m/s of closing velocity is 3.6 m
  of range error, and a few hundred metres per second of unmodelled radial
  motion is kilometres.
* It changes sign with the sweep direction.  ``K_f`` is negative on a
  down-sweep, so the same closing velocity displaces the two traces opposite
  ways: short range on the up-sweep, long range on the down-sweep, by the same
  amount.  That sign change is the entire lever C3 pulls on, and it is why
  :class:`~tss.interfaces.SweepSchedule` carries ``direction`` as a field rather
  than leaving the sign implicit in ``k_f_hz_per_s``.

What is and is not included
---------------------------
Included: platform velocity projected on the line of sight, and a reflector
velocity along the same line so that a vertical plasma drift can be modelled.
Not included: the second-order term from the geometry changing during the
sweep, and the along-track term from flying into a horizontally varying
ionosphere.  Neither is negligible for a coherent interval of a second.  The first
is quantified in :func:`residual_bias_second_order_m` so that it can at least be
bounded; the second needs a ray-traced tier and a horizontally structured L0
state, and there is no honest way to put a number on it here.

A caution that belongs in the forward model rather than only in the inversion:
a vertical plasma drift and a radial platform motion enter through the same
projection and are not separable from a single line of sight.  The simulator
can separate them because it knows the truth.  The instrument cannot, and any
retrieved drift from one line of sight is really the sum of the two.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tss.constants import C_LIGHT, doppler_range_bias
from tss.interfaces import PlatformState, SweepSchedule

__all__ = [
    "DopplerGeometry",
    "vertical_line_of_sight",
    "refractive_index_at_platform",
    "closing_velocity_mps",
    "doppler_shift_hz",
    "velocity_from_doppler_mps",
    "apparent_range_bias_m",
    "residual_bias_second_order_m",
    "coherent_interval_for_velocity_resolution_s",
    "doppler_geometry",
]


@dataclass(frozen=True)
class DopplerGeometry:
    """One frequency's Doppler picture for one echo.

    Attributes
    ----------
    freq_hz
        Sounding frequency.
    closing_velocity_mps
        Line-of-sight closing velocity, positive when the two-way path is
        shortening.
    range_rate_mps
        Its negative, for readers who prefer range rate.
    doppler_hz
        Two-way Doppler shift, positive for a closing geometry.
    range_bias_m
        Apparent-range bias on an up-sweep.  Negate it for a down-sweep.
    sweep_direction
        +1 for the up-sweep the bias was computed for, -1 for a down-sweep.
    mu_sat
        Refractive index at the platform that the shift and the bias were
        computed with.  1.0 is the free-space relation.
    """

    freq_hz: float
    closing_velocity_mps: float
    range_rate_mps: float
    doppler_hz: float
    range_bias_m: float
    sweep_direction: int = 1
    mu_sat: float = 1.0


def vertical_line_of_sight() -> np.ndarray:
    """Unit vector from platform to reflection point for a vertical sounding.

    Returned in the local north-east-down frame, so it is ``(0, 0, 1)``: down.
    This exists so that the vertical case is written once and the choice of
    frame is stated once, rather than a bare ``[0, 0, 1]`` appearing in three
    modules with three different conventions.
    """
    return np.array([0.0, 0.0, 1.0], dtype=float)


def refractive_index_at_platform(
    f_hz: float | np.ndarray,
    f_n_sat_hz: float | np.ndarray,
) -> np.ndarray:
    """Field-free refractive index at the platform, ``sqrt(1 - (f_N_sat/f)^2)``.

    This is the factor the Doppler relations of this module carry; see the
    module docstring for where it comes from.  ``f_n_sat_hz`` is the plasma
    frequency at the platform.

    NaN is returned where ``f_hz`` is at or below ``f_n_sat_hz``.  The ordinary
    wave is evanescent there, there is no echo, and a Doppler shift for it has
    no meaning.  NaN is the same convention the analytic tier uses for the
    other measurable quantities of an evanescent row.

    The magnetic field is neglected.  With the field the index at the platform
    differs between the ordinary and extraordinary modes and depends on the
    angle between the wave normal and the field (Budden 1961, eq. 6.1).  That
    belongs with :func:`tss.l1_propagation.magnetoionic.appleton_hartree`,
    which is a stub.
    """
    f = np.asarray(f_hz, dtype=float)
    f_n = np.asarray(f_n_sat_hz, dtype=float)
    if np.any(f <= 0.0):
        raise ValueError("sounding frequency must be positive")
    if np.any(f_n < 0.0):
        raise ValueError("plasma frequency at the platform must not be negative")
    one_minus_x = 1.0 - (f_n / f) ** 2
    with np.errstate(invalid="ignore"):
        return np.where(one_minus_x > 0.0, np.sqrt(np.abs(one_minus_x)), np.nan)


def _checked_mu(mu_sat: float | np.ndarray) -> np.ndarray:
    """``mu_sat`` as an array, rejecting values outside ``(0, 1]``.

    NaN passes through, because it marks an evanescent row and must propagate
    to the result rather than raise.
    """
    mu = np.asarray(mu_sat, dtype=float)
    finite = mu[np.isfinite(mu)]
    if np.any(finite <= 0.0) or np.any(finite > 1.0):
        raise ValueError(
            "mu_sat must lie in (0, 1]: it is the refractive index at the platform "
            "for a propagating wave below the upper band edge"
        )
    return mu


def closing_velocity_mps(
    platform_velocity: np.ndarray,
    u_los: np.ndarray,
    *,
    reflector_velocity: np.ndarray | None = None,
) -> float:
    """Line-of-sight closing velocity, m/s.

    ``dot(v_platform - v_reflector, u_los)`` with ``u_los`` pointing from the
    platform toward the reflection point.  Positive means the path is
    shortening.

    Both vectors must be expressed in the same frame.  The function normalises
    ``u_los`` rather than trusting the caller, because an un-normalised line of
    sight scales the answer silently and the resulting velocity error looks
    exactly like a physical one.
    """
    v = np.asarray(platform_velocity, dtype=float)
    u = np.asarray(u_los, dtype=float)
    if v.shape != (3,) or u.shape != (3,):
        raise ValueError(
            f"velocity and line of sight must both be 3-vectors; got {v.shape} and {u.shape}"
        )
    norm = float(np.linalg.norm(u))
    if norm == 0.0:
        raise ValueError("line-of-sight vector has zero length")
    u = u / norm
    if reflector_velocity is not None:
        v = v - np.asarray(reflector_velocity, dtype=float)
    return float(np.dot(v, u))


def doppler_shift_hz(
    closing_velocity_mps: float | np.ndarray,
    f_hz: float | np.ndarray,
    *,
    mu_sat: float | np.ndarray = 1.0,
) -> np.ndarray:
    """Two-way Doppler shift, Hz.  ``f_d = 2 mu_sat v f / c``, positive closing.

    ``mu_sat`` is the refractive index at the platform, from
    :func:`refractive_index_at_platform`.  The default of 1 is the free-space
    relation.

    The non-relativistic form is used.  At orbital speeds the relativistic
    correction is of order ``(v/c)^2``, about 6e-10 for 7.5 km/s, which is four
    orders of magnitude below the frequency stability INS-26 asks for and is not
    worth carrying.
    """
    v = np.asarray(closing_velocity_mps, dtype=float)
    f = np.asarray(f_hz, dtype=float)
    return 2.0 * _checked_mu(mu_sat) * v * f / C_LIGHT


def velocity_from_doppler_mps(
    doppler_hz: float | np.ndarray,
    f_hz: float | np.ndarray,
    *,
    mu_sat: float | np.ndarray = 1.0,
) -> np.ndarray:
    """Invert :func:`doppler_shift_hz`.  ``v = c f_d / (2 f mu_sat)``.

    An estimator that leaves ``mu_sat`` at 1 inside the plasma underestimates
    the velocity by the factor ``mu_sat``.
    """
    f_d = np.asarray(doppler_hz, dtype=float)
    f = np.asarray(f_hz, dtype=float)
    return C_LIGHT * f_d / (2.0 * f * _checked_mu(mu_sat))


def apparent_range_bias_m(
    closing_velocity_mps: float | np.ndarray,
    f_hz: float | np.ndarray,
    k_f_hz_per_s: float | np.ndarray,
    *,
    direction: int = 1,
    mu_sat: float | np.ndarray = 1.0,
) -> np.ndarray:
    """Apparent-range bias from Doppler in a linear sweep, metres.

    ``direction * mu_sat * Rdot * f / K_f`` with ``Rdot = -v_close``.  The
    free-space part is evaluated through
    :func:`tss.constants.doppler_range_bias` so that the relation lives in one
    place, and ``mu_sat`` is the refractive index at the platform.  ``direction`` is +1 for an up-sweep and -1 for a down-sweep,
    matching :attr:`tss.interfaces.SweepSchedule.direction`.

    The sign is the thing to read carefully.  A closing geometry makes the
    up-sweep beat frequency smaller, so the up-sweep trace appears at shorter
    range than it should and the down-sweep trace appears at longer range by
    the same amount.  That is what makes the half-sum of the two an unbiased
    range, and it is the property C3 rests on.

    ``k_f_hz_per_s`` is the magnitude of the sweep rate.  Passing a negative
    rate and a direction of -1 together would apply the sign twice, so the
    function takes the absolute value and the sign comes from ``direction``
    alone.  That is a choice, and it is made this way because the sign of the
    bias is the one thing C3 must not get wrong by accident.
    """
    if direction not in (1, -1):
        raise ValueError(f"direction must be +1 or -1; got {direction!r}")
    range_rate = -np.asarray(closing_velocity_mps, dtype=float)
    magnitude = doppler_range_bias(
        range_rate, f_hz, np.abs(np.asarray(k_f_hz_per_s, dtype=float))
    )
    return float(direction) * _checked_mu(mu_sat) * magnitude


def residual_bias_second_order_m(
    acceleration_mps2: float | np.ndarray,
    f_hz: float | np.ndarray,
    k_f_hz_per_s: float | np.ndarray,
    sweep_period_s: float,
    *,
    mu_sat: float | np.ndarray = 1.0,
) -> np.ndarray:
    """Bound on the bias left after up and down sweeps are averaged, metres.

    The half-sum of an up-sweep and a down-sweep cancels the first-order bias
    exactly, because the two biases are equal and opposite.  What survives is
    the change in closing velocity between the two sweeps.  If the line-of-sight
    acceleration is ``a`` and the two sweeps are separated by one sweep period
    ``T``, the velocities differ by ``a T`` and the residual in the half-sum is

        ``delta_R ~ mu_sat * (a * T / 2) * f / K_f``

    which is the first-order expression evaluated at half the velocity change.
    This is an order-of-magnitude bound and not an error model: it assumes the
    acceleration is constant across the pair and it ignores the variation of
    the line of sight during a sweep.

    It is here because INS-09 asks for the residual range bias to stay inside
    one range cell, and a bound that can be evaluated from an ephemeris is what
    turns that requirement into something a campaign can score before the
    estimator of WP5 exists.
    """
    a = np.asarray(acceleration_mps2, dtype=float)
    delta_v = a * float(sweep_period_s) / 2.0
    return _checked_mu(mu_sat) * doppler_range_bias(
        delta_v, f_hz, np.abs(np.asarray(k_f_hz_per_s, dtype=float))
    )


def coherent_interval_for_velocity_resolution_s(
    delta_v_mps: float | np.ndarray,
    f_hz: float | np.ndarray,
    *,
    mu_sat: float | np.ndarray = 1.0,
) -> np.ndarray:
    """Coherent interval needed to resolve a velocity difference, seconds.

    Resolving two Doppler shifts that differ by ``2 mu_sat delta_v f / c``
    needs a coherent interval of at least the reciprocal of that difference:

        ``T >= c / (2 delta_v f mu_sat)``

    With ``mu_sat`` at 1 the product ``f T = c / (2 delta_v)`` is independent
    of sounding frequency, which is the statement that the demand is on the
    velocity resolution alone.  Inside the plasma the interval lengthens by
    ``1 / mu_sat``, so the demand is hardest just above the plasma frequency at
    the platform and not simply at the low end of the band.  For the 50 m/s of SCI-23 the product is about 3.0e6, so
    1 MHz needs 3.0 s and 10 MHz needs 0.30 s.  The same relation appears in
    :mod:`tss.l3_signal.novel.motion_compensation`, where its consequence is
    that the geometry moves materially inside one velocity estimate.
    """
    dv = np.asarray(delta_v_mps, dtype=float)
    f = np.asarray(f_hz, dtype=float)
    return C_LIGHT / (2.0 * dv * f * _checked_mu(mu_sat))


def doppler_geometry(
    platform: PlatformState,
    freq_hz: float | np.ndarray,
    sweep: SweepSchedule,
    *,
    u_los: np.ndarray | None = None,
    reflector_velocity: np.ndarray | None = None,
    f_n_sat_hz: float | None = None,
) -> list[DopplerGeometry]:
    """Assemble the Doppler picture for a set of sounding frequencies.

    ``platform.velocity_ecef_mps`` and ``u_los`` must be in the same frame.  The
    interface names the field ``velocity_ecef_mps``; for a vertical sounding the
    caller usually supplies the local north-east-down components and the default
    line of sight of :func:`vertical_line_of_sight`, in which case only the
    vertical component contributes.  That is stated here rather than enforced,
    because enforcing it would need an attitude and a geodetic conversion that
    belong in the platform package and not in L1.

    ``f_n_sat_hz`` is the plasma frequency at the platform.  When it is given
    the shift and the bias carry the refractive index at the platform from
    :func:`refractive_index_at_platform`, and a frequency at or below it gets
    NaN for both.  When it is ``None`` the free-space relations are used, and
    each returned record says so through ``mu_sat = 1.0``.

    A platform with no velocity gives zero Doppler and zero bias, which is the
    correct answer for a stationary sounder and is not a masked result.
    """
    freqs = np.atleast_1d(np.asarray(freq_hz, dtype=float))
    direction = int(np.sign(sweep.direction)) or 1
    k_f = abs(float(sweep.k_f_hz_per_s))

    if platform.velocity_ecef_mps is None:
        v_close = 0.0
    else:
        los = vertical_line_of_sight() if u_los is None else np.asarray(u_los, dtype=float)
        v_close = closing_velocity_mps(
            platform.velocity_ecef_mps, los, reflector_velocity=reflector_velocity
        )

    # For a stepped sweep the bias uses the rate of the step that sounds each
    # frequency, which below 3 MHz is slower than the peak and so gives a larger
    # bias per unit velocity than k_f_hz_per_s alone would.
    from tss import sweep_law

    if sweep_law.is_stepped(sweep):
        k_by_f = np.atleast_1d(
            sweep_law.rate_hz_per_s(sweep, sweep_law.time_at_frequency_s(sweep, freqs))
        )
    else:
        k_by_f = np.full(freqs.shape, k_f)

    if f_n_sat_hz is None:
        mu_by_f = np.ones(freqs.shape)
    else:
        mu_by_f = np.atleast_1d(refractive_index_at_platform(freqs, float(f_n_sat_hz)))

    out: list[DopplerGeometry] = []
    for f, k_here, mu in zip(freqs, k_by_f, mu_by_f, strict=True):
        f_d = float(doppler_shift_hz(v_close, float(f), mu_sat=float(mu)))
        bias = float(
            apparent_range_bias_m(
                v_close, float(f), float(k_here), direction=direction, mu_sat=float(mu)
            )
        )
        out.append(
            DopplerGeometry(
                freq_hz=float(f),
                closing_velocity_mps=float(v_close),
                range_rate_mps=-float(v_close),
                doppler_hz=f_d,
                range_bias_m=bias,
                sweep_direction=direction,
                mu_sat=float(mu),
            )
        )
    return out
