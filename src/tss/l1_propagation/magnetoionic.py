"""Layer: L1
WP: WP2
Claims: C1, C5
Status: partial
Requirements: INS-01, INS-08

Magnetoionic splitting: where the O and the X mode reflect, and how far apart
their traces sit.

What is here is algebra and it is exact.  A magnetized plasma supports two
characteristic modes with different cutoffs.  The ordinary mode reflects where
the sounding frequency equals the local plasma frequency.  The extraordinary
mode reflects at the higher cutoff

    f_x = f_ce/2 + sqrt((f_ce/2)^2 + f_p^2)

which is implemented once in :func:`tss.constants.x_mode_cutoff` and imported
here rather than rewritten.  Because ``f_x > f_p`` always, the X mode reaches
its cutoff at a lower plasma density than the O mode does.  In a topside
sounding, where density increases downward toward the F2 peak, a lower density
means a higher altitude, so at a given sounding frequency the X echo comes back
from higher up and at shorter range than the O echo.  Read the other way round,
at a given virtual range the X trace sits at higher frequency, which is the
familiar statement that the X trace lies to the right of the O trace on an
ionogram.  That ordering is what makes polarization labelling of the two traces
workable at all, and it is the basis of INS-08 and of C5.

The third root, the Z mode at ``-f_ce/2 + sqrt((f_ce/2)^2 + f_p^2)``, is
included because it is the same algebra and because Z traces are visible in the
Alouette and ISIS records this work validates against.  Nothing in the current
pipeline generates Z echoes; the function is here so that a study can ask where
they would fall.

The refractive index and the group index
----------------------------------------
:func:`refractive_index_squared` is the collisionless Appleton-Hartree index
for a wave normal at an angle ``theta`` to the field, and :func:`group_index`
is its group index ``d(f mu)/df`` in closed form.  Both follow Budden (1961):
the index is his equations 12.2 and 12.3, which are 6.1 with collisions
removed, and the group index is his 12.8.  The upper sign of the square root is
the ordinary wave.  The closed form matters because what the sounder measures
is the group path, and a finite difference of ``f mu`` is ill-conditioned near
reflection, exactly where the answer matters most.

For a vertical sounding the wave normal is vertical, so ``theta`` is the angle
between the field and the vertical, and the virtual range is the integral of
the group index over height (Budden's 12.21).  That holds although the ray
itself is deviated in the magnetic meridian (his section 13.13): the lateral
deviation moves the reflection point by tens of kilometres and does not change
the virtual range.  :func:`group_range_m` evaluates the integral and
:func:`field_at` supplies the gyrofrequency and the angle from the state.

Near its reflection level the ordinary group index goes as
``cosec(theta) / sqrt(1 - X)`` (Budden 12.12), against ``1 / sqrt(1 - X)`` with
no field.  That factor is why a field-free inversion places a topside profile
tens of kilometres wrong where the field is steep, and it is what
:mod:`tss.l4_retrieval.lamination` now corrects for.

:func:`appleton_hartree` is the same index with collisions, for absorption
studies.  It takes the principal square root, which is the right branch away
from the coupling region near ``X = 1`` at small ``theta`` and is not
guaranteed inside it.

What is still not here: the branch below the gyrofrequency that reflects at
``X = 1 + Y``, which topside ionograms show as the Z trace, and oblique wave
normals.
"""

from __future__ import annotations

from typing import Literal

import numpy as np

from tss.constants import gyrofrequency, x_mode_cutoff
from tss.interfaces import IonosphereState
from tss.l1_propagation.base import (
    EVANESCENT,
    PENETRATING,
    peak_altitude_m,
    plasma_frequency_at,
)

__all__ = [
    "MODES",
    "Mode",
    "o_mode_cutoff",
    "x_mode_cutoff",
    "z_mode_cutoff",
    "cutoff_frequency",
    "gyrofrequency_at",
    "reflects",
    "reflection_altitude_m",
    "trace_separation_hz",
    "appleton_hartree",
    "refractive_index_squared",
    "refractive_index",
    "group_index",
    "field_at",
    "group_range_m",
]

Mode = Literal["O", "X"]

#: The modes the current pipeline generates.  Z is computable but not emitted.
MODES: tuple[str, ...] = ("O", "X")


def o_mode_cutoff(f_p_hz: np.ndarray | float, f_ce_hz: np.ndarray | float = 0.0):
    """Ordinary-mode cutoff, which is the plasma frequency itself.

    The gyrofrequency argument is accepted and ignored so that this has the same
    signature as :func:`x_mode_cutoff` and the two can be selected by name.  The
    O-mode cutoff genuinely does not depend on the field.
    """
    del f_ce_hz
    return np.asarray(f_p_hz, dtype=float)


def z_mode_cutoff(f_p_hz: np.ndarray | float, f_ce_hz: np.ndarray | float):
    """Z-mode cutoff, ``-f_ce/2 + sqrt((f_ce/2)^2 + f_p^2)``.

    Always below the plasma frequency, so a Z trace lies to the left of the O
    trace.  Provided for interpreting Alouette and ISIS ionograms, which show
    them; the simulator does not generate Z echoes.
    """
    f_p = np.asarray(f_p_hz, dtype=float)
    half_ce = np.asarray(f_ce_hz, dtype=float) / 2.0
    return -half_ce + np.sqrt(half_ce**2 + f_p**2)


def cutoff_frequency(
    mode: str,
    f_p_hz: np.ndarray | float,
    f_ce_hz: np.ndarray | float,
):
    """Cutoff frequency for the named mode, Hz.

    ``mode`` is ``"O"``, ``"X"`` or ``"Z"``.  This is the frequency at which a
    wave of that mode reflects from plasma of the given local plasma frequency
    and gyrofrequency.
    """
    key = str(mode).upper()
    if key == "O":
        return o_mode_cutoff(f_p_hz, f_ce_hz)
    if key == "X":
        return x_mode_cutoff(f_p_hz, f_ce_hz)
    if key == "Z":
        return z_mode_cutoff(f_p_hz, f_ce_hz)
    raise ValueError(f"unknown magnetoionic mode {mode!r}, expected 'O', 'X' or 'Z'")


def gyrofrequency_at(
    state: IonosphereState,
    alt_m: np.ndarray | float,
    *,
    default_b_tesla: float | None = None,
) -> np.ndarray:
    """Electron gyrofrequency in Hz at arbitrary altitudes.

    The field magnitude is taken from ``state.b_field_t``, which the interface
    declares as ``(3, N)`` in local north-east-down components on the state's
    own altitude grid.  A ``(N,)`` array is accepted and read as a magnitude
    that has already been taken.

    A state that carries no field and no ``default_b_tesla`` raises.  Returning
    zero instead would silently collapse the X mode onto the O mode, and the two
    traces coinciding is a result a reader would take seriously.
    """
    alt = np.asarray(state.alt_m, dtype=float)
    query = np.asarray(alt_m, dtype=float)

    if state.b_field_t is None:
        if default_b_tesla is None:
            raise ValueError(
                "IonosphereState carries no b_field_t and no default was given, so "
                "the gyrofrequency is unknown.  Supply a field or pass "
                "default_b_tesla explicitly; a field of zero would merge the X "
                "trace into the O trace without saying so."
            )
        return gyrofrequency(np.full(query.shape, float(default_b_tesla)))

    b = np.asarray(state.b_field_t, dtype=float)
    if b.ndim == 2:
        if b.shape[0] != 3 or b.shape[1] != alt.size:
            raise ValueError(
                f"IonosphereState.b_field_t must be (3, {alt.size}); got {b.shape}"
            )
        magnitude = np.sqrt(np.sum(b**2, axis=0))
    elif b.ndim == 1:
        if b.shape != alt.shape:
            raise ValueError("a one-dimensional b_field_t must be parallel to alt_m")
        magnitude = np.abs(b)
    else:
        raise ValueError(f"b_field_t must be (3, N) or (N,); got {b.shape}")

    return gyrofrequency(np.interp(query, alt, magnitude))


def reflects(
    mode: str,
    f_hz: np.ndarray | float,
    f_p_hz: np.ndarray | float,
    f_ce_hz: np.ndarray | float,
) -> np.ndarray:
    """Reflection condition: does a wave of ``f_hz`` reflect in this plasma?

    True where the sounding frequency is at or below the mode's cutoff.  For the
    O mode that is ``f <= f_p``.  For the X mode it is ``f <= f_x``, which
    rearranges to ``f_p^2 >= f * (f - f_ce)``.  Below the gyrofrequency that
    inequality holds at any density, but it does not mean an X wave reflects
    there: ``f_x`` always exceeds ``f_ce``, so there is no X = 1 - Y cutoff
    below it.  The branch that exists below the gyrofrequency is the whistler
    mode, which propagates.  Callers should treat a True result for the X mode
    at ``f < f_ce`` as "no X-mode echo", not as a reflection.
    """
    f = np.asarray(f_hz, dtype=float)
    return np.asarray(f <= cutoff_frequency(mode, f_p_hz, f_ce_hz))


def reflection_altitude_m(
    state: IonosphereState,
    f_hz: float,
    h_sat_m: float,
    mode: str = "O",
    *,
    h_floor_m: float | None = None,
    tol_m: float = 1e-3,
    interp: str = "linear",
    default_b_tesla: float | None = None,
) -> float | str:
    """Altitude at which a vertical ray of ``f_hz`` in ``mode`` reflects.

    Same contract as
    :func:`tss.l1_propagation.absorption.reflection_altitude_m`, extended to the
    X mode.  Returns an altitude in metres, or the sentinel ``"evanescent"``
    when the cutoff at the spacecraft already exceeds the sounding frequency, or
    ``"penetrating"`` when the ray passes the F2 peak.

    For the O mode this reduces to the plasma-frequency crossing and agrees with
    the absorption module's version to within the bisection tolerance.  For the
    X mode the crossing is on the X cutoff, which depends on altitude through
    both the density and the field.
    """
    f_hz = float(f_hz)
    h_sat = float(h_sat_m)
    h_floor = peak_altitude_m(state) if h_floor_m is None else float(h_floor_m)

    def g(h: float) -> float:
        f_p = float(plasma_frequency_at(state, h, interp=interp))
        f_ce = float(gyrofrequency_at(state, h, default_b_tesla=default_b_tesla))
        return float(cutoff_frequency(mode, f_p, f_ce)) - f_hz

    if g(h_sat) >= 0.0:
        return EVANESCENT
    if g(h_floor) < 0.0:
        return PENETRATING

    lo, hi = h_floor, h_sat
    span = max(h_sat - h_floor, 1.0)
    for _ in range(int(np.ceil(np.log2(span / tol_m))) + 4):
        mid = 0.5 * (lo + hi)
        if g(mid) >= 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def trace_separation_hz(
    f_p_hz: np.ndarray | float,
    f_ce_hz: np.ndarray | float,
) -> np.ndarray:
    """Frequency offset between the X and O cutoffs at the same density.

    ``f_x - f_p``.  At high frequency this tends to ``f_ce/2``, which is the
    familiar half-gyrofrequency separation of the two traces; near the bottom of
    the band, where ``f_p`` is comparable with ``f_ce``, the separation is
    larger than that and the asymptotic form should not be used.  Roughly 0.7
    to 0.8 MHz of separation is typical at F-region altitudes at mid latitudes.
    """
    f_p = np.asarray(f_p_hz, dtype=float)
    return x_mode_cutoff(f_p, f_ce_hz) - f_p


def _mode_sign(mode: str) -> float:
    key = str(mode).upper()
    if key == "O":
        return 1.0
    if key == "X":
        return -1.0
    raise ValueError(f"unknown magnetoionic mode {mode!r}, expected 'O' or 'X'")


def _denominator(x: np.ndarray, y: np.ndarray, theta: np.ndarray, sign: float):
    """Budden's ``D`` of equation 12.3, and the square root inside it.

    ``D = (1 - X) - Y_T^2/2 + sign * sqrt(Y_T^4/4 + Y_L^2 (1 - X)^2)``.  For
    the ordinary wave the two last terms nearly cancel as ``X -> 1``, so the
    difference is formed as ``b^2 / (root + a)`` instead of ``root - a``.
    """
    a = 0.5 * (y * np.sin(theta)) ** 2
    b = y * np.cos(theta) * (1.0 - x)
    root = np.sqrt(a * a + b * b)
    if sign > 0.0:
        with np.errstate(invalid="ignore", divide="ignore"):
            tail = np.where(root + a > 0.0, b * b / (root + a), 0.0)
        return (1.0 - x) + tail, root
    return (1.0 - x) - a - root, root


def refractive_index_squared(
    x: np.ndarray | float,
    y: np.ndarray | float,
    theta_rad: np.ndarray | float,
    mode: str = "O",
) -> np.ndarray:
    """Collisionless Appleton-Hartree ``mu^2`` (Budden 1961, eq. 12.2).

    ``x = (f_N / f)^2``, ``y = f_H / f`` and ``theta_rad`` is the angle between
    the wave normal and the field.  ``mu^2 = 1 - X (1 - X) / D``.  The ordinary
    wave has its zero at ``X = 1`` and the extraordinary wave, for ``y < 1``,
    at ``X = 1 - Y``.  A negative value means the wave is evanescent there.
    At ``X = 1`` exactly the ordinary value is returned as 0.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    theta = np.asarray(theta_rad, dtype=float)
    sign = _mode_sign(mode)
    d, _ = _denominator(x, y, theta, sign)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = 1.0 - x * (1.0 - x) / d
    if sign > 0.0:
        out = np.where(x == 1.0, 0.0, out)
    return out


def refractive_index(
    x: np.ndarray | float,
    y: np.ndarray | float,
    theta_rad: np.ndarray | float,
    mode: str = "O",
) -> np.ndarray:
    """``mu``, the positive root of :func:`refractive_index_squared`.

    NaN where the wave is evanescent.
    """
    m2 = refractive_index_squared(x, y, theta_rad, mode)
    with np.errstate(invalid="ignore"):
        return np.where(m2 >= 0.0, np.sqrt(np.abs(m2)), np.nan)


def group_index(
    x: np.ndarray | float,
    y: np.ndarray | float,
    theta_rad: np.ndarray | float,
    mode: str = "O",
) -> np.ndarray:
    """Collisionless group index ``mu' = d(f mu)/df`` (Budden 1961, eq. 12.8).

    Closed form, with the plasma frequency and the gyrofrequency held fixed
    while the wave frequency varies:

        ``mu mu' D = 1 - X - mu^2
                     + sign * Y_L^2 (1 - X^2)(1 - mu^2)
                       / (2 sqrt(Y_L^2 (1 - X)^2 + Y_T^4 / 4)) + mu^2 D``

    which is Budden's expression with ``D`` moved across.  With no field it
    reduces to ``1 / mu``.  NaN where the wave is evanescent.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    theta = np.asarray(theta_rad, dtype=float)
    sign = _mode_sign(mode)
    d, root = _denominator(x, y, theta, sign)
    with np.errstate(invalid="ignore", divide="ignore"):
        m2 = 1.0 - x * (1.0 - x) / d
        y_l2 = (y * np.cos(theta)) ** 2
        extra = np.where(root > 0.0, 0.5 * y_l2 * (1.0 - x * x) * (1.0 - m2) / root, 0.0)
        rhs = 1.0 - x - m2 + sign * extra
        out = (rhs / d + 1.0) / np.sqrt(m2)
    return np.where(m2 > 0.0, out, np.nan)


def field_at(
    state: IonosphereState,
    alt_m: np.ndarray | float,
    *,
    default_b_tesla: float | None = None,
    default_theta_rad: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Gyrofrequency in Hz and the angle between the field and the vertical.

    The angle is what a vertical wave normal makes with the field, in radians
    between 0 and ``pi/2``: 0 over a magnetic pole, ``pi/2`` at the magnetic
    equator.  It is the complement of the magnitude of the dip angle.

    A ``(3, N)`` field in north-east-down components gives both.  A ``(N,)``
    magnitude gives the gyrofrequency and no direction, and then
    ``default_theta_rad`` is required.  So is it, with ``default_b_tesla``, when
    the state carries no field at all.  A direction is never invented.
    """
    f_h = gyrofrequency_at(state, alt_m, default_b_tesla=default_b_tesla)
    query = np.asarray(alt_m, dtype=float)
    b = None if state.b_field_t is None else np.asarray(state.b_field_t, dtype=float)
    if b is not None and b.ndim == 2:
        alt = np.asarray(state.alt_m, dtype=float)
        magnitude = np.sqrt(np.sum(b**2, axis=0))
        with np.errstate(invalid="ignore", divide="ignore"):
            cos_theta = np.where(magnitude > 0.0, np.abs(b[2]) / magnitude, 1.0)
        theta = np.arccos(np.clip(np.interp(query, alt, cos_theta), 0.0, 1.0))
        return f_h, theta
    if default_theta_rad is None:
        raise ValueError(
            "the state carries no field direction, so the angle between the field "
            "and the vertical is unknown.  Supply b_field_t as (3, N) north-east-down "
            "components or pass default_theta_rad explicitly."
        )
    return f_h, np.full(query.shape, float(default_theta_rad))


def group_range_m(
    state: IonosphereState,
    f_hz: float,
    h_refl_m: float,
    h_sat_m: float,
    mode: str = "O",
    *,
    n_steps: int = 2000,
    interp: str = "linear",
    default_b_tesla: float | None = None,
    default_theta_rad: float | None = None,
) -> float:
    """Virtual range from the reflection altitude to the spacecraft, metres.

    The integral of :func:`group_index` over the vertical path, with the
    gyrofrequency and the field direction taken from the state at each height.
    The substitution ``h = h_refl + (h_sat - h_refl) s^2`` and the midpoint
    rule in ``s`` are the same as in the field-free
    :func:`tss.l1_propagation.analytic.group_range_m`, and for the same reason:
    the group index goes as the inverse square root of the height above the
    reflection level for both modes, and the substitution removes that.

    Samples that interpolation puts marginally inside the reflecting region
    are given the group index of the nearest propagating sample rather than a
    NaN, as the field-free version clips ``1 - X``.
    """
    h_refl = float(h_refl_m)
    h_sat = float(h_sat_m)
    if h_sat <= h_refl:
        return 0.0
    n = int(n_steps)
    s = (np.arange(n, dtype=float) + 0.5) / n
    h = h_refl + (h_sat - h_refl) * s**2
    f_p = plasma_frequency_at(state, h, interp=interp)
    f_h, theta = field_at(
        state, h, default_b_tesla=default_b_tesla, default_theta_rad=default_theta_rad
    )
    mu_g = group_index((f_p / float(f_hz)) ** 2, f_h / float(f_hz), theta, mode)
    bad = ~np.isfinite(mu_g)
    if np.all(bad):
        return float("nan")
    if np.any(bad):
        good = np.flatnonzero(~bad)
        mu_g = mu_g[good[np.searchsorted(good, np.arange(n)).clip(0, good.size - 1)]]
    return float(np.sum(mu_g * 2.0 * (h_sat - h_refl) * s) / n)


def appleton_hartree(
    f_hz: np.ndarray | float,
    f_p_hz: np.ndarray | float,
    f_ce_hz: np.ndarray | float,
    theta_rad: np.ndarray | float,
    nu_per_s: np.ndarray | float,
    mode: str = "O",
) -> np.ndarray:
    """Complex Appleton-Hartree refractive index with collisions.

    Budden (1961) equation 6.1:

        ``n^2 = 1 - X / (U - Y_T^2 / (2 (U - X))
                         +/- sqrt(Y_T^4 / (4 (U - X)^2) + Y_L^2))``

    with ``U = 1 - iZ``, ``Z = nu / omega``, and the upper sign the ordinary
    wave.  The root returned has a non-negative real part, and with the
    ``exp(i omega t)`` convention its imaginary part is non-positive, so
    ``-Im(n)`` is the attenuation.

    The principal square root is taken.  Without collisions that is the right
    branch everywhere.  With collisions it is right away from ``X = 1``, and
    inside the coupling region near ``X = 1`` at small ``theta`` the two waves
    exchange identity and no single-valued choice is correct (Budden chapter
    19).  The collisionless :func:`refractive_index_squared` and
    :func:`group_index` are what the group-delay calculation uses, and they do
    not have that problem.
    """
    f = np.asarray(f_hz, dtype=float)
    x = (np.asarray(f_p_hz, dtype=float) / f) ** 2
    y = np.asarray(f_ce_hz, dtype=float) / f
    theta = np.asarray(theta_rad, dtype=float)
    u = 1.0 - 1j * np.asarray(nu_per_s, dtype=float) / (2.0 * np.pi * f)
    sign = _mode_sign(mode)
    y_t2 = (y * np.sin(theta)) ** 2
    y_l2 = (y * np.cos(theta)) ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        half = y_t2 / (2.0 * (u - x))
        n2 = 1.0 - x / (u - half + sign * np.sqrt(half * half + y_l2))
    n = np.sqrt(n2)
    return np.where(n.real < 0.0, -n, n)
