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

What is **not** here is the full Appleton-Hartree refractive index with
collisions.  That is :func:`appleton_hartree`, and it raises.  The reason it is
a stub rather than a half-implementation is that the quantity the sounder
actually measures is the group path, and the group path for the X mode near its
cutoff is governed by the derivative of the Appleton-Hartree index with respect
to frequency, including the collision term and the angle between the ray and the
field.  A field-free group index dressed up with an X-mode cutoff altitude would
produce an X trace of about the right shape and entirely the wrong retardation
near the cutoff, and a plausible wrong number is worse here than an absent one.
The analytic tier is explicit about this: it computes the X reflection altitude
from the cutoff, computes the group delay with the field-free group index, and
flags every X row it emits so that no result quietly depends on it.
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


def appleton_hartree(
    f_hz: np.ndarray | float,
    f_p_hz: np.ndarray | float,
    f_ce_hz: np.ndarray | float,
    theta_rad: np.ndarray | float,
    nu_per_s: np.ndarray | float,
    mode: str = "O",
) -> np.ndarray:
    """Full Appleton-Hartree complex refractive index with collisions.

    Not implemented.

    Approach
    --------
    The index is

        n^2 = 1 - X / (1 - iZ - Y_T^2 / (2 (1 - X - iZ))
                       -/+ sqrt( Y_T^4 / (4 (1 - X - iZ)^2) + Y_L^2 ))

    with ``X = (f_p/f)^2``, ``Y = f_ce/f``, ``Y_L = Y cos(theta)``,
    ``Y_T = Y sin(theta)`` and ``Z = nu / omega``, the upper sign giving the
    ordinary mode.  Writing that expression down is not the work.  Three things
    are:

    1.  Branch selection.  The square root has two branches and the correct one
        swaps across the quasi-longitudinal to quasi-transverse transition.  A
        naive principal-branch evaluation produces a discontinuity in the middle
        of the band that looks like a physical feature.
    2.  The group index.  What the sounder measures is the group path, which is
        ``d(f n)/df``.  Near the X cutoff that derivative is large and the
        finite-difference evaluation of it is ill-conditioned exactly where the
        answer matters most.  An analytic derivative of the full expression is
        the defensible route and it is algebra that has to be checked against
        something.
    3.  Verification.  The check is PHaRLAP through the MATLAB bridge and the
        author's independent three-dimensional Haselgrove raytracer from the
        Chartier reproduction.  Until that comparison exists there is no way to
        tell a sign error from a physical effect.

    The open question this function will answer is how much of the O-X trace
    separation, and of the differential absorption between the modes, survives
    at the bottom of the sounding band where ``f`` is comparable with ``f_ce``
    and the quasi-longitudinal approximation fails.  That matters for INS-08 and
    for C5, because polarization labelling of the two traces is only useful over
    the band where the two traces are actually separable.

    Raises
    ------
    NotImplementedError
        Always.  WP2, docs/architecture.md section "L1 Propagation".
    """
    del f_hz, f_p_hz, f_ce_hz, theta_rad, nu_per_s, mode
    raise NotImplementedError(
        "WP2: the full Appleton-Hartree index with collisions is not implemented. "
        "See docs/architecture.md section 'L1 Propagation' and the notes in "
        "tss.l1_propagation.magnetoionic.  The analytic tier uses the field-free "
        "group index and flags every X-mode row it emits; it does not silently "
        "substitute for this function."
    )
