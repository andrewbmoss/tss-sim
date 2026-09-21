"""Layer: L1
WP: WP2
Claims: C1
Status: implemented
Requirements: INS-10

Deviative and non-deviative ionospheric absorption, integrated rather than
scaled.

The absorption coefficient is

    kappa = (omega / c) * |Im(n)|          Np/m

with the complex refractive index of a collisional, unmagnetized plasma

    n^2 = 1 - X / (1 - i Z),    X = (f_N / f)^2,    Z = nu / omega

integrated along the physical path from the spacecraft down to the reflection
height and doubled for the two-way trip, then converted to decibels with the
factor 8.685889638 nepers per decibel.

This is a reimplementation of the author's tested link budget module
``LinkBudget/tss_link_budget.py`` (functions ``collision_frequency``,
``refractive_index``, ``reflection_height`` and ``absorption_db``), rewritten
against :class:`~tss.interfaces.IonosphereState` instead of that module's
``IonosphereProfile``.  The numerical behaviour is preserved deliberately: same
Coulomb collision form, same complex refractive index, same cubic-stretched
integration grid, same bisection on the reflection height, same nepers-to-
decibels constant.  Results from the two should agree to round-off for the same
profile, and ``tests/unit/test_l1_propagation.py`` checks the properties that
matter rather than a stored table, because the stored table belongs with the
link budget's own 21 tests.

Three points about the formulation are worth keeping in front of a reader.

**It stays finite at the reflection point.**  Where X goes to one, the complex
refractive index gives ``n^2 -> -iZ / (1 - iZ)``, so ``|Im(n)| -> sqrt(Z/2)``,
which is small and finite rather than singular.  That is the whole reason for
using the complex index instead of the real one.  A real-index treatment has to
stop short of reflection and then argue about where to stop.

**Topside absorption is deviative-dominated and rises toward foF2.**  This is
the opposite of the ground-sounding intuition, and getting it backwards is easy.
For a ground sounder the absorbing D and E regions are traversed at high
frequency where X is small, the absorption is non-deviative, and it falls as
roughly f^-2: low frequencies are the ones that suffer.  A topside sounder
never sees the D region at all.  Its path lies entirely in the F region, where
the wave spends its loss budget in the region just above reflection, where X
approaches one, the group retardation is large, and the deviative term
dominates.  The consequence is that absorption grows as the sounding frequency
approaches foF2, because the reflection point moves down toward the peak where
the density, and therefore the Coulomb collision frequency, is highest and the
retarded path is longest.  An f^-2 scaling anchored at one measured value, which
is what the earlier ECE 6930 analysis used, gets both the magnitude and the
trend wrong and diverges at the bottom of the band.

**The magnetic field is neglected here.**  The absorption of the O and X modes
differs because their refractive indices differ, and that difference is not in
this module.  For a link budget the error is small compared with the two decibel
uncertainty on the external noise reference.  When it matters, it belongs with
the full Appleton-Hartree treatment stubbed in
:mod:`tss.l1_propagation.magnetoionic`, not here.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from tss.constants import C_LIGHT, plasma_frequency
from tss.interfaces import IonosphereState
from tss.l1_propagation.base import (
    EVANESCENT,
    PENETRATING,
    collision_frequency_at,
    density_at,
    peak_altitude_m,
    plasma_frequency_at,
)

__all__ = [
    "NEPER_TO_DB",
    "collision_frequency",
    "refractive_index",
    "absorption_coefficient",
    "reflection_altitude_m",
    "absorption_db",
    "absorption_profile",
]

#: Nepers to decibels.  ``20 / ln(10)``.  Spelled out to the digits the link
#: budget module uses so that the two implementations agree to round-off.
NEPER_TO_DB = 8.685889638


def collision_frequency(ne: np.ndarray | float, te: np.ndarray | float) -> np.ndarray:
    """Electron collision frequency in s^-1 from density and temperature.

    Electron-ion Coulomb collisions dominate in the topside F region, so the
    Schunk and Nagy form

        nu_ei = 54.5 * n_i[cm^-3] / Te^1.5

    is used with quasi-neutrality ``n_i = n_e``.  Electron-neutral collisions
    are negligible above roughly 300 km and are not included.  A profile
    extended downward into the E region needs them added, and this function will
    understate the collision frequency there.

    Parameters
    ----------
    ne
        Electron density, m^-3.
    te
        Electron temperature, K.
    """
    ne_cm3 = np.asarray(ne, dtype=float) * 1e-6
    te_arr = np.asarray(te, dtype=float)
    return 54.5 * ne_cm3 / np.power(te_arr, 1.5)


def refractive_index(
    f_hz: np.ndarray | float,
    ne: np.ndarray | float,
    nu: np.ndarray | float,
) -> np.ndarray:
    """Complex refractive index of a collisional, unmagnetized plasma.

        n^2 = 1 - X / (1 - i Z),    X = (f_N / f)^2,    Z = nu / omega

    The branch taken by :func:`numpy.sqrt` on complex input is the principal
    one, with non-negative real part, so the imaginary part is negative for an
    absorbing medium.  Callers take the magnitude; see
    :func:`absorption_coefficient`.
    """
    f = np.asarray(f_hz, dtype=float)
    omega = 2.0 * math.pi * f
    x = (plasma_frequency(ne) / f) ** 2
    z = np.asarray(nu, dtype=float) / omega
    n2 = 1.0 - x / (1.0 - 1j * z)
    return np.sqrt(np.asarray(n2, dtype=complex))


def absorption_coefficient(
    f_hz: float,
    ne: np.ndarray | float,
    nu: np.ndarray | float,
) -> np.ndarray:
    """Absorption coefficient in nepers per metre, ``(omega/c) * |Im(n)|``.

    Exposed separately from the integral so that the behaviour at the
    reflection point can be tested directly.  At X exactly one the coefficient
    is ``(omega/c) * sqrt(Z/2)``, which is finite for any non-zero collision
    frequency and goes to zero as the collision frequency does.
    """
    omega = 2.0 * math.pi * float(f_hz)
    n_cplx = refractive_index(f_hz, ne, nu)
    return (omega / C_LIGHT) * np.abs(np.imag(n_cplx))


def reflection_altitude_m(
    state: IonosphereState,
    f_hz: float,
    h_sat_m: float,
    *,
    h_floor_m: float | None = None,
    tol_m: float = 1e-3,
    interp: str = "linear",
) -> float | str:
    """Altitude at which a vertical ray of frequency ``f_hz`` reflects.

    The reflection condition used here is the ordinary-mode one, that the local
    plasma frequency equals the sounding frequency.  The extraordinary mode
    reflects higher, at the cutoff of
    :func:`tss.constants.x_mode_cutoff`; see
    :func:`tss.l1_propagation.magnetoionic.reflection_altitude_m` for that case.

    The search runs downward from the spacecraft to ``h_floor_m``, which
    defaults to the altitude of the density maximum on the state's own grid.

    Returns
    -------
    float
        Reflection altitude in metres, located by bisection to ``tol_m``.
    ``"evanescent"``
        The plasma frequency at the spacecraft already exceeds ``f_hz``, so the
        wave cannot propagate away from the antenna.  This is the physical floor
        on the sounding band and it is what INS-01's lower edge is set against.
    ``"penetrating"``
        The wave passes the F2 peak.  There is no topside echo.

    The bisection is not decoration.  The absorption integrand rises steeply
    just below the reflection point, so a reflection height quantised to a
    coarse search grid puts visible jitter on the absorption-versus-frequency
    curve, and that jitter then looks like structure in a trade study.
    """
    f_hz = float(f_hz)
    h_floor = peak_altitude_m(state) if h_floor_m is None else float(h_floor_m)
    h_sat = float(h_sat_m)

    def g(h: float) -> float:
        return float(plasma_frequency_at(state, h, interp=interp)) - f_hz

    if g(h_sat) >= 0.0:
        return EVANESCENT
    if g(h_floor) < 0.0:
        return PENETRATING

    # g increases as altitude decreases toward the peak.
    lo, hi = h_floor, h_sat          # g(lo) >= 0, g(hi) < 0
    span = max(h_sat - h_floor, 1.0)
    n_iter = int(np.ceil(np.log2(span / tol_m))) + 4
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        if g(mid) >= 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def absorption_db(
    state: IonosphereState,
    f_hz: float,
    h_sat_m: float,
    *,
    h_refl_m: float | str | None = None,
    n_steps: int = 3000,
    interp: str = "linear",
) -> float:
    """Two-way ionospheric absorption in dB for a vertical sounding.

    Parameters
    ----------
    state
        Truth ionosphere.  Needs ``ne`` and either ``nu`` or ``te``.
    f_hz
        Sounding frequency.
    h_sat_m
        Spacecraft altitude.
    h_refl_m
        Reflection altitude, if already known.  Pass the sentinel strings
        ``"evanescent"`` or ``"penetrating"`` straight through and they are
        handled.  When ``None`` the reflection altitude is found here.
    n_steps
        Number of samples on the stretched integration grid.  The default of
        3000 matches the link budget module.  Convergence is checked in
        ``tests/unit/test_l1_propagation.py``.

    Returns
    -------
    float
        Two-way absorption in decibels, a positive loss.  ``nan`` for an
        evanescent case, where there is no propagation and absorption is not
        defined.  For a penetrating ray the integral is taken down to the peak,
        which gives the one-way absorption to the peak doubled; that is not a
        real echo and the caller is expected to have flagged the row
        ``"penetrating"`` already.

    Notes
    -----
    The integration grid is non-uniform, clustered toward the reflection point:

        ``s = linspace(0, 1, n_steps) ** 3``,  ``h = h_refl + (h_sat - h_refl) * s``

    Near reflection ``Im(n)`` varies as the square root of the distance from the
    reflection height, so a uniform grid converges slowly and, worse,
    non-monotonically in ``n_steps``, which makes a convergence check read as
    noise.  The cubic stretch fixes that at the same cost.  The stretch is in
    the sampling only; the trapezoidal rule is applied in ``h``, so no Jacobian
    appears.
    """
    f_hz = float(f_hz)
    h_sat = float(h_sat_m)

    if h_refl_m is None:
        h_refl_m = reflection_altitude_m(state, f_hz, h_sat, interp=interp)

    if h_refl_m == EVANESCENT:
        return float("nan")
    if h_refl_m == PENETRATING:
        h_refl = peak_altitude_m(state)
    else:
        h_refl = float(h_refl_m)  # type: ignore[arg-type]

    if h_refl >= h_sat:
        return 0.0

    s = np.linspace(0.0, 1.0, int(n_steps)) ** 3
    h = h_refl + (h_sat - h_refl) * s
    ne = density_at(state, h, interp=interp)
    nu = collision_frequency_at(state, h, interp=interp)
    kappa = absorption_coefficient(f_hz, ne, nu)
    one_way_np = float(np.trapezoid(kappa, h))
    return 2.0 * abs(one_way_np) * NEPER_TO_DB


def absorption_profile(
    state: IonosphereState,
    freq_hz: np.ndarray,
    h_sat_m: float,
    *,
    n_steps: int = 3000,
    interp: str = "linear",
) -> dict[str, Any]:
    """Absorption and reflection altitude across a set of frequencies.

    A convenience for studies and for the analytic tier, which needs both
    quantities for every row it emits.  Returns a plain dict with parallel
    arrays ``absorption_db``, ``reflection_alt_m`` (``nan`` where there is no
    reflection) and ``flags`` (an object array of ``""``, ``"evanescent"`` or
    ``"penetrating"``).
    """
    freqs = np.atleast_1d(np.asarray(freq_hz, dtype=float))
    a_db = np.empty(freqs.shape, dtype=float)
    h_refl = np.empty(freqs.shape, dtype=float)
    flags = np.empty(freqs.shape, dtype=object)

    for i, f in enumerate(freqs):
        state_or_alt = reflection_altitude_m(state, float(f), h_sat_m, interp=interp)
        if isinstance(state_or_alt, str):
            flags[i] = state_or_alt
            h_refl[i] = np.nan
        else:
            flags[i] = ""
            h_refl[i] = state_or_alt
        a_db[i] = absorption_db(
            state,
            float(f),
            h_sat_m,
            h_refl_m=state_or_alt,
            n_steps=n_steps,
            interp=interp,
        )

    return {
        "freq_hz": freqs,
        "absorption_db": a_db,
        "reflection_alt_m": h_refl,
        "flags": flags,
    }
