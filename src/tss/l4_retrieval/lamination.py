"""Layer: L4
WP: WP6
Claims: C2
Status: implemented
Requirements: INS-05, SCI-05, SCI-07

Virtual height to true height inversion by lamination.

The classical case is implemented and it is implemented properly: a smooth,
monotonic topside profile inverted downward from the spacecraft, with the
plasma frequency taken as linear in height inside each lamina, which gives a
closed-form group-path integral and a lower-triangular system that solves by
forward substitution.  This is the method Jackson used for the Alouette and
ISIS topside ionograms and it is the standard against which anything else has
to be measured.

Two forms are here.  The field-free form is the classical closed-form
lamination.  The magnetoionic form uses the same laminae and the same forward
substitution with the ordinary-wave group index of the magnetised plasma in
place of the field-free one, which is the construction of Budden (1961),
section 12.13.  The magnetoionic form is the one a profile is scored from,
because the field-free form places a topside profile tens of kilometres wrong
where the field is steep.  See "The magnetoionic form" below.

The geometry
------------
A topside sounder looks down.  The spacecraft sits at altitude ``h_s`` in
plasma of local plasma frequency ``f_0``.  A sounding at frequency ``f`` passes
downward through increasing density until it reaches the level where ``f_N = f``
and reflects.  What the instrument measures is the group range

    R'(f) = integral from h_r(f) to h_s of mu'(f, h) dh

where ``mu' = 1 / sqrt(1 - f_N(h)^2 / f^2)`` is the group refractive index in
the absence of a magnetic field.  The integrand is ``1 / mu`` at the
spacecraft, where ``mu`` is the refractive index there, and diverges at the reflection point, which is why the integral needs a closed form
inside each lamina rather than a quadrature rule.

The lamination
--------------
Order the detected sounding frequencies ``f_1 < f_2 < ... < f_n``, all above
``f_0``.  Each one reflects at a level that is deeper than the last, so the
frequencies partition the profile into ``n`` laminae, lamina ``j`` running from
the level where ``f_N = f_{j-1}`` down to the level where ``f_N = f_j``, with
``f_0`` at the spacecraft.  Let ``t_j`` be the thickness of lamina ``j``.

Inside lamina ``j``, take ``f_N`` linear in height.  Then for a sounding at
``f_k >= f_j`` the contribution of that lamina to the group range is

    integral over lamina j = t_j * f_k * [asin(f_j / f_k) - asin(f_{j-1} / f_k)]
                                       / (f_j - f_{j-1})

which is the substitution ``dh = -t_j df_N / (f_j - f_{j-1})`` applied to
``dh / sqrt(1 - f_N^2/f_k^2)``.  Call that coefficient ``G[k, j]``.  It is
finite at ``j = k``, where ``asin(1) = pi/2``, because the divergence of the
integrand at the reflection point is integrable.  Then

    R'(f_k) = sum over j <= k of G[k, j] * t_j

which is a lower-triangular system.  It solves by forward substitution, one
lamina at a time, each thickness following from the measurement at its own
frequency and the thicknesses already found:

    t_k = ( R'(f_k) - sum over j < k of G[k, j] t_j ) / G[k, k]

The reflection height for ``f_k`` is then ``h_k = h_s - (t_1 + ... + t_k)`` and
the density there is ``Ne = (f_k / 8.98)^2`` in the usual units.

Three properties of this system are worth holding on to.

* It is exact for a profile that really is piecewise linear in ``f_N`` against
  height.  The only error for such a profile is arithmetic.  That makes the
  unit test unambiguous: build such a profile, generate the group ranges from
  the same matrix, and the inversion must return the thicknesses it started
  from to machine precision.
* Errors accumulate downward.  ``t_k`` is measured from ``h_{k-1}``, so an
  error in an early lamina displaces everything below it and never washes out.
  This is the dominant structural feature of the error budget and it is why
  :mod:`tss.l4_retrieval.error_budget` cannot be a per-point uncertainty.
* The starting value is the first lamina's upper boundary.  Getting ``f_0``
  wrong displaces the whole profile.  :mod:`tss.l4_retrieval.base` states the
  three ways of obtaining it and why they are not interchangeable.

The magnetoionic form
---------------------
With a magnetic field the ordinary wave still reflects where ``f_N = f``, so
the laminae and their boundaries are unchanged.  What changes is the group
index, which becomes ``mu'(f, f_N, f_H, theta)`` of
:func:`tss.l1_propagation.magnetoionic.group_index`, with ``theta`` the angle
between the field and the vertical.  The coefficient is then

    G[k, j] = 1 / (f_j - f_{j-1}) * integral from f_{j-1} to f_j of
              mu'(f_k, f_N, f_H_j, theta_j) df_N

which is Budden's equation 12.46 with his 12.54.  It has no closed form.  The
substitution ``f_N = f_k sin(phi)`` removes the singularity at the reflection
level, where ``mu'`` goes as ``cosec(theta) / cos(phi)``, and the integral is
taken by Gauss-Legendre quadrature in ``phi``.  With no field the integrand is
``f_k`` and the coefficient reduces to the arcsine form above.

The gyrofrequency and the angle belong to the lamina, so they are taken at its
mid height from a field model evaluated at the spacecraft's position.  The
height of a lamina is not known until its thickness is, so the diagonal
coefficient is computed with the field at the top of the lamina and corrected
once.  The field varies by about one per cent in 20 km, so one correction is
enough.

The field model is the International Geomagnetic Reference Field through
:mod:`tss.l0_ionosphere.geomag`, which is what an instrument would use.  In the
simulator it is the same model the truth was built with, so the retrieval is
handed a field with no error.  The error of a real field model is a few tenths
of a per cent in magnitude and is not in this budget.

The size of what this corrects: for an exponential topside with 1 MHz at the
spacecraft and a 150 km scale height, the field-free form is in error by up to
34 km when the field is 10 degrees from the vertical, 26 km at 30 degrees and
14 km at 50 degrees.  SCI-07 asks for 20 km.  The error is largest at high
magnetic latitude, where ``theta`` is small.

Limits
------
**The field-free form ignores the magnetic field.**  It records
``flags["magnetoionic"] = False`` so that no result from it is read as a
field-corrected retrieval, and L5 does not score SCI-05 or SCI-07 from it.

**The magnetoionic form inverts the ordinary trace only.**  The extraordinary
trace reflects at a different density for the same frequency and needs the
variable of Budden's section 12.17.  A column labelled as the extraordinary
mode is refused.

**Monotonicity is assumed.**  Every lamina thickness is assumed positive, which
is the statement that density increases downward all the way from the
spacecraft to the F2 peak.  A depletion crossing breaks that, and the algebra
does not notice: it returns a thickness, the thickness is wrong, and the
profile below it is displaced.  A non-positive thickness is the one symptom the
algebra does produce, and it is counted and flagged rather than clamped,
because clamping would remove the only evidence that the assumption failed.
:mod:`tss.l4_retrieval.depletion` is where the non-monotonic case belongs.

**Vertical incidence is assumed.**  Every echo is taken to have come from
directly below.  An off-vertical return from a depletion wall did not, and
inverting it as though it did produces a profile that is wrong and looks
plausible.  That is stated in :mod:`tss.l1_propagation.offvertical` from the
other side.

**Nothing below the lowest reflection is retrieved.**  The profile stops at the
deepest sounding frequency that produced an echo.  Above foF2 the wave
penetrates and there is no echo at all, so the F2 peak and everything below it
are outside a topside sounding by construction.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from tss.interfaces import DensityProfile, Ionogram, PlatformState
from tss.l4_retrieval.base import (
    RetrievalBase,
    density_from_trace,
    starting_plasma_frequency_hz,
    trace_points,
)
from tss.registry import register

__all__ = [
    "LaminationRetrieval",
    "MagnetoionicLaminationRetrieval",
    "FieldModel",
    "constant_field",
    "geomagnetic_field_model",
    "lamination_matrix",
    "lamina_thicknesses_m",
    "virtual_range_from_laminae_m",
    "magnetoionic_lamination_row",
    "lamina_thicknesses_magnetoionic_m",
    "virtual_range_from_laminae_magnetoionic_m",
]

#: Gyrofrequency in Hz and the angle between the field and the vertical in
#: radians, as functions of altitude in metres.
FieldModel = Callable[[np.ndarray], tuple[np.ndarray, np.ndarray]]

_GL_CACHE: dict[int, tuple[np.ndarray, np.ndarray]] = {}


def _gauss_legendre(n: int) -> tuple[np.ndarray, np.ndarray]:
    if n not in _GL_CACHE:
        _GL_CACHE[n] = np.polynomial.legendre.leggauss(n)
    return _GL_CACHE[n]


def constant_field(f_h_hz: float, theta_rad: float) -> FieldModel:
    """A field model with the same gyrofrequency and angle at every height."""

    def model(alt_m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        alt = np.asarray(alt_m, dtype=float)
        return np.full(alt.shape, float(f_h_hz)), np.full(alt.shape, float(theta_rad))

    return model


def geomagnetic_field_model(
    lat_deg: float,
    lon_deg: float,
    epoch: str,
    alt_max_m: float,
    *,
    prefer: str = "auto",
    alt_step_m: float = 1.0e4,
) -> tuple[FieldModel, dict[str, Any]]:
    """Field model along the vertical below a spacecraft, with its provenance.

    Evaluates :func:`tss.l0_ionosphere.geomag.field_ned_t` on a grid from 80 km
    to ``alt_max_m`` and interpolates the gyrofrequency and the cosine of the
    angle between the field and the vertical.  ``prefer`` is passed through:
    ``"auto"`` tries IGRF and falls back to the centred dipole, and the
    provenance says which answered.
    """
    from tss.constants import gyrofrequency
    from tss.l0_ionosphere import geomag

    top = max(float(alt_max_m), 9.0e4)
    grid = np.arange(8.0e4, top + float(alt_step_m), float(alt_step_m))
    b, provenance = geomag.field_ned_t(lat_deg, lon_deg, grid, epoch, prefer=prefer)
    magnitude = np.sqrt(np.sum(b**2, axis=0))
    f_h_grid = gyrofrequency(magnitude)
    cos_grid = np.abs(b[2]) / magnitude

    def model(alt_m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        alt = np.asarray(alt_m, dtype=float)
        f_h = np.interp(alt, grid, f_h_grid)
        theta = np.arccos(np.clip(np.interp(alt, grid, cos_grid), 0.0, 1.0))
        return f_h, theta

    return model, dict(provenance)


def magnetoionic_lamination_row(
    f_k_hz: float,
    lower_hz: np.ndarray,
    upper_hz: np.ndarray,
    f_h_hz: np.ndarray,
    theta_rad: np.ndarray,
    *,
    n_nodes: int = 24,
    n_nodes_reflection: int = 64,
) -> np.ndarray:
    """Group-path coefficients of the laminae a sounding at ``f_k_hz`` crosses.

    ``lower_hz`` and ``upper_hz`` are the plasma frequencies at the top and the
    bottom of each lamina, with the last ``upper_hz`` equal to ``f_k_hz``, the
    reflection level.  ``f_h_hz`` and ``theta_rad`` are the gyrofrequency and
    the angle between the field and the vertical in each lamina.  Returns one
    coefficient per lamina, dimensionless: the group path a lamina of unit
    thickness contributes.

    The last lamina holds the reflection level, where the integrand has its
    structure, and gets ``n_nodes_reflection`` quadrature nodes.  The others
    get ``n_nodes``.
    """
    from tss.l1_propagation.magnetoionic import group_index

    f_k = float(f_k_hz)
    lower = np.asarray(lower_hz, dtype=float)
    upper = np.asarray(upper_hz, dtype=float)
    f_h = np.asarray(f_h_hz, dtype=float)
    theta = np.asarray(theta_rad, dtype=float)
    out = np.empty(lower.size, dtype=float)

    def integrate(sel: slice, n: int) -> None:
        lo, hi = lower[sel], upper[sel]
        if lo.size == 0:
            return
        phi_lo = np.arcsin(np.clip(lo / f_k, -1.0, 1.0))
        phi_hi = np.arcsin(np.clip(hi / f_k, -1.0, 1.0))
        nodes, weights = _gauss_legendre(n)
        half = 0.5 * (phi_hi - phi_lo)[:, None]
        phi = 0.5 * (phi_hi + phi_lo)[:, None] + half * nodes[None, :]
        x = np.sin(phi) ** 2
        mu_g = group_index(x, (f_h[sel] / f_k)[:, None], theta[sel][:, None], "O")
        integrand = mu_g * f_k * np.cos(phi)
        out[sel] = np.sum(integrand * weights[None, :], axis=1) * half[:, 0] / (hi - lo)

    integrate(slice(0, lower.size - 1), int(n_nodes))
    integrate(slice(lower.size - 1, lower.size), int(n_nodes_reflection))
    return out


def _checked_grid(freq_hz: np.ndarray, f_start_hz: float) -> tuple[np.ndarray, np.ndarray]:
    f = np.asarray(freq_hz, dtype=float)
    f0 = float(f_start_hz)
    if f.ndim != 1 or f.size == 0:
        raise ValueError(f"freq_hz must be a non-empty 1-D array; got shape {f.shape}")
    if not np.all(np.diff(f) > 0.0):
        raise ValueError("freq_hz must be strictly increasing")
    if f[0] <= f0:
        raise ValueError(
            f"lowest sounding frequency {f[0]:.6g} Hz is at or below the plasma "
            f"frequency at the spacecraft {f0:.6g} Hz; that wave is evanescent and "
            "cannot have produced an echo"
        )
    return f, np.concatenate(([f0], f))


def lamina_thicknesses_magnetoionic_m(
    virtual_range_m: np.ndarray,
    freq_hz: np.ndarray,
    f_start_hz: float,
    h_sat_m: float,
    field: FieldModel,
) -> np.ndarray:
    """Lamina thicknesses from an ordinary-mode trace, with the field, metres.

    Forward substitution as in :func:`lamina_thicknesses_m`, with the
    coefficients of :func:`magnetoionic_lamination_row`.  The field of each
    lamina is taken at its mid height.  For the lamina being solved that height
    is not yet known, so its coefficient is formed with the field at its top
    and corrected once with the field at the mid height that results.

    No clamping, for the reason given in :func:`lamina_thicknesses_m`.
    """
    r = np.asarray(virtual_range_m, dtype=float)
    f, edges = _checked_grid(freq_hz, f_start_hz)
    n = f.size
    if r.shape != (n,):
        raise ValueError(f"virtual_range_m must have {n} entries; got {r.shape}")
    lower, upper = edges[:-1], edges[1:]
    h_s = float(h_sat_m)

    t = np.zeros(n, dtype=float)
    f_h = np.zeros(n, dtype=float)
    theta = np.zeros(n, dtype=float)
    top = h_s
    for k in range(n):
        f_h[k], theta[k] = (float(np.ravel(v)[0]) for v in field(np.array([top])))
        for _ in range(2):
            row = magnetoionic_lamination_row(
                f[k], lower[: k + 1], upper[: k + 1], f_h[: k + 1], theta[: k + 1]
            )
            diagonal = float(row[k])
            if not np.isfinite(diagonal) or diagonal <= 0.0:
                raise ValueError(
                    f"lamina {k} has a non-positive or undefined diagonal coefficient "
                    f"{diagonal!r}; the frequency grid or the field model is degenerate"
                )
            accumulated = float(row[:k] @ t[:k]) if k else 0.0
            t[k] = (float(r[k]) - accumulated) / diagonal
            f_h[k], theta[k] = (
                float(np.ravel(v)[0]) for v in field(np.array([top - 0.5 * t[k]]))
            )
        top -= t[k]
    return t


def virtual_range_from_laminae_magnetoionic_m(
    thickness_m: np.ndarray,
    freq_hz: np.ndarray,
    f_start_hz: float,
    h_sat_m: float,
    field: FieldModel,
) -> np.ndarray:
    """Forward operator of the magnetoionic lamination.

    The group range each sounding frequency would show for a profile made of
    these laminae, with the field of each lamina taken at its mid height.  A
    round trip through this and :func:`lamina_thicknesses_magnetoionic_m`
    returns the thicknesses to the accuracy of the one field correction.
    """
    t = np.asarray(thickness_m, dtype=float)
    f, edges = _checked_grid(freq_hz, f_start_hz)
    if t.shape != (f.size,):
        raise ValueError(f"thickness_m must have {f.size} entries; got {t.shape}")
    lower, upper = edges[:-1], edges[1:]
    mid = float(h_sat_m) - (np.cumsum(t) - 0.5 * t)
    f_h, theta = field(mid)
    out = np.empty(f.size, dtype=float)
    for k in range(f.size):
        row = magnetoionic_lamination_row(
            f[k], lower[: k + 1], upper[: k + 1], f_h[: k + 1], theta[: k + 1]
        )
        out[k] = float(row @ t[: k + 1])
    return out


def lamination_matrix(freq_hz: np.ndarray, f_start_hz: float) -> np.ndarray:
    """Lower-triangular group-path coefficients ``G[k, j]``, dimensionless.

    ``G[k, j]`` is the group path a unit-thickness lamina ``j`` contributes to
    the sounding at ``f_k``, for ``j <= k``:

        ``G[k, j] = f_k (asin(f_j / f_k) - asin(f_{j-1} / f_k)) / (f_j - f_{j-1})``

    with ``f_0 = f_start_hz`` the plasma frequency at the spacecraft.  Entries
    with ``j > k`` are zero: a lamina deeper than the reflection level does not
    contribute, because the wave never reaches it.

    ``freq_hz`` must be strictly increasing and every entry must exceed
    ``f_start_hz``.  A frequency at or below the local plasma frequency at the
    spacecraft is evanescent and cannot have produced an echo, so its presence
    means the trace or the starting value is wrong, and this raises rather than
    quietly producing a coefficient from a complex arcsine.
    """
    f = np.asarray(freq_hz, dtype=float)
    f0 = float(f_start_hz)
    if f.ndim != 1 or f.size == 0:
        raise ValueError(f"freq_hz must be a non-empty 1-D array; got shape {f.shape}")
    if not np.all(np.diff(f) > 0.0):
        raise ValueError("freq_hz must be strictly increasing")
    if f[0] <= f0:
        raise ValueError(
            f"lowest sounding frequency {f[0]:.6g} Hz is at or below the plasma "
            f"frequency at the spacecraft {f0:.6g} Hz; that wave is evanescent and "
            "cannot have produced an echo"
        )

    edges = np.concatenate(([f0], f))          # f_0 .. f_n
    lower = edges[:-1]                         # f_{j-1}, j = 1..n
    upper = edges[1:]                          # f_j
    width = upper - lower
    if np.any(width <= 0.0):
        raise ValueError("lamina boundaries must be strictly increasing in frequency")

    # ratios[k, j] = f_j / f_k, clipped at 1 for the j == k diagonal where the
    # ratio is exactly 1 and floating point can land a few ulp above it.
    with np.errstate(invalid="ignore"):
        upper_ratio = np.clip(upper[None, :] / f[:, None], -1.0, 1.0)
        lower_ratio = np.clip(lower[None, :] / f[:, None], -1.0, 1.0)

    g = f[:, None] * (np.arcsin(upper_ratio) - np.arcsin(lower_ratio)) / width[None, :]
    return np.tril(g)


def virtual_range_from_laminae_m(
    thickness_m: np.ndarray,
    freq_hz: np.ndarray,
    f_start_hz: float,
) -> np.ndarray:
    """Forward operator: group range from a set of lamina thicknesses.

    ``R' = G t``.  This is the exact forward model of the lamination inversion,
    so a round trip through it and :func:`lamina_thicknesses_m` is exact to
    machine precision, and that round trip is what the unit test checks.  It is
    not the exact forward model of the ionosphere, which is a different claim
    and is what the comparison against a closed-form profile tests instead.
    """
    t = np.asarray(thickness_m, dtype=float)
    g = lamination_matrix(freq_hz, f_start_hz)
    if t.shape != (g.shape[0],):
        raise ValueError(f"thickness_m must have {g.shape[0]} entries; got {t.shape}")
    return g @ t


def lamina_thicknesses_m(
    virtual_range_m: np.ndarray,
    freq_hz: np.ndarray,
    f_start_hz: float,
) -> np.ndarray:
    """Solve the lower-triangular system for the lamina thicknesses, metres.

    Forward substitution, one lamina at a time.  No clamping: a non-positive
    thickness means the measured trace is not consistent with a monotonic
    profile, and that is information the caller needs rather than something to
    be tidied away.  :class:`LaminationRetrieval` counts them and flags the
    result.
    """
    r = np.asarray(virtual_range_m, dtype=float)
    g = lamination_matrix(freq_hz, f_start_hz)
    n = g.shape[0]
    if r.shape != (n,):
        raise ValueError(f"virtual_range_m must have {n} entries; got {r.shape}")

    t = np.zeros(n, dtype=float)
    for k in range(n):
        accumulated = float(g[k, :k] @ t[:k]) if k else 0.0
        diagonal = float(g[k, k])
        if diagonal <= 0.0:
            raise ValueError(
                f"lamina {k} has a non-positive diagonal coefficient {diagonal:.6g}; "
                "the frequency grid is degenerate"
            )
        t[k] = (float(r[k]) - accumulated) / diagonal
    return t


@register("l4", "lamination")
class LaminationRetrieval(RetrievalBase):
    """Laminar topside inversion, field-free or magnetoionic.

    Parameters
    ----------
    mode_index
        Which mode column of the ionogram to invert.  The baseline L3 tier
        produces one column and it is unlabelled; see
        :mod:`tss.l3_signal.ionogram`.
    start_plasma_hz
        Plasma frequency at the spacecraft, Hz.  When given it is treated as an
        in-situ measurement.  When ``None`` the lowest detected sounding
        frequency is used as a proxy, which is biased; see
        :mod:`tss.l4_retrieval.base`.
    include_spacecraft_point
        Whether to include the spacecraft altitude and its starting density as
        the top point of the returned profile.  On by default, because a
        profile that omits it hides the initial condition the rest of it was
        measured from.
    magnetoionic
        False inverts with the field-free group index.  True inverts the
        ordinary trace with the magnetoionic group index, which is what a scored
        profile needs.  The registered tier ``lamination_magnetoionic`` is this
        class with it set.
    geomagnetic_model
        Where the magnetoionic form gets its field: ``"auto"``, ``"igrf"`` or
        ``"dipole"``, passed to :func:`geomagnetic_field_model` at the
        platform's position and epoch.
    field
        A field model to use instead, as a callable from altitude to
        gyrofrequency and angle.  For tests and for a field known some other
        way.
    """

    def __init__(
        self,
        *,
        mode_index: int = 0,
        start_plasma_hz: float | None = None,
        include_spacecraft_point: bool = True,
        magnetoionic: bool = False,
        geomagnetic_model: str = "auto",
        field: FieldModel | None = None,
    ) -> None:
        super().__init__(mode_index=mode_index, start_plasma_hz=start_plasma_hz)
        self.include_spacecraft_point = bool(include_spacecraft_point)
        self.magnetoionic = bool(magnetoionic)
        if geomagnetic_model not in ("auto", "igrf", "dipole"):
            raise ValueError(
                f"unknown geomagnetic_model {geomagnetic_model!r}, expected auto, igrf or dipole"
            )
        self.geomagnetic_model = geomagnetic_model
        self.field = field

    def _settings(self) -> dict[str, Any]:
        return {
            "mode_index": self.mode_index,
            "start_plasma_hz": self.start_plasma_hz,
            "include_spacecraft_point": self.include_spacecraft_point,
            "magnetoionic": self.magnetoionic,
            "geomagnetic_model": self.geomagnetic_model if self.magnetoionic else None,
            "field_supplied": self.field is not None,
            "lamina_law": "plasma frequency linear in height",
            "status": "implemented",
        }

    def invert(
        self,
        ionogram: Ionogram,
        platform: PlatformState,
    ) -> DensityProfile:
        """Invert an ionogram to electron density against true height.

        Returns a new :class:`~tss.interfaces.DensityProfile` with altitude
        ascending.  ``ionogram`` is not modified.

        An ionogram with no detected trace returns an empty profile with
        ``flags["no_trace"]`` set, rather than raising.  A sounding that
        detected nothing is a result: it happens wherever the wave was
        evanescent across the whole band or penetrated across the whole band,
        and contract rule 3 in :mod:`tss.interfaces` says a tier reports it.
        """
        mode_label = str(tuple(ionogram.modes)[self.mode_index])
        if self.magnetoionic and mode_label.upper() == "X":
            raise ValueError(
                "the magnetoionic lamination inverts the ordinary trace; this column "
                "is labelled as the extraordinary mode, which reflects at a different "
                "density and needs the treatment of Budden (1961) section 12.17"
            )

        freq_hz, virtual_range_m, _ = trace_points(ionogram, mode_index=self.mode_index)
        h_s = float(platform.alt_m)

        if freq_hz.size == 0:
            return DensityProfile(
                alt_m=np.empty(0, dtype=float),
                ne=np.empty(0, dtype=float),
                flags={
                    "no_trace": True,
                    "reason": "the ionogram carries no detected trace point in this mode column",
                    "spacecraft_alt_m": h_s,
                    "mode": str(tuple(ionogram.modes)[self.mode_index]),
                    "magnetoionic": self.magnetoionic,
                },
                source=f"tss.l4_retrieval.lamination.{type(self).__name__}",
            )

        f0, start_source = starting_plasma_frequency_hz(
            ionogram,
            given_hz=self.start_plasma_hz,
            mode_index=self.mode_index,
        )

        # A frequency at or below the starting plasma frequency cannot have
        # produced an echo.  When the starting value came from the trace itself
        # the lowest point is exactly f0, so it is dropped; the profile then
        # starts at the second detected frequency and the first lamina runs
        # from the spacecraft to it.
        usable = freq_hz > f0
        dropped = int(freq_hz.size - np.count_nonzero(usable))
        freq_hz = freq_hz[usable]
        virtual_range_m = virtual_range_m[usable]

        if freq_hz.size == 0:
            return DensityProfile(
                alt_m=np.array([h_s], dtype=float),
                ne=density_from_trace(np.array([f0])),
                flags={
                    "no_trace": True,
                    "reason": (
                        "every detected frequency is at or below the assumed plasma "
                        "frequency at the spacecraft, so no lamina can be formed"
                    ),
                    "starting_value": start_source,
                    "start_plasma_hz": float(f0),
                    "spacecraft_alt_m": h_s,
                    "magnetoionic": self.magnetoionic,
                },
                source=f"tss.l4_retrieval.lamination.{type(self).__name__}",
            )

        field_record: dict[str, Any] = {}
        if self.magnetoionic:
            if self.field is not None:
                field = self.field
                field_record = {"model": "supplied by the caller"}
            else:
                field, field_record = geomagnetic_field_model(
                    platform.lat_deg,
                    platform.lon_deg,
                    platform.epoch,
                    h_s,
                    prefer=self.geomagnetic_model,
                )
            thickness_m = lamina_thicknesses_magnetoionic_m(
                virtual_range_m, freq_hz, f0, h_s, field
            )
            f_h_top, theta_top = field(np.array([h_s]))
            field_record["gyrofrequency_at_spacecraft_hz"] = float(np.ravel(f_h_top)[0])
            field_record["field_angle_from_vertical_deg"] = float(
                np.degrees(np.ravel(theta_top)[0])
            )
        else:
            thickness_m = lamina_thicknesses_m(virtual_range_m, freq_hz, f0)
        depth_m = np.cumsum(thickness_m)
        alt_m = h_s - depth_m
        ne = density_from_trace(freq_hz)

        if self.include_spacecraft_point:
            alt_m = np.concatenate(([h_s], alt_m))
            ne = np.concatenate((density_from_trace(np.array([f0])), ne))

        n_nonpositive = int(np.count_nonzero(thickness_m <= 0.0))
        order = np.argsort(alt_m)

        flags: dict[str, Any] = {
            "starting_value": start_source,
            "start_plasma_hz": float(f0),
            "spacecraft_alt_m": h_s,
            "mode": str(tuple(ionogram.modes)[self.mode_index]),
            "n_laminae": int(thickness_m.size),
            "n_nonpositive_laminae": n_nonpositive,
            "heights_monotonic": bool(n_nonpositive == 0),
            "dropped_below_start": dropped,
            "lowest_retrieved_alt_m": float(np.min(alt_m)),
            "lamina_law": "plasma frequency linear in height",
            "magnetoionic": self.magnetoionic,
            "magnetoionic_note": (
                "inverted with the magnetoionic ordinary-wave group index "
                "(Budden 1961, section 12.13)"
                if self.magnetoionic
                else "inverted with the no-field group refractive index, which places "
                "a profile tens of kilometres wrong where the field is steep; not "
                "scored against SCI-05 or SCI-07"
            ),
            "vertical_incidence_assumed": True,
        }
        if self.magnetoionic:
            flags["field_model"] = field_record
        if n_nonpositive:
            flags["monotonicity_warning"] = (
                f"{n_nonpositive} lamina thickness(es) came out non-positive, so the "
                "measured trace is not consistent with a monotonic profile.  The "
                "values are returned unclamped; see tss.l4_retrieval.depletion (WP10)."
            )

        return DensityProfile(
            alt_m=alt_m[order],
            ne=ne[order],
            ne_sigma=None,
            alt_sigma_m=None,
            flags=flags,
            source=f"tss.l4_retrieval.lamination.{type(self).__name__}",
        )


@register("l4", "lamination_magnetoionic")
class MagnetoionicLaminationRetrieval(LaminationRetrieval):
    """:class:`LaminationRetrieval` with the magnetoionic group index.

    This is the tier a scored campaign names.  It differs from the field-free
    tier only in its default, so the two cannot drift apart.
    """

    def __init__(self, *, magnetoionic: bool = True, **kwargs: Any) -> None:
        super().__init__(magnetoionic=magnetoionic, **kwargs)

