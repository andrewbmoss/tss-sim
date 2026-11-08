"""Layer: L1
WP: WP4
Claims: C1, C2
Status: implemented
Requirements: SCI-05, SCI-07

The stratified group-path model: the virtual range a laminated profile
returns, field-free and with the magnetoionic ordinary-wave group index
(Budden 1961, section 12.13), and the forward substitution that inverts it
for the lamina thicknesses; and the field model along the vertical below a
spacecraft.  The theory is set out in :mod:`tss.l4_retrieval.lamination`.

These are propagation physics, used by L4 (the retrieval tiers, which add
the policy: the starting value, the monotonicity record, the flags) and by
L3 (the trace tracker, which labels O and X by the self-consistency of an
inversion, and the field model it uses).  They lived in
:mod:`tss.l4_retrieval.lamination` until 2026-11-08, which made the two
L3 uses imports from L3 into L4 (code review F4); that module re-exports
them.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

__all__ = [
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
