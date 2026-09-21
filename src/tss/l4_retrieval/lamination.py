"""Layer: L4
WP: WP6
Claims: C2
Status: partial
Requirements: INS-05

Virtual height to true height inversion by lamination.

The classical case is implemented and it is implemented properly: a smooth,
monotonic topside profile inverted downward from the spacecraft, with the
plasma frequency taken as linear in height inside each lamina, which gives a
closed-form group-path integral and a lower-triangular system that solves by
forward substitution.  This is the method Jackson used for the Alouette and
ISIS topside ionograms and it is the standard against which anything else has
to be measured.

What is not implemented is the magnetoionic correction.  See "Limits" below.
That is why the status is partial rather than implemented, and it is stated
here rather than in a footnote because the size of the correction is not small.

The geometry
------------
A topside sounder looks down.  The spacecraft sits at altitude ``h_s`` in
plasma of local plasma frequency ``f_0``.  A sounding at frequency ``f`` passes
downward through increasing density until it reaches the level where ``f_N = f``
and reflects.  What the instrument measures is the group range

    R'(f) = integral from h_r(f) to h_s of mu'(f, h) dh

where ``mu' = 1 / sqrt(1 - f_N(h)^2 / f^2)`` is the group refractive index in
the absence of a magnetic field.  The integrand is 1 at the spacecraft and
diverges at the reflection point, which is why the integral needs a closed form
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

Limits
------
**The magnetic field is ignored.**  The group refractive index used is the
no-field form.  In the real topside the ordinary and extraordinary waves have
different group indices and reflect at different levels, and the correction to
the O-mode group path from the field is not negligible at the low end of the
band, where the gyrofrequency is a substantial fraction of the sounding
frequency.  The proper treatment uses the magnetoionic group index and needs
the field direction along the ray, which means it needs
:mod:`tss.l1_propagation.magnetoionic` to expose a group index and the L0 state
to carry ``b_field_t``.  Until then this tier inverts as though the field were
absent and records ``flags["magnetoionic"] = False`` so that no result from it
is read as a field-corrected retrieval.  The correction is WP6 and it is the
single largest known bias in this module.

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
    "lamination_matrix",
    "lamina_thicknesses_m",
    "virtual_range_from_laminae_m",
]

_MAGNETOIONIC_MESSAGE = (
    "WP6: the magnetoionic lamination is not implemented.  This tier inverts with "
    "the no-field group refractive index and records flags['magnetoionic'] = False.  "
    "The field-corrected inversion needs a magnetoionic group index from "
    "tss.l1_propagation.magnetoionic and the field direction along the ray from the "
    "L0 state, neither of which is available yet.  See the 'Limits' section of the "
    "module docstring of tss.l4_retrieval.lamination and docs/work-packages.md row WP6."
)


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
    """Classical laminar topside inversion, no-field group index.

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
        Reserved.  Setting it True raises :class:`NotImplementedError` naming
        WP6.  The field-corrected inversion is the largest known bias in this
        module and it is not written.
    """

    def __init__(
        self,
        *,
        mode_index: int = 0,
        start_plasma_hz: float | None = None,
        include_spacecraft_point: bool = True,
        magnetoionic: bool = False,
    ) -> None:
        super().__init__(mode_index=mode_index, start_plasma_hz=start_plasma_hz)
        self.include_spacecraft_point = bool(include_spacecraft_point)
        self.magnetoionic = bool(magnetoionic)

    def _settings(self) -> dict[str, Any]:
        return {
            "mode_index": self.mode_index,
            "start_plasma_hz": self.start_plasma_hz,
            "include_spacecraft_point": self.include_spacecraft_point,
            "magnetoionic": self.magnetoionic,
            "lamina_law": "plasma frequency linear in height",
            "status": "partial",
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
        if self.magnetoionic:
            raise NotImplementedError(_MAGNETOIONIC_MESSAGE)

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
                    "magnetoionic": False,
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
                    "magnetoionic": False,
                },
                source=f"tss.l4_retrieval.lamination.{type(self).__name__}",
            )

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
            "magnetoionic": False,
            "magnetoionic_note": (
                "inverted with the no-field group refractive index; the field "
                "correction is WP6 and is the largest known bias in this tier"
            ),
            "vertical_incidence_assumed": True,
        }
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
