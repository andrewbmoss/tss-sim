"""Layer: L1
WP: WP2
Claims: C1, C2, C3
Status: implemented
Requirements: INS-01, INS-10

The analytic L1 tier: vertical incidence, no ray bending.

This is the fast tier.  It assumes a plane-stratified ionosphere and a ray that
goes straight up and straight back down, which for a topside sounder looking at
a horizontally uniform F region is a good assumption and is exactly the case
that V2 checks against the ray-traced tier of
:mod:`tss.l1_propagation.raytrace`.  It is fast enough to sweep a trade study,
and that is what it is for.

What it computes, per sounding frequency and per magnetoionic mode:

**Reflection altitude** by bisection on the mode's cutoff.  For the O mode that
is the plasma-frequency crossing; for the X mode it is the crossing of
:func:`tss.constants.x_mode_cutoff`.

**Group delay** by integrating the group refractive index from the reflection
altitude up to the spacecraft and doubling.  In the field-free, collisionless
limit the phase index is ``n = sqrt(1 - X)`` and the group index is its
reciprocal, so the virtual range is

    P' = integral from h_r to h_sat of  dh / sqrt(1 - X(h))

which diverges as an inverse square root at the reflection point and is
therefore integrated after the substitution ``h = h_r + (h_sat - h_r) s^2``.
That substitution turns the integrand into a bounded function of ``s``, because
near reflection ``n`` grows as ``s`` while the Jacobian ``2 (h_sat - h_r) s``
falls as ``s``, and the two cancel.  The midpoint rule is then used in ``s`` so
that the endpoint is never evaluated at all.  This is the standard treatment
and it converges as the square of the step count; the alternative of stopping
the integral short of reflection converges as the square root of wherever it
was stopped, which is to say not usefully.

**Amplitude** from free-space spreading plus the absorption integral.  A smooth
overdense layer reflects specularly and does not focus, so the geometric term
is one-way free-space spreading over a path of twice the equivalent range,
``20 log10(lambda / (4 pi * 2 P'))``.  The equivalent range used is the group
path by default, because that is the path the wavefront has actually travelled.
The author's link budget module uses the geometric range for the same term, so
the two differ by ``20 log10(P'/R)``, which is a fraction of a decibel except
close to the critical frequency where the retardation runs away.  Pass
``spreading_range="reflection"`` to reproduce the link budget exactly.  That
seam is deliberate and it is recorded here rather than left for a reader to
discover from a discrepancy.

**Magnetoionic splitting** into O and X.  Honest statement of what this costs:
the reflection altitudes are exact algebra, but the group delays are both
computed with the field-free group index, which is exact for the O mode only in
the limit of a negligible field and is an approximation for the X mode near its
cutoff, where the magnetoionic group index genuinely differs.  Every X row
carries the flag ``"x_group_delay_field_free"`` so that no downstream result can
depend on that approximation without saying so, and
:func:`tss.l1_propagation.magnetoionic.appleton_hartree` is the stub that will
replace it.  A campaign that cares about the X trace shape should use the
ray-traced tier.

**Doppler** from :mod:`tss.l1_propagation.doppler`, which is the forward model
C3 must invert.  The tier reports the two-way shift in ``EchoSet.doppler_hz``,
carrying the field-free refractive index at the spacecraft, and leaves the apparent-range bias to the layers that know the sweep rate,
because :class:`~tss.interfaces.EchoSet` has no sweep to refer to.

Band edges are results, not errors.  Below the local plasma frequency at the
spacecraft the row is flagged ``"evanescent"``; above the critical frequency it
is flagged ``"penetrating"``.  In both cases the measurable quantities are NaN
and nothing raises.  Those two flags are the band-edge study of INS-01 and they
are the reason contract rule 3 in :mod:`tss.interfaces` exists.

The tier is deterministic.  It accepts the ``rng`` the protocol requires and
does not draw from it.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np

from tss.constants import C_LIGHT
from tss.interfaces import EchoSet, IonosphereState, PlatformState
from tss.l1_propagation import magnetoionic
from tss.l1_propagation.absorption import absorption_db, reflection_altitude_m
from tss.l1_propagation.base import PropagationBase, plasma_frequency_at
from tss.l1_propagation.doppler import (
    closing_velocity_mps,
    doppler_shift_hz,
    refractive_index_at_platform,
    vertical_line_of_sight,
)
from tss.registry import register

__all__ = ["AnalyticPropagation", "group_range_m"]

#: Flag put on every X-mode row, recording that its group delay was computed
#: with the field-free group index rather than the Appleton-Hartree one.
X_GROUP_FLAG = "x_group_delay_field_free"

#: Flag put on an X-mode row when the state carries no magnetic field, so the
#: X cutoff cannot be evaluated at all.
X_NO_FIELD_FLAG = "x_mode_no_field"


def group_range_m(
    state: IonosphereState,
    f_hz: float,
    h_refl_m: float,
    h_sat_m: float,
    *,
    n_steps: int = 2000,
    interp: str = "linear",
) -> float:
    """Virtual range from the reflection altitude to the spacecraft, metres.

    The integral of the field-free group refractive index ``1 / sqrt(1 - X)``
    over the vertical path, evaluated with the substitution
    ``h = h_refl + (h_sat - h_refl) s^2`` and the midpoint rule in ``s``.

    Parameters
    ----------
    n_steps
        Number of midpoint cells.  The error falls as ``n_steps^-2``.  2000 is
        comfortably converged for the profiles in use and costs nothing.

    Notes
    -----
    ``1 - X`` is clipped at a small positive floor before the square root.  On a
    strictly monotonic topside that clip never bites, because the reflection
    altitude was located where ``X = 1`` and ``X`` falls above it.  It bites
    only when interpolation between grid nodes puts a sample marginally inside
    the reflecting region, and clipping there is better than a NaN propagating
    through the whole row.  A profile where the clip bites often is a profile
    whose altitude grid is too coarse for the sounding frequency, and that is an
    L0 resolution problem rather than something this function can fix.
    """
    h_refl = float(h_refl_m)
    h_sat = float(h_sat_m)
    if h_sat <= h_refl:
        return 0.0

    n = int(n_steps)
    s = (np.arange(n, dtype=float) + 0.5) / n
    h = h_refl + (h_sat - h_refl) * s**2
    f_p = plasma_frequency_at(state, h, interp=interp)
    x = (f_p / float(f_hz)) ** 2
    n_phase = np.sqrt(np.clip(1.0 - x, 1e-12, None))
    integrand = (1.0 / n_phase) * 2.0 * (h_sat - h_refl) * s
    return float(np.sum(integrand) / n)


@register("l1", "analytic")
class AnalyticPropagation(PropagationBase):
    """Vertical-incidence analytic propagation tier.

    Parameters
    ----------
    modes
        Magnetoionic modes to emit, in the order the rows appear.  The default
        emits both.  ``("O",)`` gives an O-only echo set, which is what a
        campaign that has no field model should ask for.
    n_absorption_steps
        Samples on the absorption integration grid.  Passed through to
        :func:`tss.l1_propagation.absorption.absorption_db`.
    n_group_steps
        Midpoint cells on the group-delay integral.
    spreading_range
        ``"group"`` uses the virtual range for the free-space spreading term,
        which is the physically correct equivalent path.  ``"reflection"`` uses
        the geometric range to the reflection point, which reproduces the
        author's link budget module.
    default_b_tesla
        Field magnitude to fall back on when the state carries none.  Leaving it
        ``None`` means an X-mode row on a fieldless state is masked and flagged
        rather than computed against an invented field.
    interp
        Profile interpolation, ``"linear"`` or ``"log"``.  See
        :func:`tss.l1_propagation.base.density_at`.
    """

    def __init__(
        self,
        *,
        modes: Sequence[str] = ("O", "X"),
        n_absorption_steps: int = 3000,
        n_group_steps: int = 2000,
        spreading_range: str = "group",
        default_b_tesla: float | None = None,
        interp: str = "linear",
    ) -> None:
        super().__init__(interp=interp)
        modes_upper = tuple(str(m).upper() for m in modes)
        for m in modes_upper:
            if m not in ("O", "X"):
                raise ValueError(f"unknown mode {m!r}; this tier emits 'O' and 'X'")
        if spreading_range not in ("group", "reflection"):
            raise ValueError(
                f"unknown spreading_range {spreading_range!r}, "
                "expected 'group' or 'reflection'"
            )
        self.modes = modes_upper
        self.n_absorption_steps = int(n_absorption_steps)
        self.n_group_steps = int(n_group_steps)
        self.spreading_range = spreading_range
        self.default_b_tesla = default_b_tesla

    # -- manifest -----------------------------------------------------------

    def _settings(self) -> dict[str, Any]:
        return {
            "interp": self.interp,
            "modes": list(self.modes),
            "n_absorption_steps": self.n_absorption_steps,
            "n_group_steps": self.n_group_steps,
            "spreading_range": self.spreading_range,
            "default_b_tesla": self.default_b_tesla,
            "geometry": "vertical incidence, no ray bending",
            "group_index": "field-free, collisionless",
            "doppler_index": "field-free refractive index at the platform",
        }

    # -- the protocol -------------------------------------------------------

    def echoes(
        self,
        state: IonosphereState,
        platform: PlatformState,
        freq_hz: np.ndarray,
        rng: np.random.Generator,
    ) -> EchoSet:
        """Echoes for a set of sounding frequencies.  Satisfies ``L1Propagation``.

        ``rng`` is accepted because the protocol requires it and is not drawn
        from: this tier is deterministic and two calls on the same inputs give
        identical output whatever generator is passed.
        """
        del rng  # deterministic tier

        freqs = np.atleast_1d(np.asarray(freq_hz, dtype=float))
        h_sat = float(platform.alt_m)
        v_close = self._closing_velocity(platform)
        f_n_sat = float(self.plasma_frequency_at(state, h_sat))

        n_rows = freqs.size * len(self.modes)
        out_freq = np.empty(n_rows, dtype=float)
        out_mode = np.empty(n_rows, dtype=object)
        out_delay = np.full(n_rows, np.nan, dtype=float)
        out_amp = np.full(n_rows, np.nan, dtype=float)
        out_abs = np.full(n_rows, np.nan, dtype=float)
        out_dopp = np.full(n_rows, np.nan, dtype=float)
        out_alt = np.full(n_rows, np.nan, dtype=float)
        out_flags = np.empty(n_rows, dtype=object)

        row = 0
        for mode in self.modes:
            for f in freqs:
                f = float(f)
                out_freq[row] = f
                out_mode[row] = mode
                # NaN below the plasma frequency at the spacecraft, where the
                # row is evanescent and there is no echo to carry a shift.
                mu_sat = float(refractive_index_at_platform(f, f_n_sat))
                out_dopp[row] = float(doppler_shift_hz(v_close, f, mu_sat=mu_sat))

                h_refl, flag = self._reflection(state, f, h_sat, mode)
                if isinstance(h_refl, str):
                    out_flags[row] = h_refl
                    row += 1
                    continue
                if flag == X_NO_FIELD_FLAG:
                    out_flags[row] = X_NO_FIELD_FLAG
                    row += 1
                    continue

                out_alt[row] = h_refl
                p_group = group_range_m(
                    state,
                    f,
                    h_refl,
                    h_sat,
                    n_steps=self.n_group_steps,
                    interp=self.interp,
                )
                out_delay[row] = 2.0 * p_group / C_LIGHT

                a_db = absorption_db(
                    state,
                    f,
                    h_sat,
                    h_refl_m=h_refl,
                    n_steps=self.n_absorption_steps,
                    interp=self.interp,
                )
                out_abs[row] = a_db

                r_spread = p_group if self.spreading_range == "group" else (h_sat - h_refl)
                out_amp[row] = self._amplitude_db(f, r_spread, a_db)
                out_flags[row] = flag
                row += 1

        return EchoSet(
            freq_hz=out_freq,
            mode=out_mode,
            group_delay_s=out_delay,
            amplitude_db=out_amp,
            absorption_db=out_abs,
            doppler_hz=out_dopp,
            aoa_deg=None,
            reflection_alt_m=out_alt,
            flags=out_flags,
            source=(
                f"tss.l1_propagation.analytic (modes={'+'.join(self.modes)}, "
                f"spreading={self.spreading_range})"
            ),
        )

    # -- internals ----------------------------------------------------------

    def _closing_velocity(self, platform: PlatformState) -> float:
        """Line-of-sight closing velocity for a vertical sounding, m/s.

        Zero when the platform carries no velocity, which is the right answer
        for a stationary sounder rather than a masked one.  The line of sight is
        local down, so only the vertical component of the platform velocity
        contributes; see the frame note in
        :func:`tss.l1_propagation.doppler.doppler_geometry`.
        """
        if platform.velocity_ecef_mps is None:
            return 0.0
        return closing_velocity_mps(
            platform.velocity_ecef_mps, vertical_line_of_sight()
        )

    def _reflection(
        self,
        state: IonosphereState,
        f_hz: float,
        h_sat_m: float,
        mode: str,
    ) -> tuple[float | str, str]:
        """Reflection altitude and the flag the row should carry.

        Returns the sentinel string in place of an altitude for the evanescent
        and penetrating cases, and a second string that is the row flag: empty
        for a clean O row, :data:`X_GROUP_FLAG` for a clean X row, and
        :data:`X_NO_FIELD_FLAG` when the X cutoff could not be evaluated.
        """
        if mode == "O":
            h = reflection_altitude_m(state, f_hz, h_sat_m, interp=self.interp)
            return h, ""

        try:
            h = magnetoionic.reflection_altitude_m(
                state,
                f_hz,
                h_sat_m,
                "X",
                interp=self.interp,
                default_b_tesla=self.default_b_tesla,
            )
        except ValueError:
            # No field and no default.  Mask the row rather than invent one.
            return float("nan"), X_NO_FIELD_FLAG
        return h, X_GROUP_FLAG

    @staticmethod
    def _amplitude_db(f_hz: float, range_m: float, absorption_loss_db: float) -> float:
        """Two-way path gain in dB relative to the transmitted EIRP.

        Free-space spreading over a path of twice ``range_m``, less the
        absorption.  Negative, because it is a loss expressed as a gain.
        ``absorption_loss_db`` is the positive two-way loss that
        :func:`tss.l1_propagation.absorption.absorption_db` returns.
        """
        if not np.isfinite(range_m) or range_m <= 0.0:
            return float("nan")
        lam = C_LIGHT / f_hz
        spreading_db = 20.0 * math.log10(lam / (4.0 * math.pi * 2.0 * range_m))
        return float(spreading_db - absorption_loss_db)
