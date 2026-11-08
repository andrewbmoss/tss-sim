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

from typing import Any

import numpy as np

from tss.interfaces import DensityProfile, Ionogram, PlatformState
from tss.l1_propagation.lamination import (  # the kernels moved to L1 (code review F4)
    FieldModel,
    constant_field,
    geomagnetic_field_model,
    lamina_thicknesses_m,
    lamina_thicknesses_magnetoionic_m,
    lamination_matrix,
    magnetoionic_lamination_row,
    virtual_range_from_laminae_m,
    virtual_range_from_laminae_magnetoionic_m,
)
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
    start_policy
        ``"onboard"`` (the default since WP4 Phase 6) or ``"lowest_trace"``
        (the WP4 Phase 4 behaviour).  ``"onboard"``: when
        no value is given, use the tracker's estimate of the plasma frequency
        at the spacecraft (``ionogram.flags["tracker"]["priors"]["f_n_hz"]``)
        where it lies below the lowest detected frequency (WP4 Phase 6).
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
        start_policy: str = "onboard",
    ) -> None:
        super().__init__(mode_index=mode_index, start_plasma_hz=start_plasma_hz)
        if start_policy not in ("lowest_trace", "onboard"):
            raise ValueError(f"unknown start_policy {start_policy!r}, expected lowest_trace or onboard")
        self.start_policy = start_policy
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
            "start_policy": self.start_policy,
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

        onboard = None
        if self.start_policy == "onboard":
            pri = ((ionogram.flags or {}).get("tracker") or {}).get("priors") or {}
            v = pri.get("f_n_hz")
            onboard = float(v) if isinstance(v, (int, float)) else None
        f0, start_source = starting_plasma_frequency_hz(
            ionogram,
            given_hz=self.start_plasma_hz,
            mode_index=self.mode_index,
            onboard_hz=onboard,
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

