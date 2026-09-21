"""Layer: L0
WP: WP10
Claims: C2
Status: stub
Requirements: SCI-11

Equatorial plasma bubble depletions superposed on a background L0 state.  Not
implemented.  Goal level, WP10.

What the structure is
---------------------
An equatorial plasma bubble is a region of depleted plasma that rises out of the
post-sunset bottomside F region and extends upward along the geomagnetic field.
Four measured statistics set the shape, and they are the parameters this module
is written around rather than free knobs:

**Zonal envelope near 300 km.**  The across-field width of a depletion, which is
what a spacecraft crossing one actually samples.  It sets how long a crossing
lasts and therefore whether the sounder sees the wall, the interior or both
within one sweep.  At 7.5 km/s a 300 km envelope is a 40 second crossing, which
is about twelve sweeps at the 3.3 second baseline sweep period.

**Depletion depth of two to three orders of magnitude.**  The density inside a
mature bubble falls by a factor of 100 to 1000 against the ambient.  That is a
factor of 10 to 32 in plasma frequency, which is the number that matters here:
it moves the local plasma frequency at the spacecraft down by an order of
magnitude, and the bottom of the sounding band with it.  A bubble is therefore
the one common case where the band edge of INS-01 is set by something other than
the ambient topside.

**Field-aligned scale length near 550 km.**  Depletions are elongated along the
field by roughly this much, which is what makes them nearly two-dimensional in a
field-aligned frame and is what will let the generator work in that frame rather
than in three unconstrained dimensions.

**Eastward drift of 75 to 250 m/s.**  The structure moves, and at the fast end
its drift is a few per cent of the orbital velocity.  That matters for C3: the
Doppler decoupling assumes the reflecting surface is not itself moving fast, and
a 250 m/s drift with a component along the line of sight is a bias the
algorithm has to be shown to survive.

Why it is not implemented
-------------------------
Three decisions have to be made first and none of them is a detail.

1.  Whether the wall is a density discontinuity, a steep gradient with a scale
    length, or a statistical scattering surface.  The three give different
    off-vertical returns and therefore different answers to the question WP10
    exists to ask.  :mod:`tss.l1_propagation.offvertical` cannot be written
    until this is settled, and it is the consumer.
2.  Whether the generator perturbs a background state or replaces it.  Perturbing
    keeps the ambient topside from whichever L0 tier produced it, which is what a
    controlled comparison wants.  Replacing is easier and makes the comparison
    meaningless.
3.  What the truth is for scoring.  A retrieval that runs through a depletion
    produces a profile that is wrong in a way that looks plausible, so L5 has to
    be able to say which samples were inside the structure.  That means the
    generator has to return a mask as well as a state, and
    :class:`~tss.interfaces.IonosphereState` has no slot for one.  It would go in
    ``extra``, and no tier may depend on an ``extra`` key being present, so the
    scoring path needs designing rather than assuming.

Raises
------
NotImplementedError
    On every call.  WP10, ``docs/work-packages.md`` row WP10.  The interface
    below is the one the module will keep when it is written, so a campaign that
    names it can be validated now and will not have to be rewritten later.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from tss.interfaces import IonosphereState

__all__ = [
    "BubbleParameters",
    "bubble_depletion_factor",
    "perturb_state",
    "DEFAULT_ZONAL_ENVELOPE_M",
    "DEFAULT_FIELD_ALIGNED_SCALE_M",
    "DEFAULT_DEPLETION_DECADES",
    "DEFAULT_EASTWARD_DRIFT_MPS",
]

#: Zonal envelope near 300 km, metres.  The across-field width a crossing samples.
DEFAULT_ZONAL_ENVELOPE_M = 300.0e3

#: Field-aligned scale length, metres.
DEFAULT_FIELD_ALIGNED_SCALE_M = 550.0e3

#: Depletion depth as a range in decades of density.
DEFAULT_DEPLETION_DECADES = (2.0, 3.0)

#: Eastward drift range, metres per second.
DEFAULT_EASTWARD_DRIFT_MPS = (75.0, 250.0)

_MESSAGE = (
    "WP10: the equatorial plasma bubble generator is not implemented.  It is "
    "goal level, and three modelling decisions listed in the module docstring "
    "have to be made before it can be written: the wall representation, whether "
    "it perturbs or replaces the background state, and how L5 is told which "
    "samples were inside the structure.  See docs/work-packages.md row WP10."
)


class BubbleParameters:
    """The measured statistics one depletion is drawn from.  Not implemented.

    Parameters
    ----------
    zonal_envelope_m
        Across-field width near 300 km.  Default
        :data:`DEFAULT_ZONAL_ENVELOPE_M`.
    field_aligned_scale_m
        Elongation along the geomagnetic field.  Default
        :data:`DEFAULT_FIELD_ALIGNED_SCALE_M`.
    depletion_decades
        Depth of the depletion in decades of electron density, as a
        ``(low, high)`` range.  Default :data:`DEFAULT_DEPLETION_DECADES`, that
        is a factor of 100 to 1000.
    eastward_drift_mps
        Zonal drift speed range.  Default :data:`DEFAULT_EASTWARD_DRIFT_MPS`.
    wall_scale_m
        Gradient scale length across the wall.  Left without a default on
        purpose: this is decision 1 in the module docstring and putting a number
        on it here would settle by default a question that has not been settled.

    Raises
    ------
    NotImplementedError
        Always.  WP10.
    """

    def __init__(
        self,
        zonal_envelope_m: float = DEFAULT_ZONAL_ENVELOPE_M,
        field_aligned_scale_m: float = DEFAULT_FIELD_ALIGNED_SCALE_M,
        depletion_decades: Sequence[float] = DEFAULT_DEPLETION_DECADES,
        eastward_drift_mps: Sequence[float] = DEFAULT_EASTWARD_DRIFT_MPS,
        wall_scale_m: float | None = None,
    ) -> None:
        del (
            zonal_envelope_m,
            field_aligned_scale_m,
            depletion_decades,
            eastward_drift_mps,
            wall_scale_m,
        )
        raise NotImplementedError(_MESSAGE)


def bubble_depletion_factor(
    lat_deg: float,
    lon_deg: float,
    alt_m: np.ndarray,
    epoch: str,
    parameters: BubbleParameters,
    rng: np.random.Generator,
) -> np.ndarray:
    """Multiplicative depletion factor on the ambient density.  Not implemented.

    The intended return is an array parallel to ``alt_m`` with values in
    ``(0, 1]``, one meaning ambient and 1e-3 meaning a three-decade depletion, so
    that a caller multiplies rather than substitutes and the ambient topside from
    whichever L0 tier produced it survives the perturbation.

    Randomness comes from ``rng`` and from nowhere else, per contract rule 4 in
    :mod:`tss.interfaces`: the drift speed and the depth are drawn from the
    ranges in ``parameters``, so a run reproduces from its recorded seed.

    Raises
    ------
    NotImplementedError
        Always.  WP10.
    """
    del lat_deg, lon_deg, alt_m, epoch, parameters, rng
    raise NotImplementedError(_MESSAGE)


def perturb_state(
    state: IonosphereState,
    parameters: BubbleParameters,
    rng: np.random.Generator,
    *,
    crossing_offset_m: float = 0.0,
    extra: Mapping[str, Any] | None = None,
) -> IonosphereState:
    """A background state with one depletion superposed.  Not implemented.

    ``crossing_offset_m`` is the across-field distance from the centre of the
    depletion to the sub-satellite point, so that a campaign can sweep a whole
    crossing by stepping it rather than by moving the structure.

    The input state is not modified; a new
    :class:`~tss.interfaces.IonosphereState` is returned, per contract rule 1.

    Raises
    ------
    NotImplementedError
        Always.  WP10.
    """
    del state, parameters, rng, crossing_offset_m, extra
    raise NotImplementedError(_MESSAGE)
