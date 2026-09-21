"""Layer: L3
WP: WP5
Claims: C2, C3
Status: stub
Requirements: SCI-09

Along-track aperture synthesis.  Not implemented.  Exploratory extension, WP5.

The idea
--------
The spacecraft moves 7.6 to 15 km during the 1 to 2 s coherent interval INS-15
asks for.  That motion is a synthetic aperture whether or not it is used as
one.  If the returns from successive sweeps stay coherent over that path, the
phase history across the path carries along-track information, and processing
it coherently would give an along-track resolution finer than the footprint of
the electrically small antenna, which has essentially no directivity at all.

The connection to the requirement is direct.  SCI-09 asks for 50 km along-track
sampling at threshold and 25 km at goal, and that number comes from the
horizontal structure scales seen by GOLD, COSMIC-2 and Swarm.  Sampling at
50 km with a footprint much larger than 50 km does not resolve 50 km structure;
it averages over it and reports the average at 50 km intervals.  Aperture
synthesis is the only mechanism available on this platform that would make the
resolution match the sampling.

What is already known and what is not
---------------------------------------
Ganguly and colleagues observed that the orbital motion of a topside sounder is
in principle sufficient to synthesise a useful aperture, and did not develop
the processing to do it.  That is the state of the published art as far as the
author has found it: the observation exists, the implementation does not.  This
module is where the implementation would go, and it is listed as exploratory
rather than as a claim because the questions below are not answered and a
dissertation that promised aperture synthesis on top of C3 would be promising
two unstarted things.

Open questions
--------------
1.  **Coherence.**  Synthetic aperture processing needs the scene to stay
    coherent across the aperture.  The reflecting surface here is a plasma
    layer that is moving, evolving and turbulent, not a stationary rough
    surface.  Whether the topside echo is coherent over 7.6 to 15 km of
    along-track path is unmeasured, and nothing in the Alouette and ISIS
    archive bears on it, because a pulse sounder never integrated across
    sweeps.  If the coherence length is much shorter than the aperture then
    nothing else in this module matters.
2.  **What the scatterer is.**  Synthetic aperture radar images point
    scatterers on a surface.  A topside sounding reflects from an extended,
    smoothly varying refractive index profile, and the specular return from a
    stratified medium does not have the along-track phase history a point
    scatterer would.  The correct signal model is not the usual one and
    borrowing the usual one would produce an image of something that is not
    there.
3.  **Interaction with C3.**  The along-track phase history and the Doppler
    history are the same measurement seen two ways.  Using the motion to
    separate Doppler from range, which is C3, and using the same motion to
    synthesise an aperture are not independent uses of the same information.
    Whether they can both be had, or whether one spends the coherence the other
    needs, is unresolved and is the first thing to settle.
4.  **Geometry.**  The reflection point for a vertical sounding is below the
    spacecraft and moves with it, so the along-track baseline does not view a
    fixed point from different angles the way a side-looking radar does.  It is
    closer to a forward-looking geometry with very poor along-track angular
    diversity, which is the geometry synthetic aperture processing handles
    worst.

Raises
------
NotImplementedError
    On every call.  WP5, exploratory.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from tss.interfaces import Ionogram, PlatformState, ReceivedSignal

__all__ = [
    "along_track_aperture_m",
    "synthesize_aperture",
    "along_track_resolution_m",
]

_MESSAGE = (
    "WP5: along-track aperture synthesis is not implemented.  It is an "
    "exploratory extension, not a claim.  Ganguly et al. observed that the "
    "orbital motion suffices in principle to synthesise a useful aperture but did "
    "not develop the processing, and the open questions in the module docstring "
    "of tss.l3_signal.novel.aperture_synthesis, starting with whether the topside "
    "echo is coherent over 7.6 to 15 km of along-track path, are unanswered.  "
    "See docs/work-packages.md row WP5."
)


def along_track_aperture_m(
    ground_speed_mps: float | np.ndarray,
    interval_s: float | np.ndarray,
) -> np.ndarray:
    """Length of the aperture the motion sweeps out.  Not implemented.

    The arithmetic is a product.  It is left unimplemented because a function
    that returns an aperture length invites a resolution to be quoted from it,
    and the resolution such a formula would give assumes a coherence that has
    not been demonstrated.  Quoting an along-track resolution from an aperture
    length alone is the specific overstatement this module exists to avoid.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del ground_speed_mps, interval_s
    raise NotImplementedError(_MESSAGE)


def synthesize_aperture(
    signals: Sequence[ReceivedSignal],
    platform: Sequence[PlatformState],
) -> Ionogram:
    """Coherently combine successive sweeps along track.  Not implemented.

    Intended to return an ionogram whose along-track cell is set by the
    synthetic aperture rather than by the antenna footprint, with the achieved
    cell recorded in ``flags`` so that SCI-09 is scored against a measured
    resolution and not against the sampling interval.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del signals, platform
    raise NotImplementedError(_MESSAGE)


def along_track_resolution_m(
    aperture_m: float | np.ndarray,
    f_hz: float | np.ndarray,
    range_m: float | np.ndarray,
    *,
    coherence_length_m: float | None = None,
) -> np.ndarray:
    """Achievable along-track resolution.  Not implemented.

    ``coherence_length_m`` is in the signature because the answer is bounded by
    the smaller of the aperture and the coherence length, and a function that
    did not take it would be quoting the aperture-limited answer as though the
    coherence were unlimited.  There is no measured value for it.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del aperture_m, f_hz, range_m, coherence_length_m
    raise NotImplementedError(_MESSAGE)
