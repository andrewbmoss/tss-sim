"""Layer: L4
WP: WP10
Claims: C2
Status: stub
Requirements: SCI-09

Inversion across a depletion crossing.  Not implemented.  Goal level, WP10.

Why the classical method fails here
------------------------------------
:mod:`tss.l4_retrieval.lamination` assumes two things that an equatorial plasma
bubble breaks at once.

It assumes the profile is monotonic, density increasing all the way down from
the spacecraft to the F2 peak.  Inside a depletion the density falls with
depth over some interval, and a sounding frequency that would have reflected
inside that interval passes through it instead and reflects deeper.  The
lamination has no lamina to put the interval in.  What it produces is a
thickness for the lamina that spans the depletion, and that thickness is wrong
in a way that displaces every height beneath it, because the inversion
accumulates downward.

It also assumes every echo came from directly below.  A depletion has steep
walls, and a wall returns energy from a direction that is not straight down, at
a shorter slant range than a vertical echo of the same delay.  On the ionogram
that appears as a second trace above the vertical one, and a trace extractor
that takes the first leading edge will take the wall return.  That geometry is
:mod:`tss.l1_propagation.offvertical` and it is a stub for the same reason.

The failure mode is the dangerous kind.  The inversion does not raise, does not
return NaN and does not produce anything obviously wrong.  It returns a smooth
plausible profile that is displaced, and the displacement is largest exactly
where the science is.

What would have to exist first
-------------------------------
1.  An L0 depletion model, planned as ``src/tss/l0_ionosphere/bubbles.py``.
    There is nothing to invert across until there is something to fly past.
2.  The off-vertical echo geometry of :mod:`tss.l1_propagation.offvertical`,
    which needs the ray-traced L1 tier of WP2.
3.  Arrival angle from the vector receive of :mod:`tss.l3_signal.vector`,
    WP11, if the detection of an off-vertical echo is to be made from the
    measurement rather than from a model.  If it is, then this work depends on
    C5 and it moves from goal level to conditional on another goal.

Open questions
--------------
* **Detect or model.**  Detecting a depletion crossing from the ionogram alone
  means recognising a signature in the trace, and the signature is not
  established.  Modelling one means assuming a bubble geometry and fitting it,
  which answers a different question.
* **Reject or ingest.**  Rejecting the affected frequencies loses the crossing,
  which is the science case for looking at bubbles at all.  Ingesting them
  corrupts the profile below.  A third option, inverting the profile above and
  below the crossing separately and leaving the crossing unretrieved with a
  stated gap, is probably the honest one, and it costs the continuity that the
  downward accumulation depends on.
* **What a retrieved profile through a depletion even means.**  A vertical
  profile is a fiction inside a structure whose horizontal scale is comparable
  to the sounding geometry.  The retrieval might have to return something other
  than a one-dimensional profile, and :class:`~tss.interfaces.DensityProfile`
  cannot represent that.  If it does, the interface changes, and that is worth
  knowing before WP10 rather than during it.

Raises
------
NotImplementedError
    On every call.  WP10.
"""

from __future__ import annotations

import numpy as np

from tss.interfaces import DensityProfile, Ionogram, PlatformState

__all__ = [
    "detect_depletion_crossing",
    "invert_across_depletion",
    "segment_trace",
]

_MESSAGE = (
    "WP10: inversion across a depletion crossing is not implemented.  It is goal "
    "level and it depends on the L0 depletion model, the off-vertical geometry of "
    "tss.l1_propagation.offvertical and, for a measured rather than modelled "
    "detection, the vector receive of tss.l3_signal.vector.  The classical "
    "lamination assumes a smooth monotonic profile and vertical incidence, and a "
    "depletion breaks both without raising.  See docs/work-packages.md row WP10."
)


def detect_depletion_crossing(
    ionogram: Ionogram,
    *,
    aoa_threshold_deg: float | None = None,
) -> np.ndarray:
    """Boolean mask of trace points affected by a depletion.  Not implemented.

    Intended to use arrival angle where the vector receive supplies it and a
    trace-shape signature where it does not.  The signature is not established,
    which is open question one in the module docstring.

    Raises
    ------
    NotImplementedError
        Always.  WP10.
    """
    del ionogram, aoa_threshold_deg
    raise NotImplementedError(_MESSAGE)


def segment_trace(
    ionogram: Ionogram,
    affected: np.ndarray,
) -> list[tuple[int, int]]:
    """Split a trace into segments the classical method can invert.  Not implemented.

    Intended to return index ranges either side of a crossing so that each can
    be inverted separately with its own starting value, leaving the crossing
    unretrieved and the gap stated.  The cost is that the segment below the gap
    has no reliable starting value, which is the continuity the downward
    accumulation depends on.

    Raises
    ------
    NotImplementedError
        Always.  WP10.
    """
    del ionogram, affected
    raise NotImplementedError(_MESSAGE)


def invert_across_depletion(
    ionogram: Ionogram,
    platform: PlatformState,
) -> DensityProfile:
    """Retrieve a profile through a depletion crossing.  Not implemented.

    Raises
    ------
    NotImplementedError
        Always.  WP10.
    """
    del ionogram, platform
    raise NotImplementedError(_MESSAGE)
