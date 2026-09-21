"""Layer: L1
WP: WP10
Claims: C2
Status: stub
Requirements: -

Off-vertical echo geometry for equatorial plasma depletion walls.  Not
implemented.  Goal level, WP10.

The problem
-----------
An equatorial plasma bubble is a region of depleted density that extends along
the magnetic field and has steep walls across it.  A sounder flying past one
sees echoes that do not come from straight below.  The wall acts as a near
specular reflector when the ray meets it close to normal incidence, and as a
scatterer when it does not, and the returns arrive at a shorter slant range than
the vertical echo would for the same delay, so on the ionogram they appear as a
second trace above the vertical one.

That is the whole difficulty for retrieval.  The lamination inversion of
:mod:`tss.l4_retrieval.lamination` assumes every echo came from directly below
and that the profile it is inverting is smooth and monotonic.  An off-vertical
echo violates the first assumption and a depletion violates the second.  A
retrieval that ingests both without noticing produces a profile that is wrong
in a way that looks plausible, which is the worst failure mode available.

Approach
--------
The geometry is a ray from the platform to a point on the wall, with the wall
represented as a surface rather than as a stratified medium.  Three pieces have
to exist first:

1.  An L0 depletion model, planned as ``src/tss/l0_ionosphere/bubbles.py``.
    Without a wall there is nothing to reflect from.  Whether the wall is
    modelled as a density discontinuity, as a steep gradient with a scale
    length, or as a statistical scattering surface changes the answer and is
    not a detail.
2.  The ray-traced tier of :mod:`tss.l1_propagation.raytrace`, because finding
    the ray that meets the wall at the right angle is a search over a fan and
    not a closed-form calculation.
3.  A decision about what an off-vertical echo means in
    :class:`~tss.interfaces.EchoSet`.  The interface already carries
    ``aoa_deg``, which is the right place for it, but nothing populates that
    field yet and nothing downstream reads it.

Open questions
--------------
* Whether the returns are strong enough to matter at a sub-watt transmit power.
  A wall return is not a specular reflection from an overdense layer and does
  not get the free-space image geometry the vertical echo gets; if it is 30 dB
  down it is a curiosity rather than a problem.
* Whether an off-vertical echo is separable from the vertical one without the
  vector receive of WP11.  If separation needs arrival angle, then this work
  depends on C5 and it moves from goal level to conditional on another goal.
* What the retrieval should do when it detects one.  Rejecting the echo loses
  the depletion crossing entirely, which is the science case for looking at
  them; ingesting it corrupts the profile.  That is the question
  :mod:`tss.l4_retrieval.depletion` exists to answer, and neither module can
  answer it alone.

Raises
------
NotImplementedError
    On every call.  WP10, ``docs/work-packages.md`` row WP10 and
    ``docs/architecture.md`` section "Known seams", subsection "The analytic L1
    tier cannot produce off-vertical echoes".
"""

from __future__ import annotations

from typing import Any

import numpy as np

from tss.interfaces import EchoSet, IonosphereState, PlatformState

__all__ = ["wall_echo_geometry", "off_vertical_echoes"]

_MESSAGE = (
    "WP10: off-vertical echo geometry for depletion walls is not implemented.  "
    "It is goal level and it depends on the L0 depletion model and on the "
    "ray-traced L1 tier, neither of which exists.  See docs/work-packages.md "
    "row WP10 and docs/architecture.md section 'Known seams'."
)


def wall_echo_geometry(
    platform: PlatformState,
    wall_normal: np.ndarray,
    wall_point_m: np.ndarray,
    f_hz: float | np.ndarray,
) -> dict[str, Any]:
    """Slant range, arrival angle and incidence angle for a wall return.

    Not implemented.  The intended return is a mapping carrying the slant range
    to the specular point, the arrival angle in the platform frame, and the
    angle of incidence on the wall, which is what decides whether the return is
    specular at all.

    Raises
    ------
    NotImplementedError
        Always.  WP10.
    """
    del platform, wall_normal, wall_point_m, f_hz
    raise NotImplementedError(_MESSAGE)


def off_vertical_echoes(
    state: IonosphereState,
    platform: PlatformState,
    freq_hz: np.ndarray,
    rng: np.random.Generator,
) -> EchoSet:
    """Echo set containing the off-vertical returns only.  Not implemented.

    Intended to be added to the vertical echo set rather than to replace it, so
    that a study can run with and without depletion returns and compare the two
    retrievals.  That comparison is the point of WP10.

    Raises
    ------
    NotImplementedError
        Always.  WP10.
    """
    del state, platform, freq_hz, rng
    raise NotImplementedError(_MESSAGE)
