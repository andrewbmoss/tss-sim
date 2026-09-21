"""Layer: L1
WP: WP2
Claims: C1, C2
Status: stub
Requirements: INS-01

The ray-traced L1 tier.  Not implemented.

Approach
--------
Three-dimensional numerical ray tracing is delegated to PHaRLAP 4.7.4 through
the file bridge described in ``docs/matlab-python-bridge.md``.  PHaRLAP is a
MATLAB toolkit with no Python binding, and reimplementing Haselgrove ray tracing
in Python purely to avoid a language boundary would produce a second thing to
verify and no new science.  So the rays are generated in MATLAB by
``matlab/tss_raytrace_vertical.m``, written as ``.mat`` v7.3 into
``data/bridge/``, and this tier reads them back and presents them through the
ordinary :class:`~tss.interfaces.L1Propagation` protocol.  From the pipeline's
point of view it is a tier like any other.  The difference is that it requires a
generation step to have been run first, and when the ray set it needs is not on
disk it must say so with a clear error rather than falling back to the analytic
tier.  A silent fallback would attribute analytic results to the ray-traced
tier in the manifest, which is exactly the failure the registry's duplicate-name
check exists to prevent elsewhere.

The cross-check is not PHaRLAP against itself.  The author has an independent
three-dimensional Haselgrove ray tracer written for the Chartier reproduction,
and it is available to run the same geometries.  Two independent
implementations disagreeing is a finding; one implementation agreeing with
itself is not evidence.  That comparison is the V2 case for L1, and the
analytic tier of :mod:`tss.l1_propagation.analytic` is the third leg: at
vertical incidence in a horizontally uniform ionosphere all three must give the
same group path, and any disagreement there is a bug in one of them rather than
physics.

This tier is also where off-vertical echoes are generated for WP10.  A
plasma depletion wall returns energy from a direction that is not straight down,
and the analytic tier cannot produce such an echo at all because it has no
geometry in which to put one.  The geometry itself is
:mod:`tss.l1_propagation.offvertical`; this tier is what traces rays into it.

Open questions
--------------
1.  What ray set is needed.  A vertical sounding needs one ray per frequency per
    mode, but an off-vertical search needs a fan, and the fan's angular extent
    and spacing decide both the run time and whether a wall echo is found at
    all.  There is no useful guess to make here before the depletion geometry of
    WP10 exists.
2.  How the ray set is keyed.  A cached ray set is only reusable if the
    ionosphere state it was traced through can be identified exactly.  That is a
    provenance question, and :mod:`tss.provenance` is where the answer belongs.
3.  Whether the group path from PHaRLAP and the group path from the analytic
    tier agree at vertical incidence to better than one range cell.  If they do
    not, the analytic tier cannot be used for the trade study, and the whole
    campaign cost changes.

Raises
------
NotImplementedError
    On every call.  WP2, ``docs/matlab-python-bridge.md`` section "Which tiers
    are MATLAB-side, and why", and ``docs/architecture.md`` section "L1
    Propagation".
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from tss.interfaces import EchoSet, IonosphereState, PlatformState
from tss.l1_propagation.base import PropagationBase
from tss.registry import register

__all__ = ["RaytracePropagation"]

_MESSAGE = (
    "WP2: the ray-traced L1 tier is not implemented.  Rays come from PHaRLAP "
    "through the MATLAB file bridge; see docs/matlab-python-bridge.md section "
    "'Which tiers are MATLAB-side, and why' and docs/architecture.md section "
    "'L1 Propagation'.  The analytic tier is not a substitute for it and this "
    "tier will not silently fall back to one."
)


@register("l1", "raytrace")
class RaytracePropagation(PropagationBase):
    """Ray-traced propagation through PHaRLAP.  Stub.

    Registering an unimplemented tier is deliberate.  It lets
    ``tss config <campaign>`` validate a campaign that names ``l1: raytrace``
    before the tier exists, which is the point made in ``CONTRIBUTING.md``:
    a stub raises when it is called and not when it is built.

    Parameters
    ----------
    bridge_dir
        Where to look for the ``.mat`` ray sets the MATLAB side writes.
        Recorded in :meth:`describe` so that a run manifest says which ray sets
        a result came from, once there are any.
    fan_half_angle_deg
        Angular half-extent of the ray fan.  Zero means vertical only.  Carried
        now so that the setting is in the manifest from the first run rather
        than added later and absent from earlier ones.
    """

    def __init__(
        self,
        *,
        bridge_dir: str | Path = "data/bridge",
        fan_half_angle_deg: float = 0.0,
        fan_step_deg: float = 1.0,
        interp: str = "linear",
    ) -> None:
        super().__init__(interp=interp)
        self.bridge_dir = Path(bridge_dir)
        self.fan_half_angle_deg = float(fan_half_angle_deg)
        self.fan_step_deg = float(fan_step_deg)

    def _settings(self) -> dict[str, Any]:
        return {
            "interp": self.interp,
            "bridge_dir": str(self.bridge_dir),
            "fan_half_angle_deg": self.fan_half_angle_deg,
            "fan_step_deg": self.fan_step_deg,
            "engine": "PHaRLAP 4.7.4 via the MATLAB file bridge",
            "status": "stub",
        }

    def echoes(
        self,
        state: IonosphereState,
        platform: PlatformState,
        freq_hz: np.ndarray,
        rng: np.random.Generator,
    ) -> EchoSet:
        """Read a PHaRLAP ray set and present it as an ``EchoSet``.  Not implemented."""
        del state, platform, freq_hz, rng
        raise NotImplementedError(_MESSAGE)
