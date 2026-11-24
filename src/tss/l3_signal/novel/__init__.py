"""Layer: L3
WP: WP5
Claims: C3
Status: partial
Requirements: INS-09, INS-15, SCI-15, SCI-23

The novel L3 algorithms of claim C3: spacecraft motion separates plasma drift
from platform motion.  The estimators here are implemented at the parameter
level and tested by their studies; ``sweep_pairing`` is a stub.

They are separated into their own subpackage because the boundary between the
conventional processing of WP4 and what is claimed as a contribution should be
visible in the directory listing and not only in the prose.

============================  =====  =========================================
Module                        Claim  What it does
============================  =====  =========================================
``route_one_estimator``       C3     Route one: drift from the Doppler zero
                                     crossing and the power peak of one facet
                                     echo, with the extended-reflector refusal.
                                     Parameter level; ``studies/c3_verify``.
``field_aligned_estimator``   C3     Route one for a field-aligned sheet, its
                                     orientation from the field model.
                                     Parameter level;
                                     ``studies/field_aligned_verify``.
``second_look_estimator``     C3     Route two, two nodes: drift and tilt from
                                     two passes.  Parameter level;
                                     ``studies/route_two``.
``network_tilt_estimator``    C3     Route two, oblique links in a
                                     constellation.  Parameter level;
                                     ``studies/route_two_network``.
``sweep_pairing``             C3     Stub.  The ionogram-mode INS-09 bias
                                     check from an up and down sweep pair
                                     (WP5 Phase 5).  Not the contribution.
============================  =====  =========================================

Retired 2026-11-24 (WP5 decision D9): ``motion_compensation`` (the ephemeris
platform-Doppler removal moves into the drift-mode L3 tier of WP5 Phase 3),
``aperture_synthesis`` (exploratory, superseded by the route-one analysis of
TV section 8.8) and ``leakage_cancel`` (claim C4 closed in simulation in WP3
with :class:`tss.l2_instrument.leakage.RFCanceller`).  The composing tier
``motion_separated`` and ``config/campaigns/c3_doppler.yaml`` were retired with
them; the drift-mode tier of WP5 Phase 3 replaces both.

Every function of a stub raises :class:`NotImplementedError` with a message
naming its work package.  A stub that returned a plausible number would let a
campaign produce a result that looks like evidence for C3, and the first
governing principle of Section 3.1 of the proposal is that no figure is scored
until the thing it is scored against exists.
"""

from __future__ import annotations

from tss.l3_signal.novel import (  # noqa: F401
    route_one_estimator,
    sweep_pairing,
)
from tss.l3_signal.novel.sweep_pairing import (  # noqa: F401
    SweepPair,
    associate_traces,
    pair_sweeps,
    range_rate_from_pair_mps,
    residual_bias_cells,
    unbiased_range_m,
)

__all__ = [
    # sweep_pairing, the INS-09 bias check (stub, WP5 Phase 5)
    "SweepPair",
    "pair_sweeps",
    "associate_traces",
    "unbiased_range_m",
    "range_rate_from_pair_mps",
    "residual_bias_cells",
]
