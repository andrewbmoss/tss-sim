"""Layer: L3
WP: WP5
Claims: C3, C4, C2
Status: stub
Requirements: INS-09, INS-15, INS-23, INS-24

The novel L3 algorithms.  Every module in this subpackage is a stub.

These four modules are the research the dissertation proposes, as opposed to
the baseline chain of WP4 which is conventional processing written so that the
novel work has something to be measured against.  They are separated into their
own subpackage for exactly that reason: the boundary between what is standard
and what is claimed as a contribution should be visible in the directory
listing and not only in the prose.

======================  ======  =====  ==============================================
Module                  WP      Claim  What it would do
======================  ======  =====  ==============================================
``sweep_pairing``       WP5     C3     Unbiased range from the half-sum of an up and
                                       down sweep pair, radial velocity from the
                                       half-difference.  Principal contribution.
``motion_compensation`` WP5     C3     Coherent integration over the known ephemeris
                                       across the interval INS-15 demands.
``aperture_synthesis``  WP5     C2     Along-track aperture from the orbital motion.
                                       Exploratory extension, not a claim.
``leakage_cancel``      WP5     C4     Adaptive cancellation of transmit-to-receive
                                       leakage.  The claim the concept now turns on.
======================  ======  =====  ==============================================

Every function raises :class:`NotImplementedError` with a message naming its
work package.  That is the point.  A stub that returned a plausible number
would let a campaign produce a result that looks like evidence for C3, and the
first governing principle of Section 3.1 of the proposal is that no figure is
scored until the thing it is scored against exists.  A fabricated intermediate
is worse than a missing one because it cannot be found by looking at what is
missing.

The tier that composes these, :mod:`tss.l3_signal.motion_separated`, is
registered so that ``config/campaigns/c3_doppler.yaml`` validates before the
algorithms exist.  It raises when it is run.
"""

from __future__ import annotations

from tss.l3_signal.novel import (  # noqa: F401
    aperture_synthesis,
    leakage_cancel,
    motion_compensation,
    sweep_pairing,
)
from tss.l3_signal.novel.aperture_synthesis import (  # noqa: F401
    along_track_resolution_m,
    synthesize_aperture,
)
from tss.l3_signal.novel.leakage_cancel import (  # noqa: F401
    cancel_leakage,
    estimate_leakage,
    residual_floor_db,
)
from tss.l3_signal.novel.motion_compensation import (  # noqa: F401
    compensate_motion,
    integrate_coherently,
    phase_history_rad,
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
    # sweep_pairing, C3
    "SweepPair",
    "pair_sweeps",
    "associate_traces",
    "unbiased_range_m",
    "range_rate_from_pair_mps",
    "residual_bias_cells",
    # motion_compensation, C3
    "phase_history_rad",
    "compensate_motion",
    "integrate_coherently",
    # aperture_synthesis, C2, exploratory
    "synthesize_aperture",
    "along_track_resolution_m",
    # leakage_cancel, C4
    "estimate_leakage",
    "cancel_leakage",
    "residual_floor_db",
]
