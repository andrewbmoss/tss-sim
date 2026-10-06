"""Layer: L1
WP: WP2
Claims: C1, C2, C3
Status: partial
Requirements: INS-01, INS-09, INS-10

L1, propagation: turn a truth ionosphere state and a platform state into the
echoes that return to the platform.

The package is thin, as ``CONTRIBUTING.md`` requires of a layer package.  It
imports its tier modules so that registration happens, and re-exports the names
another layer or a study legitimately needs.  There is no physics here.

What is delivered and what is not:

======================  ========  ==========================================
Module                  Status    Note
======================  ========  ==========================================
``absorption``          done      Deviative and non-deviative absorption,
                                  reimplemented from the author's tested link
                                  budget module.
``magnetoionic``        partial   Cutoffs and reflection conditions are real
                                  algebra.  Full Appleton-Hartree with
                                  collisions raises.
``doppler``             done      The forward model C3 has to invert.
``analytic``            done      Vertical-incidence tier, registered as
                                  ``l1/analytic``.
``raytrace``            stub      PHaRLAP through the MATLAB bridge.
``offvertical``         stub      Depletion wall geometry, WP10, goal level.
======================  ========  ==========================================

Both stubs are registered or importable so that a campaign naming them can be
validated before they exist.  They raise :class:`NotImplementedError` when
called, with a message naming the work package.
"""

from __future__ import annotations

from tss.l1_propagation import (  # noqa: F401
    absorption,
    analytic,
    doppler,
    magnetoionic,
    offvertical,
    raytrace,
)
from tss.l1_propagation.absorption import (  # noqa: F401
    absorption_coefficient,
    absorption_db,
    absorption_profile,
    collision_frequency,
    refractive_index,
)
from tss.l1_propagation.analytic import AnalyticPropagation, group_range_m  # noqa: F401
from tss.l1_propagation.base import (  # noqa: F401
    EVANESCENT,
    PENETRATING,
    PropagationBase,
    critical_frequency_hz,
    density_at,
    peak_altitude_m,
    plasma_frequency_at,
)
from tss.l1_propagation.doppler import (  # noqa: F401
    DopplerGeometry,
    apparent_range_bias_m,
    closing_velocity_mps,
    coherent_interval_for_velocity_resolution_s,
    doppler_geometry,
    doppler_shift_hz,
    refractive_index_at_platform,
    residual_bias_second_order_m,
    velocity_from_doppler_mps,
)
from tss.l1_propagation.magnetoionic import (  # noqa: F401
    appleton_hartree,
    cutoff_frequency,
    gyrofrequency_at,
    o_mode_cutoff,
    trace_separation_hz,
    x_mode_cutoff,
    z_mode_cutoff,
)
from tss.l1_propagation.raytrace import RaytracePropagation  # noqa: F401

__all__ = [
    # tiers
    "AnalyticPropagation",
    "RaytracePropagation",
    "PropagationBase",
    # absorption
    "absorption_db",
    "absorption_coefficient",
    "absorption_profile",
    "collision_frequency",
    "refractive_index",
    # magnetoionic
    "cutoff_frequency",
    "o_mode_cutoff",
    "x_mode_cutoff",
    "z_mode_cutoff",
    "gyrofrequency_at",
    "trace_separation_hz",
    "appleton_hartree",
    # doppler
    "DopplerGeometry",
    "doppler_shift_hz",
    "refractive_index_at_platform",
    "velocity_from_doppler_mps",
    "apparent_range_bias_m",
    "residual_bias_second_order_m",
    "coherent_interval_for_velocity_resolution_s",
    "closing_velocity_mps",
    "doppler_geometry",
    # profile helpers and sentinels
    "density_at",
    "plasma_frequency_at",
    "peak_altitude_m",
    "critical_frequency_hz",
    "group_range_m",
    "EVANESCENT",
    "PENETRATING",
]
