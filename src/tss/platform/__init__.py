"""Layer: platform
WP: WP0, WP1
Claims: C3
Status: implemented
Requirements: MIS-03

The platform: where the spacecraft is and how fast it is going.

This is not one of the L0 to L5 layers.  It sits beside them, because the
platform state is an input to L1 and to the novel L3 algorithms rather than
something any single layer owns.  :class:`~tss.interfaces.PlatformState` is
defined in :mod:`tss.interfaces` for the same reason.

The package is thin.  It imports its one module and re-exports the names another
layer or a study needs.  There is no physics here.

======================  ========  ==========================================
Module                  Status    Note
======================  ========  ==========================================
``orbits``              done      Circular orbits with secular J2 nodal
                                  precession, the O1 and O2 orbit classes,
                                  the sun-synchronous inclination MIS-03
                                  excludes, and the platform state with an
                                  analytic Earth-fixed velocity.
======================  ========  ==========================================

Nothing here registers a tier, because the platform is not a layer and the
registry has no slot for it.  A campaign names an orbit by its parameters, not
by a tier name.
"""

from __future__ import annotations

from tss.platform import orbits  # noqa: F401
from tss.platform.orbits import (  # noqa: F401
    EARTH_ROTATION_RAD_PER_S,
    J2,
    MU,
    RE,
    SIDEREAL_DAY,
    SOLAR_DAY,
    SUN_DRIFT_DEG_PER_DAY,
    CircularOrbit,
    orbit_classes,
    platform_state_at,
    sun_synchronous_inclination,
)

__all__ = [
    "CircularOrbit",
    "orbit_classes",
    "sun_synchronous_inclination",
    "platform_state_at",
    "MU",
    "RE",
    "J2",
    "SOLAR_DAY",
    "SIDEREAL_DAY",
    "SUN_DRIFT_DEG_PER_DAY",
    "EARTH_ROTATION_RAD_PER_S",
]
