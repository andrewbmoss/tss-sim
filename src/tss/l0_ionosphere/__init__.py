"""Layer: L0
WP: WP1
Claims: C1, C2
Status: partial
Requirements: INS-01

L0, ionosphere state: produce the truth ionosphere at a location and epoch on a
given altitude grid.  Nothing in this layer knows that a radar exists.

The package is thin, as ``CONTRIBUTING.md`` requires of a layer package.  It
imports its tier modules so that registration happens, and re-exports the names
another layer or a study legitimately needs.  There is no physics here.

What is delivered and what is not:

======================  ========  ==========================================
Module                  Status    Note
======================  ========  ==========================================
``chapman``             done      Closed-form alpha-Chapman parameterised by
                                  the plasma scale height.  The default tier
                                  and the only one that runs with no optional
                                  dependency installed.
``backends``            partial   The migrated model backends.  ``analytic``
                                  runs anywhere; ``iri_wrapper`` and
                                  ``matlab_iri`` run where their model does;
                                  ``iri2020`` is registered with its density
                                  call deliberately unbound.
``geomag``              done      IGRF through ppigrf when it is importable,
                                  a tilted centred dipole when it is not.
``base``                done      Shared helpers and the legacy profile to
                                  ``IonosphereState`` conversion.
``perturbations``       stub      Plasma bubble depletions, WP10, goal level.
======================  ========  ==========================================

The stub is importable so that a campaign naming it can be validated before it
exists.  It raises :class:`NotImplementedError` when called, with a message
naming the work package.
"""

from __future__ import annotations

from tss.l0_ionosphere import (  # noqa: F401
    backends,
    base,
    chapman,
    geomag,
    perturbations,
)
from tss.l0_ionosphere.backends import (  # noqa: F401
    EXTRA_COLUMNS,
    IRI_TOPSIDE_OPTIONS,
    PHARLAP_DEFAULT,
    AnalyticBackend,
    Backend,
    IRIWrapperBackend,
    MatlabIRIBackend,
    PyIRIBackend,
    make_backend,
)
from tss.l0_ionosphere.base import (  # noqa: F401
    IonosphereBase,
    collision_frequency,
    default_te_k,
    parse_epoch,
    r12_to_f107,
    state_from_profile,
    validate_alt_grid,
)
from tss.l0_ionosphere.chapman import (  # noqa: F401
    ChapmanIonosphere,
    chapman_density,
    chapman_scale_height_m,
)
from tss.l0_ionosphere.geomag import (  # noqa: F401
    dipole_field_ned_t,
    field_ned_t,
    geomagnetic_latitude_deg,
)
from tss.l0_ionosphere.perturbations import (  # noqa: F401
    BubbleParameters,
    bubble_depletion_factor,
    perturb_state,
)

__all__ = [
    # tiers
    "ChapmanIonosphere",
    "AnalyticBackend",
    "PyIRIBackend",
    "IRIWrapperBackend",
    "MatlabIRIBackend",
    "Backend",
    "IonosphereBase",
    # construction
    "make_backend",
    # profile helpers
    "chapman_density",
    "chapman_scale_height_m",
    "state_from_profile",
    "validate_alt_grid",
    "default_te_k",
    "collision_frequency",
    "parse_epoch",
    "r12_to_f107",
    # geomagnetic field
    "field_ned_t",
    "dipole_field_ned_t",
    "geomagnetic_latitude_deg",
    # WP10 stub
    "BubbleParameters",
    "bubble_depletion_factor",
    "perturb_state",
    # constants
    "IRI_TOPSIDE_OPTIONS",
    "EXTRA_COLUMNS",
    "PHARLAP_DEFAULT",
]
