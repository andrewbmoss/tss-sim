"""Layer: L4
WP: WP6, WP10
Claims: C2
Status: partial
Requirements: INS-05

L4, retrieval: invert virtual height to true height, giving electron density
against true altitude with an uncertainty.

The package is thin, as ``CONTRIBUTING.md`` requires of a layer package.  It
imports its tier modules so that registration happens, and re-exports the names
another layer or a study legitimately needs.  There is no inversion physics
here.

What is delivered and what is not:

================  ========  ==================================================
Module            Status    Note
================  ========  ==================================================
``base``          done      Trace reading, the starting-value policy and the
                            reflection-density relation.
``lamination``    done      ``l4/lamination`` is the classical no-field case.
                            ``l4/lamination_magnetoionic`` inverts the
                            ordinary trace with the magnetoionic group index
                            and is the tier a profile is scored from.
``error_budget``  stub      Propagated uncertainty, WP6.  Until it exists,
                            ``ne_sigma`` and ``alt_sigma_m`` stay ``None``
                            and no retrieval accuracy can be quoted.
``depletion``     stub      Non-monotonic crossings, WP10, goal level.
================  ========  ==================================================

The dependency worth stating: L4 is gated by L3.  Until WP4 delivers an
ionogram from a simulated receiver there is nothing to invert except synthetic
ionograms, and inverting a synthetic ionogram tests the inversion arithmetic
and nothing about the instrument.  The arithmetic test is still worth having
and it is what ``tests/unit/test_l4_retrieval.py`` does.
"""

from __future__ import annotations

from tss.l4_retrieval import (  # noqa: F401
    base,
    depletion,
    error_budget,
    lamination,
)
from tss.l4_retrieval.base import (  # noqa: F401
    STARTING_VALUE_SOURCES,
    RetrievalBase,
    density_from_trace,
    starting_plasma_frequency_hz,
    trace_points,
)
from tss.l4_retrieval.error_budget import with_uncertainty  # noqa: F401
from tss.l4_retrieval.lamination import (  # noqa: F401
    LaminationRetrieval,
    MagnetoionicLaminationRetrieval,
    constant_field,
    geomagnetic_field_model,
    lamina_thicknesses_m,
    lamina_thicknesses_magnetoionic_m,
    lamination_matrix,
    magnetoionic_lamination_row,
    virtual_range_from_laminae_m,
    virtual_range_from_laminae_magnetoionic_m,
)

__all__ = [
    # tiers
    "LaminationRetrieval",
    "MagnetoionicLaminationRetrieval",
    "RetrievalBase",
    # lamination algebra
    "lamination_matrix",
    "lamina_thicknesses_m",
    "virtual_range_from_laminae_m",
    "magnetoionic_lamination_row",
    "lamina_thicknesses_magnetoionic_m",
    "virtual_range_from_laminae_magnetoionic_m",
    "constant_field",
    "geomagnetic_field_model",
    # helpers
    "trace_points",
    "starting_plasma_frequency_hz",
    "density_from_trace",
    "STARTING_VALUE_SOURCES",
    # error budget
    "with_uncertainty",
]
