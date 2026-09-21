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
``lamination``    partial   Registered as ``l4/lamination``.  The classical
                            no-field case is complete.  The magnetoionic
                            correction is not written and is the largest
                            known bias in the tier.
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
    lamina_thicknesses_m,
    lamination_matrix,
    virtual_range_from_laminae_m,
)

__all__ = [
    # tiers
    "LaminationRetrieval",
    "RetrievalBase",
    # lamination algebra
    "lamination_matrix",
    "lamina_thicknesses_m",
    "virtual_range_from_laminae_m",
    # helpers
    "trace_points",
    "starting_plasma_frequency_hz",
    "density_from_trace",
    "STARTING_VALUE_SOURCES",
    # error budget
    "with_uncertainty",
]
