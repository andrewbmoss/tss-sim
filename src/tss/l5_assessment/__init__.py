"""Layer: L5
WP: WP0, WP7
Claims: C1, C2, C3, C4, C5
Status: partial
Requirements: INS-02, INS-04, INS-05, INS-07, INS-08, INS-09, INS-10, INS-15, SCI-09, SCI-23

L5, assessment: score the retrieval against the L0 truth it started from, and
score the run against the requirement set, per configuration.

The package is thin, as ``CONTRIBUTING.md`` requires of a layer package.  It
imports its tier modules so that registration happens, and re-exports the names
another layer or a study legitimately needs.

What is delivered and what is not:

======================  ========  ==========================================
Module                  Status    Note
======================  ========  ==========================================
``base``                done      Requirement-set loading and the verdict
                                  summary.
``scoring``             done      Truth-versus-retrieval metrics and the
                                  achieved instrument numbers.
``requirements_check``  done      The metric-to-requirement table and the
                                  scoring call.  A metric that was not
                                  computed gives ``not_evaluated`` with a
                                  reason, never a silent pass.
``standard``            done      Registered as ``l5/standard``.
``report``              done      Text and JSON rendering.
``feasibility``         stub      Trade-space surfaces, WP7.
======================  ========  ==========================================

A note on module names.  ``docs/architecture.md`` lists this layer's tier
module as ``score.py``.  The delivered layout splits it into ``scoring.py``,
which measures, and ``requirements_check.py``, which judges, with
``standard.py`` as the registered tier that composes the two.  The split is
deliberate: changing what is measured and changing what counts as passing are
different changes and should be visible as different diffs.  The architecture
document is stale on this point and, by its own rule, the code is right.
"""

from __future__ import annotations

from tss.l5_assessment import (  # noqa: F401
    base,
    feasibility,
    report,
    requirements_check,
    scoring,
    standard,
)
from tss.l5_assessment.base import (  # noqa: F401
    VERDICTS,
    AssessmentBase,
    summarise_verdicts,
)
from tss.l5_assessment.report import (  # noqa: F401
    render_json,
    render_text,
    to_dict,
    write_report,
)
from tss.l5_assessment.requirements_check import (  # noqa: F401
    MAPPING_TABLE,
    METRIC_FOR_REQUIREMENT,
    NOT_SCORED_HERE,
    achieved_from_metrics,
    check_requirements,
)
from tss.l5_assessment.scoring import (  # noqa: F401
    compute_metrics,
    height_error_at_matched_density_m,
    ionogram_metrics,
    profile_metrics,
)
from tss.l5_assessment.standard import StandardAssessment  # noqa: F401

__all__ = [
    # tiers
    "StandardAssessment",
    "AssessmentBase",
    # scoring
    "compute_metrics",
    "profile_metrics",
    "ionogram_metrics",
    "height_error_at_matched_density_m",
    # requirement mapping
    "check_requirements",
    "achieved_from_metrics",
    "MAPPING_TABLE",
    "METRIC_FOR_REQUIREMENT",
    "NOT_SCORED_HERE",
    # reporting
    "render_text",
    "render_json",
    "write_report",
    "to_dict",
    "summarise_verdicts",
    "VERDICTS",
]
