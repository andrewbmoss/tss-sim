"""Layer: L5
WP: WP7
Claims: C1, C2
Status: stub
Requirements: INS-05, INS-10, INS-14, SCI-09

Feasibility surfaces over the configuration space.  Not implemented.  WP7.

What it would produce
----------------------
``config/campaigns/trade_space.yaml`` sweeps antenna length, transmit power and
sweep rate across four configurations and runs a full pipeline at every point.
Each point yields an :class:`~tss.interfaces.Assessment`.  This module is what
turns that pile of assessments into the surfaces Section 1.3 of the proposal
reports: for each configuration, the region of the design space in which every
threshold requirement closes, and inside that, the region in which the goals
close.

The output is per configuration and not collapsed to a single verdict.  A
single number over the whole space would hide the thing the trade study exists
to show, which is which requirement binds where.  The useful picture is the
boundary of the feasible region annotated with the requirement that forms each
segment of it: sensitivity at the small-antenna end, range resolution at the
low-sweep-rate end, along-track sampling at the long-period end.

Why it is a stub rather than an aggregation of assessments
------------------------------------------------------------
The aggregation is easy and it is not the difficulty.  Three questions have to
be settled first and none of them is settled:

1.  **What counts as feasible when requirements are ``not_evaluated``.**  Most
    points in the current trade space will have eight or more requirements
    unevaluated, because the tiers that would measure them do not exist.  A
    feasible region computed from the requirements that happen to be scorable
    today would grow as the codebase grows, which makes it useless as a record.
    The region has to be defined against a stated subset and the subset has to
    be part of the result.
2.  **How to combine configurations.**  The four configurations are not four
    samples of one space; they are four different instruments.  Overlaying
    their feasible regions on one axis pair would suggest a continuum that is
    not there.
3.  **What the axes are.**  Antenna length, transmit power and sweep rate are
    three axes and a surface has two.  Which pair is shown and what the third
    is held at is a presentation decision that changes what the reader
    concludes, and it should be made once, deliberately, and recorded.

Until those are answered, the aggregation would produce a picture that looks
authoritative and is not, which is the specific failure the stub policy of this
package exists to prevent.

Raises
------
NotImplementedError
    On every call.  WP7.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from tss.interfaces import Assessment

__all__ = [
    "feasible_mask",
    "feasibility_surface",
    "binding_requirement",
]

_MESSAGE = (
    "WP7: feasibility surfaces are not implemented.  Aggregating assessments is "
    "easy; what is not settled is what counts as feasible when most requirements "
    "are still not_evaluated, how four different configurations share an axis, and "
    "which pair of the three trade axes a surface shows.  See the module docstring "
    "of tss.l5_assessment.feasibility and docs/work-packages.md row WP7."
)


def feasible_mask(
    assessments: Sequence[Assessment],
    *,
    level: str = "threshold",
    required_ids: Sequence[str] | None = None,
) -> np.ndarray:
    """Boolean mask of which runs closed at the stated level.  Not implemented.

    ``required_ids`` is in the signature because a feasible region has to be
    defined against a stated subset of requirements; see question one in the
    module docstring.  A default of "all of them" would make the region grow as
    the codebase grows.

    Raises
    ------
    NotImplementedError
        Always.  WP7.
    """
    del assessments, level, required_ids
    raise NotImplementedError(_MESSAGE)


def binding_requirement(
    assessments: Sequence[Assessment],
) -> np.ndarray:
    """Which requirement binds at each point of the space.  Not implemented.

    This is the annotation that makes a feasibility plot informative rather
    than decorative: the boundary of a feasible region is made of segments and
    each segment is a different requirement.

    Raises
    ------
    NotImplementedError
        Always.  WP7.
    """
    del assessments
    raise NotImplementedError(_MESSAGE)


def feasibility_surface(
    assessments: Sequence[Assessment],
    axes: Mapping[str, Sequence[float]],
    *,
    configuration: str | None = None,
) -> dict[str, Any]:
    """Feasible region over a pair of trade axes, per configuration.  Not implemented.

    Raises
    ------
    NotImplementedError
        Always.  WP7.
    """
    del assessments, axes, configuration
    raise NotImplementedError(_MESSAGE)
