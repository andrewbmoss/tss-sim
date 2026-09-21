"""Feasibility surfaces across the configuration space.

Layer:        study
WP:           WP7
Claims:       C2
Status:       stub
Requirements: INS-01, INS-02, INS-05, INS-07, INS-14, SCI-09

Section 1.3 of the proposal says the work is a design and not a point design, and
that results are reported per configuration.  This study is what makes that true
in practice.  It sweeps the block named ``sweep_over`` in
``config/campaigns/trade_space.yaml`` (antenna length, transmit power, sweep
rate) across the named configurations, runs each point through the full L0 to L5
pipeline, and produces a surface showing what closes where and with what margin.

It is a stub because it is downstream of everything.  A feasibility surface built
on an L3 tier that cannot yet separate Doppler from range would be a surface of
the wrong function.  WP4 has to land, then WP5, before this produces a number
worth plotting.

The design is not complicated.  The value is in reporting per configuration
rather than collapsing to a single verdict, and in naming the configuration
alongside every result so that a result without its configuration never escapes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from tss.studies import StudyResult, register_study


@register_study(
    "trade_space",
    description="Feasibility surfaces per configuration for WP7",
    work_packages=["WP7"],
    claims=["C2"],
    status="stub",
    requirements=["INS-01", "INS-02", "INS-05", "INS-07", "INS-14", "SCI-09"],
    blocked_on="WP4",
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    raise NotImplementedError(
        "WP7 is downstream of WP4 and WP5.  The campaign shape is already fixed in "
        "config/campaigns/trade_space.yaml, so this study becomes a matter of "
        "running it once the baseline signal chain and the novel algorithms exist.  "
        "See docs/work-packages.md for the sequencing."
    )
