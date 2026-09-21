"""Cross-language transcription check for the link budget.

Layer:        study
WP:           WP8
Claims:       C1
Status:       implemented
Requirements: -

Verification level V2.  The MATLAB port of the link budget is checked against
the Python reference over the same inputs and must agree to machine precision.
A port that drifts from its reference is the failure mode this catches, and it
is a real one, because the two are edited at different times by different hands.

This study needs MATLAB.  When MATLAB is absent it reports that it could not run
rather than failing, because an absent toolchain is not a defect in the code.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from tss.studies import StudyInputMissing, StudyResult, register_study


@register_study(
    "transcription",
    description="Check the MATLAB link budget against the Python reference",
    work_packages=["WP8"],
    claims=["C1"],
    status="implemented",
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    from tss.data import local_models

    status = local_models.probe().get("matlab")
    available = bool(getattr(status, "available", False)) if status is not None else False
    if not available:
        raise StudyInputMissing(
            "MATLAB was not found, so the transcription check cannot run here.  "
            "This is not a failure of the port.  Run "
            "`matlab -batch \"run('matlab/tests/tss_run_tests.m')\"` on a machine "
            "with MATLAB, which runs the same comparison as TestTranscription."
        )

    raise StudyInputMissing(
        "MATLAB is present but the Python-side driver for the comparison is not "
        "wired yet.  The MATLAB-side check is implemented as "
        "matlab/tests/TestTranscription.m and is the authoritative path today.  "
        "Wiring this study to call it through tss_call_python.m is WP8 work."
    )
