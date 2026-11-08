"""Cross-language transcription check for the link budget.

Layer:        study
WP:           WP8
Claims:       C1
Status:       stub
Requirements: -

Verification level V2.  The MATLAB port of the link budget is to be checked
against the Python reference over the same inputs and must agree to machine
precision.  A port that drifts from its reference is the failure mode this
catches, and it is a real one, because the two are edited at different times by
different hands.  The procedure and the tolerance standard are in
``docs/matlab-python-bridge.md``.

Status: stub, owned by WP8.  The MATLAB-side verifier that procedure names,
``matlab/tss_verify_transcription.m``, and a Python driver that calls it are
not yet in the repository, so the check cannot run whether or not MATLAB is
installed.  Until 2026-11-08 this study was registered ``implemented``, probed
for MATLAB and reported a missing input that pointed at MATLAB test files that
did not exist either (code review F6).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from tss.studies import StudyResult, register_study


@register_study(
    "transcription",
    description="Check the MATLAB link budget against the Python reference",
    work_packages=["WP8"],
    claims=["C1"],
    status="stub",
    blocked_on="WP8",
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    raise NotImplementedError(
        "WP8: the MATLAB-side transcription verifier (matlab/tss_verify_transcription.m, "
        "described in docs/matlab-python-bridge.md) and a Python driver for it are not "
        "yet in the repository, so the link-budget transcription check cannot run."
    )
