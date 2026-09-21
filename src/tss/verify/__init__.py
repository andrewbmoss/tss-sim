"""Layer: support
WP: WP8
Claims: C1, C2, C3, C4, C5
Status: partial
Requirements: -

The verification suite, V1 to V4.

``docs/verification.md`` is the authority on what each level answers and, more
importantly, on what each level cannot reach.  This package is the executable
half of it: :mod:`tss.verify.levels` holds the check registry, the runner and
the checks themselves, and :mod:`tss.verify.report` writes the result as JSON
and Markdown under ``results/verification/``.

``tss verify --level all --write`` runs everything available on the machine and
writes the report.  Checks that could not run appear with their reason.  They
are never silently omitted, and a check whose subject is an unwritten algorithm
reports ``not_run`` with the work package named rather than passing vacuously.

The present state, which the report states for itself every time it runs: the
V1 closed-form cases are real and they cover the beat-to-range mapping, the
per-cell noise bandwidth, the transform, the Doppler bias antisymmetry, the
range-correlation limit, the lamination inversion and the requirement scoring.
Everything at V2, V3 and V4 is registered and not run, because the
implementations or the reference data they would compare against do not exist
yet.  That is a weak verification position and the report says so in the
header.
"""

from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path

from tss.verify import levels, report  # noqa: F401
from tss.verify.levels import (  # noqa: F401
    LEVELS,
    VERDICTS,
    CheckResult,
    NotRun,
    Skipped,
    check,
    registered,
    run_all,
    run_level,
)
from tss.verify.report import (  # noqa: F401
    render_markdown,
    summarise,
    to_dict,
    write_report,
)

__all__ = [
    "run_verification",
    "run_level",
    "run_all",
    "registered",
    "check",
    "NotRun",
    "Skipped",
    "CheckResult",
    "LEVELS",
    "VERDICTS",
    "summarise",
    "to_dict",
    "render_markdown",
    "write_report",
]


def run_verification(
    level: str = "all",
    write: bool = False,
    *,
    out: str | Path | None = None,
    stream=None,
) -> int:
    """Run the verification suite and print a summary.  Called by ``tss verify``.

    Parameters
    ----------
    level
        One of ``"V1"``, ``"V2"``, ``"V3"``, ``"V4"`` or ``"all"``.
    write
        Write the report under ``results/verification/``.
    out
        Override the output directory.
    stream
        Where to print.  Defaults to standard output.

    Returns
    -------
    int
        0 when no check failed, 1 when at least one did.  A check that did not
        run is not a failure, and a run in which nothing ran and nothing failed
        still returns 0.  The exit status alone is therefore not evidence of
        verification, which is why the summary line reports how many checks
        actually ran and the written report says so in its header.
    """
    out_stream = stream if stream is not None else sys.stdout
    levels_to_run: Sequence[str] = LEVELS if level == "all" else (level,)
    for name in levels_to_run:
        if name not in LEVELS:
            raise ValueError(f"unknown verification level {name!r}; expected {LEVELS}")

    results = {name: run_level(name) for name in levels_to_run}
    counts = summarise(results)
    total = sum(counts.values())
    ran = counts["pass"] + counts["fail"] + counts["partial"]

    for name in levels_to_run:
        print(f"{name}:", file=out_stream)
        for record in results[name]:
            note = record.reason or record.detail
            observed = "" if record.observed is None else f"  ({record.observed:.4g})"
            print(f"  {record.verdict:9s} {record.name}{observed}", file=out_stream)
            if note and record.verdict != "pass":
                print(f"            {note}", file=out_stream)

    print(
        f"\n{total} checks: {counts['pass']} pass, {counts['fail']} fail, "
        f"{counts['partial']} partial, {counts['skipped']} skipped, "
        f"{counts['not_run']} not run.  {ran} of {total} actually ran.",
        file=out_stream,
    )

    if write:
        json_path, md_path = write_report(results, out)
        print(f"wrote {json_path}\nwrote {md_path}", file=out_stream)

    return 1 if counts["fail"] else 0
