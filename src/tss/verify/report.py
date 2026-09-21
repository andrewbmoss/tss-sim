"""Layer: support
WP: WP8
Claims: C1, C2, C3, C4, C5
Status: implemented
Requirements: -

Write the verification report as JSON and as Markdown.

The layout follows ``docs/verification.md``: one section per level, one row per
check, carrying the verdict, the quantity compared, the observed discrepancy
and the tolerance.  Two properties of the report are deliberate.

**A weak report says so at the top.**  ``docs/verification.md`` states that a
report with many skips is a weak report and must say so.  The header line
carries the counts in the order that makes that visible, and when more checks
did not run than ran, the header says that in words.  A reader who sees only
"no failures" has been misled by omission.

**Nothing is omitted.**  Every registered check appears, including the ones
that cannot run and the one at V4 that records a permanent limitation rather
than a gap.  The report is the full registry every time.

The report carries the requirement set revision and the git state, so that a
verification result can be traced to the code and the requirement numbers that
produced it.  A report from a dirty working tree records that it was.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from tss.verify.levels import LEVELS, VERDICTS, CheckResult

__all__ = ["summarise", "to_dict", "render_markdown", "write_report", "default_out_dir"]

_LEVEL_TITLE = {
    "V1": "V1  Internal consistency",
    "V2": "V2  Cross-fidelity",
    "V3": "V3  Literature replication",
    "V4": "V4  Archived real measurements",
}


def default_out_dir() -> Path:
    """``results/verification`` under the repository root."""
    from tss.config import repo_root

    return repo_root() / "results" / "verification"


def summarise(results: Mapping[str, Sequence[CheckResult]]) -> dict[str, int]:
    """Verdict counts across every level, with zero counts kept."""
    counts = {verdict: 0 for verdict in VERDICTS}
    for records in results.values():
        for record in records:
            counts[record.verdict] = counts.get(record.verdict, 0) + 1
    return counts


def to_dict(results: Mapping[str, Sequence[CheckResult]]) -> dict[str, Any]:
    """The whole report as a JSON-serialisable dictionary."""
    from tss.provenance import git_state

    try:
        from tss.requirements import load_requirements

        revision = load_requirements().revision
    except Exception as exc:  # pragma: no cover - only when the file is missing
        revision = f"unavailable: {exc}"

    try:
        git = git_state()
    except Exception as exc:  # pragma: no cover - not a git checkout
        git = {"error": str(exc)}

    counts = summarise(results)
    total = sum(counts.values())
    return {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "requirements_revision": revision,
        "git": git,
        "summary": counts,
        "total_checks": total,
        "weak": counts["not_run"] + counts["skipped"] > total - counts["not_run"] - counts["skipped"],
        "levels": {
            level: [record.as_dict() for record in results.get(level, ())]
            for level in LEVELS
        },
    }


def render_markdown(results: Mapping[str, Sequence[CheckResult]]) -> str:
    """The report as Markdown."""
    payload = to_dict(results)
    counts = payload["summary"]
    total = int(payload["total_checks"])
    ran = counts["pass"] + counts["fail"] + counts["partial"]

    lines: list[str] = []
    lines.append("# Verification report")
    lines.append("")
    lines.append(f"Generated {payload['generated']}.")
    lines.append(f"Requirement set: {payload['requirements_revision']}.")
    git = payload.get("git") or {}
    if git.get("commit"):
        dirty = " (dirty working tree)" if git.get("dirty") else ""
        lines.append(f"Git: {git['commit']}{dirty}.")
    lines.append("")
    lines.append(
        f"{total} checks: {counts['pass']} pass, {counts['fail']} fail, "
        f"{counts['partial']} partial, {counts['skipped']} skipped, "
        f"{counts['not_run']} not run."
    )
    if payload["weak"]:
        lines.append("")
        lines.append(
            f"**This is a weak report.**  Only {ran} of {total} checks actually ran.  "
            "The rest could not, either because the physics they would compare is "
            "not written or because an external dependency is absent.  Read the "
            "reasons before reading the verdicts."
        )
    lines.append("")

    for level in LEVELS:
        records = list(results.get(level, ()))
        lines.append(f"## {_LEVEL_TITLE.get(level, level)}")
        lines.append("")
        if not records:
            lines.append("No checks registered at this level.")
            lines.append("")
            continue
        lines.append("| Check | Verdict | Quantity | Observed | Tolerance |")
        lines.append("| --- | --- | --- | --- | --- |")
        for record in records:
            lines.append(
                f"| `{record.name}` | {record.verdict} | {record.quantity or '-'} | "
                f"{_fmt(record.observed)} | {_fmt(record.tolerance)} |"
            )
        lines.append("")
        for record in records:
            note = record.reason or record.detail
            if note:
                lines.append(f"* `{record.name}`: {note}")
        lines.append("")

    return "\n".join(lines)


def _fmt(value: float | None) -> str:
    if value is None:
        return "-"
    if value == 0.0:
        return "0"
    if abs(value) >= 1.0e5 or abs(value) < 1.0e-2:
        return f"{value:.3e}"
    return f"{value:.4g}"


def write_report(
    results: Mapping[str, Sequence[CheckResult]],
    out_dir: str | Path | None = None,
    *,
    stem: str = "verification",
) -> tuple[Path, Path]:
    """Write the JSON and Markdown renderings and return their paths."""
    out = Path(out_dir) if out_dir is not None else default_out_dir()
    out.mkdir(parents=True, exist_ok=True)
    json_path = out / f"{stem}.json"
    md_path = out / f"{stem}.md"
    json_path.write_text(
        json.dumps(to_dict(results), indent=2, default=str) + "\n", encoding="utf-8"
    )
    md_path.write_text(render_markdown(results) + "\n", encoding="utf-8")
    return json_path, md_path
