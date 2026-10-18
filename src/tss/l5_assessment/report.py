"""Layer: L5
WP: WP0, WP7
Claims: C1, C2, C3, C4, C5
Status: implemented
Requirements: -

Render an :class:`~tss.interfaces.Assessment` as a text table and as JSON.

Simple on purpose.  The rendering is the last thing between a result and a
reader, so it should have no logic in it that could change what the result
says.  Everything shown here is read straight off the assessment; nothing is
recomputed, rounded into a different verdict, or omitted for brevity.

Two rules the layout follows.

**The full requirement set is shown every time.**  Requirements that were not
evaluated are listed with their reason, not left out.  A table that showed only
what was scored would read as a clean sweep, which is the failure
``docs/architecture.md`` calls out in the L5 section and accepts the extra
verbosity to avoid.

**The summary line leads with what is missing.**  A run that met nine
requirements and could not evaluate eight is not a run that met nine
requirements.  The counts appear in the order goal, threshold, fail,
not_evaluated, and the not_evaluated count is called out in the header line
rather than being left at the end of a row of numbers.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from tss.interfaces import Assessment
from tss.l5_assessment.base import summarise_verdicts

__all__ = ["render_text", "render_json", "write_report", "to_dict"]

_SYMBOL = {
    "goal": "GOAL",
    "threshold": "THRESH",
    "fail": "FAIL",
    "not_evaluated": "-",
}


def _fmt(value: Any) -> str:
    """Format a number compactly, and anything else as itself."""
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int,)):
        return str(value)
    if isinstance(value, float):
        if value == 0.0:
            return "0"
        magnitude = abs(value)
        if magnitude >= 1.0e5 or magnitude < 1.0e-2:
            return f"{value:.3e}"
        return f"{value:.4g}"
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(_fmt(v) for v in value) + "]"
    return str(value)


def to_dict(assessment: Assessment) -> dict[str, Any]:
    """The assessment as a plain JSON-serialisable dictionary."""
    return {
        "configuration": assessment.configuration,
        "source": assessment.source,
        "summary": summarise_verdicts(assessment.requirements),
        "requirements": {
            rid: dict(record) for rid, record in assessment.requirements.items()
        },
        "errors": dict(assessment.errors),
        "flags": dict(assessment.flags),
    }


def render_json(assessment: Assessment, *, indent: int = 2) -> str:
    """The assessment as a JSON string.

    ``default=str`` is used so that a numpy scalar or a Path that found its way
    into a flag renders as text rather than raising.  A report that fails to
    render because of a stray type is worse than one that renders a repr.
    """
    return json.dumps(to_dict(assessment), indent=indent, default=str)


def render_text(assessment: Assessment, *, width: int = 96) -> str:
    """The assessment as a readable text table."""
    counts = summarise_verdicts(assessment.requirements)
    total = sum(counts.values())

    lines: list[str] = []
    rule = "=" * width
    lines.append(rule)
    lines.append(f"Assessment: {assessment.configuration}")
    revision = assessment.flags.get("requirements_revision", "(unknown revision)")
    lines.append(f"Requirement set: {revision}")
    lines.append(
        f"{total} requirements: {counts['goal']} goal, {counts['threshold']} threshold, "
        f"{counts['fail']} fail, {counts['not_evaluated']} NOT EVALUATED"
    )
    coverage = assessment.flags.get("coverage") or {}
    if coverage.get("n_catalog"):
        held = ", ".join(
            f"{n} {cat}" for cat, n in coverage.get("not_scored_by_category", {}).items() if n
        )
        lines.append(
            f"Catalogue Rev {coverage.get('catalog_revision', '?')}: "
            f"{coverage['n_catalog']} requirements, {coverage.get('n_scored', 0)} scored here, "
            f"{coverage.get('n_not_scored', 0)} NOT SCORED in simulation ({held})"
        )
    lines.append(rule)

    header = (
        f"{'ID':9s} {'verdict':8s} {'achieved':>13s} {'threshold':>13s} "
        f"{'goal':>13s} {'unit':<10s} title"
    )
    lines.append(header)
    lines.append("-" * width)

    for rid in sorted(assessment.requirements):
        record = assessment.requirements[rid]
        verdict = str(record.get("verdict", "not_evaluated"))
        lines.append(
            f"{rid:9s} {_SYMBOL.get(verdict, verdict):8s} "
            f"{_fmt(record.get('achieved')):>13s} "
            f"{_fmt(record.get('threshold')):>13s} "
            f"{_fmt(record.get('goal')):>13s} "
            f"{str(record.get('unit', ''))[:10]:<10s} "
            f"{str(record.get('title', ''))}"
        )
        if verdict == "not_evaluated" and record.get("reason"):
            lines.append(f"{'':9s} {'':8s} reason: {record['reason']}")

    not_scored = assessment.flags.get("not_scored") or {}
    if not_scored:
        lines.append("-" * width)
        lines.append("Not scored in simulation, by catalogue order")
        for rid, record in not_scored.items():
            through = record.get("scored_through") or []
            wps = record.get("work_packages") or []
            extra = (
                f" via {', '.join(through)}" if through
                else (f" ({', '.join(wps)})" if wps else "")
            )
            lines.append(f"{rid:9s} {str(record.get('category', '')):12s}{extra}")
            lines.append(f"{'':9s} reason: {record.get('reason', '')}")

    if assessment.errors:
        lines.append("-" * width)
        lines.append("Errors and achieved numbers")
        for key in sorted(assessment.errors):
            lines.append(f"  {key:36s} {_fmt(assessment.errors[key])}")

    notes = _notable_flags(assessment.flags)
    if notes:
        lines.append("-" * width)
        lines.append("Flags that change how the numbers should be read")
        for key, value in notes:
            lines.append(f"  {key:36s} {_fmt(value)}")

    lines.append(rule)
    return "\n".join(lines)


def _notable_flags(flags: Mapping[str, Any]) -> list[tuple[str, Any]]:
    """The flags a reader needs in order to not misread the table.

    Not every flag: the list is the ones that change the meaning of a number
    rather than describing how it was computed.  Which starting value the
    inversion used, whether the modes were labelled, whether Doppler was
    separated, and whether the trace came from a non-monotonic profile all
    qualify.
    """
    wanted = (
        "l4.starting_value",
        "l4.heights_monotonic",
        "l4.n_nonpositive_laminae",
        "l4.magnetoionic",
        "l4.vertical_incidence_assumed",
        "l3.mode_discrimination",
        "l3.doppler_separated",
        "l3.n_detected",
        "l3.record_truncated_to_one_sweep",
    )
    return [(key, flags[key]) for key in wanted if key in flags]


def write_report(
    assessment: Assessment,
    out_dir: str | Path,
    *,
    stem: str | None = None,
) -> tuple[Path, Path]:
    """Write the text and JSON renderings and return their paths.

    The stem defaults to the configuration name, so a campaign of several
    configurations writes one pair of files each and they do not overwrite one
    another.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    name = stem or assessment.configuration or "assessment"
    text_path = out / f"{name}.txt"
    json_path = out / f"{name}.json"
    text_path.write_text(render_text(assessment) + "\n", encoding="utf-8")
    json_path.write_text(render_json(assessment) + "\n", encoding="utf-8")
    return text_path, json_path
