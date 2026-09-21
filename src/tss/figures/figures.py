"""Layer: support
WP: WP8
Claims: C1, C2, C3, C4, C5
Status: partial
Requirements: -

The named figure registry.

Every figure in the dissertation is produced by ``tss figure <name>`` from a
stored result, which is the reproducibility promise of Section 4.2 made
operational.  A figure function loads a result file, draws, and writes; it
never computes a result.  If a figure needs a number the stored result does not
carry, the fix is to put the number in the result, not to recompute it here,
because a number recomputed at plotting time is not the number that was scored.

Declared but not written
------------------------
Most of the dissertation's figures are registered here with a description, an
input path and ``status="declared"``.  They are declared rather than omitted so
that the figure inventory is a list in the code rather than a list in someone's
head, and so that ``tss figure <name>`` on an unwritten figure says which study
has to run first instead of "unknown figure".

Asking for a declared figure returns exit status 2 and prints what is missing.
Asking for an implemented figure whose input is not on disk returns 3 and names
the file and the study that writes it.  Neither draws anything, which is the
point: a figure drawn from data that is not there would be a figure of nothing.

The two implemented figures both read files this package writes itself, so
their data shape is fixed by :mod:`tss.verify.report` and
:mod:`tss.l5_assessment.report` rather than guessed at.  The rest wait on the
studies that produce their inputs.
"""

from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

__all__ = [
    "Figure",
    "FIGURES",
    "register_figure",
    "available_figures",
    "default_figure_dir",
    "results_dir",
]


@dataclass(frozen=True)
class Figure:
    """One named figure: what it shows, what it reads, and whether it exists."""

    name: str
    description: str
    inputs: tuple[str, ...]
    status: str                       # "implemented" or "declared"
    produced_by: str = ""             # which study or command writes the input
    draw: Callable[..., list[Path]] | None = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "inputs": list(self.inputs),
            "status": self.status,
            "produced_by": self.produced_by,
        }


FIGURES: dict[str, Figure] = {}


def results_dir() -> Path:
    """``results/`` under the repository root."""
    from tss.config import repo_root

    return repo_root() / "results"


def default_figure_dir() -> Path:
    """``results/figures/``, the path ``config/default.yaml`` names."""
    return results_dir() / "figures"


def register_figure(
    name: str,
    *,
    description: str,
    inputs: tuple[str, ...] = (),
    status: str = "declared",
    produced_by: str = "",
) -> Callable[[Callable[..., list[Path]]], Callable[..., list[Path]]]:
    """Register a figure.  Used as a decorator for the implemented ones."""

    def _decorate(fn: Callable[..., list[Path]]) -> Callable[..., list[Path]]:
        if name in FIGURES:
            raise ValueError(f"figure {name!r} is already registered")
        FIGURES[name] = Figure(
            name=name,
            description=description,
            inputs=inputs,
            status=status,
            produced_by=produced_by,
            draw=fn,
        )
        return fn

    return _decorate


def _declare(name: str, description: str, inputs: tuple[str, ...], produced_by: str) -> None:
    """Register a figure that is described but not yet written."""
    if name in FIGURES:
        raise ValueError(f"figure {name!r} is already registered")
    FIGURES[name] = Figure(
        name=name,
        description=description,
        inputs=inputs,
        status="declared",
        produced_by=produced_by,
        draw=None,
    )


def available_figures() -> dict[str, dict[str, Any]]:
    """Every registered figure, for ``tss figure list`` and for tests."""
    return {name: fig.as_dict() for name, fig in sorted(FIGURES.items())}


# ---------------------------------------------------------------------------
# Implemented figures
# ---------------------------------------------------------------------------


@register_figure(
    "verification_summary",
    description=(
        "Verdict counts per verification level, showing how much of the suite "
        "actually ran.  The not-run bar is the point of the figure."
    ),
    inputs=("results/verification/verification.json",),
    status="implemented",
    produced_by="tss verify --level all --write",
)
def _fig_verification_summary(out: Path) -> list[Path]:
    """Stacked bars of verification verdicts by level.

    The figure is drawn so that ``not_run`` and ``skipped`` are visually
    dominant when they dominate the counts, because a verification report whose
    checks mostly did not run should not be able to look like a clean sweep.
    """
    import numpy as np

    from tss.figures import style

    payload = _load_json(results_dir() / "verification" / "verification.json")
    levels = list(payload["levels"])
    verdict_order = ("pass", "partial", "fail", "skipped", "not_run")

    counts = {
        verdict: np.array(
            [
                sum(1 for r in payload["levels"][lev] if r["verdict"] == verdict)
                for lev in levels
            ],
            dtype=float,
        )
        for verdict in verdict_order
    }

    style.use()
    fig, ax = style.figure(figsize=(7.0, 3.6))
    bottom = np.zeros(len(levels))
    for verdict in verdict_order:
        values = counts[verdict]
        if not values.any():
            continue
        ax.bar(
            levels,
            values,
            bottom=bottom,
            label=verdict.replace("_", " "),
            color=style.VERDICT_COLOURS.get(verdict, style.INK_MUTED),
            width=0.62,
            zorder=3,
        )
        bottom += values

    ax.set_ylabel("checks")
    ax.set_title("Verification coverage by level")
    ax.legend(ncol=len(verdict_order), loc="upper right")
    style.grid(ax, axis="y")

    ran = payload["summary"]["pass"] + payload["summary"]["fail"] + payload["summary"]["partial"]
    style.caption(
        fig,
        f"{ran} of {payload['total_checks']} registered checks ran.  "
        f"Requirement set {payload['requirements_revision']}.  "
        "A check reported not run has no implementation to test; a skipped one "
        "has no MATLAB, PHaRLAP or cached data on this machine.",
    )
    return style.savefig(fig, Path(out) / "verification_summary")


@register_figure(
    "requirements_scorecard",
    description=(
        "Per-requirement verdict for one configuration, with the achieved value "
        "against threshold and goal.  Unevaluated requirements are shown, not "
        "dropped."
    ),
    inputs=("results/assessments/<configuration>.json",),
    status="implemented",
    produced_by="tss run <campaign>, via tss.l5_assessment.report.write_report",
)
def _fig_requirements_scorecard(out: Path, *, configuration: str | None = None) -> list[Path]:
    """Horizontal verdict bar per requirement for one assessment.

    Every requirement in the set appears, in ID order, including the ones that
    were not evaluated.  ``docs/architecture.md`` accepts the extra verbosity
    for the reason this figure makes visual: a scorecard that showed only what
    was scored would read as a clean sweep.
    """
    from tss.figures import style

    directory = results_dir() / "assessments"
    if configuration:
        path = directory / f"{configuration}.json"
    else:
        candidates = sorted(directory.glob("*.json")) if directory.exists() else []
        if not candidates:
            raise FileNotFoundError(
                f"no assessment JSON under {directory}; run a campaign first"
            )
        path = candidates[0]
    payload = _load_json(path)

    records = payload["requirements"]
    ids = sorted(records)
    verdicts = [str(records[rid].get("verdict", "not_evaluated")) for rid in ids]
    colours = [style.VERDICT_COLOURS.get(v, style.INK_MUTED) for v in verdicts]

    style.use()
    height = max(3.0, 0.28 * len(ids) + 1.4)
    fig, ax = style.figure(figsize=(7.4, height))
    positions = range(len(ids))
    ax.barh(list(positions), [1.0] * len(ids), color=colours, height=0.62, zorder=3)
    ax.set_yticks(list(positions))
    ax.set_yticklabels(ids)
    ax.set_ylim(len(ids) - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_xlim(0.0, 1.0)
    for spine in ("bottom", "left"):
        ax.spines[spine].set_visible(False)

    for y, rid in enumerate(ids):
        record = records[rid]
        achieved = record.get("achieved")
        label = record.get("verdict", "not_evaluated")
        if achieved is not None:
            label = f"{label}  ({_short(achieved)} {record.get('unit', '')})".strip()
        ax.text(0.02, y, label, va="center", ha="left", fontsize=9, color=style.SURFACE)

    summary = payload["summary"]
    ax.set_title(f"Requirements: {payload['configuration']}")
    style.caption(
        fig,
        f"{summary['goal']} goal, {summary['threshold']} threshold, "
        f"{summary['fail']} fail, {summary['not_evaluated']} not evaluated.  "
        "Not evaluated means the metric was not computed by this run, not that "
        "the requirement was met.  Source: " + str(path.name),
    )
    return style.savefig(fig, Path(out) / f"requirements_{payload['configuration']}")


def _load_json(path: Path) -> dict[str, Any]:
    if not Path(path).exists():
        raise FileNotFoundError(str(path))
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _short(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if value == 0.0:
            return "0"
        return f"{value:.3e}" if abs(value) >= 1.0e5 or abs(value) < 1.0e-2 else f"{value:.4g}"
    return str(value)


# ---------------------------------------------------------------------------
# Declared, not yet written
# ---------------------------------------------------------------------------
#
# Each row names the study that has to run before the figure can be drawn.  The
# inventory is here rather than in a document so that `tss figure <name>` can
# answer the question directly.

_declare(
    "band_edge_coverage",
    "Profile coverage against lower band edge over the mission envelope, with "
    "the 95 and 99 per cent points that INS-01's threshold and goal come from.",
    ("results/min_frequency/*.json",),
    "tss study min_frequency",
)
_declare(
    "link_budget_snr",
    "SNR per range cell against sounding frequency for each configuration, "
    "against the INS-10 threshold and goal.  The figure C1 rests on.",
    ("results/link_budget/*.json",),
    "tss study link_budget",
)
_declare(
    "sky_noise_reference",
    "External noise density against frequency, showing where the Franklin and "
    "MacLean measurement ends at about 0.6 MHz and the extrapolation begins.",
    ("results/sky_noise/*.json",),
    "tss study sky_noise",
)
_declare(
    "antenna_scale",
    "Radiation efficiency and achieved SNR against deployed antenna length, "
    "which is the axis the 6U volume constraint binds on.",
    ("results/antenna_scale/*.json",),
    "tss study antenna_scale",
)
_declare(
    "isis1_validation",
    "Retrieved profiles against the hand-scaled ISIS-1 archive, with the "
    "fractional density difference distribution.  The V4 evidence for C1.",
    ("results/isis_validation/*.json",),
    "tss study isis_validation",
)
_declare(
    "baseline_ionogram",
    "A simulated ionogram from the baseline chain, with the extracted trace "
    "overlaid on the range-frequency power surface.",
    ("results/<campaign>/ionogram.npz",),
    "tss run baseline, once the pipeline stores intermediates",
)
_declare(
    "doppler_range_bias",
    "Apparent-range bias against sounding frequency for a range of "
    "line-of-sight velocities, and the residual after the up and down sweep "
    "pairing.  The figure C3 rests on.",
    ("results/c3_doppler/*.json",),
    "tss run c3_doppler, blocked on WP5",
)
_declare(
    "leakage_residual",
    "De-chirped leakage residual against isolation, with the INS-24 margin "
    "over the weakest echo INS-10 admits.  The figure C4 rests on.",
    ("results/leakage/*.json",),
    "WP5 processing and WP9 bench measurement",
)
_declare(
    "feasibility_surface",
    "Feasible region over antenna length and transmit power, per "
    "configuration, annotated with the requirement binding each boundary "
    "segment.",
    ("results/trade_space/*.json",),
    "tss run trade_space, blocked on WP7",
)
_declare(
    "retrieval_error_budget",
    "Height and density error against altitude, with the random and the bias "
    "terms shown separately rather than combined.",
    ("results/<campaign>/error_budget.json",),
    "WP6",
)
