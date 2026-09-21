"""Layer: support
WP: WP8
Claims: none
Status: implemented
Requirements: -

One place to change how every figure looks.

Adapted from ``src/nisc/figures/style.py`` in the author's Chartier
verification work, with one structural change: matplotlib is imported inside
the functions that need it and never at module scope.  ``import tss`` imports
every layer package in order to populate the tier registry, and
``CONTRIBUTING.md`` requires that import to stay cheap and to work on a machine
with none of the optional dependencies installed.  matplotlib is optional, so
importing it here would break that on a headless machine with no plotting
stack.  ``tests/test_import_cost.py`` guards the rule.

The palette is the validated categorical default, assigned in fixed slot order
and never cycled.  Fixed slot order matters for a dissertation: the same series
should be the same colour in every figure it appears in, and a cycling palette
gives it whatever colour its position in that particular plot happened to
allocate.

Sequential ramps are single hue and are for magnitude encodings: signal loss,
electron density, reflection height.  They are built inside
:func:`sequential` rather than at module scope, because constructing a
:class:`matplotlib.colors.LinearSegmentedColormap` needs matplotlib.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

__all__ = [
    "SERIES",
    "SURFACE",
    "INK",
    "INK_2",
    "INK_MUTED",
    "GRID",
    "RC",
    "VERDICT_COLOURS",
    "use",
    "sequential",
    "grid",
    "caption",
    "annotate",
    "savefig",
    "figure",
]

#: Categorical slots, in fixed order.  Assigned by slot, never cycled.
SERIES: tuple[str, ...] = (
    "#2a78d6",
    "#eb6834",
    "#1baf7a",
    "#eda100",
    "#e87ba4",
    "#008300",
    "#4a3aa7",
    "#e34948",
)

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8983"
GRID = "#e4e3de"

#: Verdict colours, used by the requirement scorecard and the verification
#: summary so that the two read the same way.  ``not_evaluated`` is muted grey
#: on purpose: it is neither good news nor bad news and should not look like
#: either.
VERDICT_COLOURS: dict[str, str] = {
    "goal": "#1baf7a",
    "threshold": "#2a78d6",
    "pass": "#1baf7a",
    "partial": "#eda100",
    "fail": "#e34948",
    "skipped": "#eb6834",
    "not_run": "#8a8983",
    "not_evaluated": "#8a8983",
}

#: rcParams applied by :func:`use`.  Kept as a plain dict so that it can be
#: read and diffed without importing matplotlib.
RC: dict[str, Any] = {
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.size": 10.5,
    "font.family": "DejaVu Sans",
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "axes.titleweight": "semibold",
    "axes.titlelocation": "left",
    "axes.titlepad": 10,
    "axes.edgecolor": GRID,
    "axes.labelcolor": INK_2,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "text.color": INK,
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "grid.color": GRID,
    "grid.linewidth": 0.7,
    "legend.frameon": False,
    "legend.fontsize": 9.5,
    "lines.linewidth": 2.0,
    "lines.markersize": 5,
    "figure.dpi": 130,
    "savefig.dpi": 190,
    "savefig.bbox": "tight",
}

_SEQUENTIAL_STOPS: dict[str, tuple[str, ...]] = {
    "blue": ("#e8f0fb", "#a9c8ee", "#5f9ae2", "#2a78d6", "#18509a", "#0d2f5c"),
    "warm": ("#fdf0e8", "#f6c4a6", "#f0956a", "#eb6834", "#a8451f", "#5e2611"),
}


def use(*, backend: str = "Agg") -> None:
    """Apply the style globally and select a non-interactive backend.

    The Agg backend is selected by default because every figure in this package
    is written to a file rather than shown, and a figure regenerated on a
    headless machine must not try to open a window.  Pass ``backend=None`` to
    leave the current backend alone.
    """
    import matplotlib as mpl

    if backend:
        mpl.use(backend)
    mpl.rcParams.update(RC)


def sequential(name: str = "blue"):
    """Return a single-hue sequential colormap for a magnitude encoding."""
    from matplotlib.colors import LinearSegmentedColormap

    key = str(name).lower()
    if key not in _SEQUENTIAL_STOPS:
        raise ValueError(
            f"unknown sequential ramp {name!r}; expected one of "
            f"{tuple(_SEQUENTIAL_STOPS)}"
        )
    return LinearSegmentedColormap.from_list(f"seq_{key}", list(_SEQUENTIAL_STOPS[key]))


def figure(*args: Any, **kwargs: Any):
    """``matplotlib.pyplot.subplots`` with the style already applied."""
    use()
    import matplotlib.pyplot as plt

    return plt.subplots(*args, **kwargs)


def grid(ax, axis: str = "both") -> None:
    """Recessive grid, behind the data."""
    ax.grid(True, axis=axis, alpha=0.85, zorder=0)
    ax.set_axisbelow(True)


def caption(fig, text: str, y: float = -0.02) -> None:
    """A small note under a figure, in muted ink.

    Used for the thing a reader has to know in order to not misread the plot:
    which tier produced it, what was assumed, what is missing.
    """
    fig.text(0.0, y, text, ha="left", va="top", fontsize=8.5, color=INK_MUTED, wrap=True)


def annotate(ax, x: float, y: float, text: str, *, color: str = INK_2,
             dx: float = 6, dy: float = 6, **kwargs: Any) -> None:
    """Direct label.  Used selectively, never on every point."""
    ax.annotate(
        text,
        (x, y),
        textcoords="offset points",
        xytext=(dx, dy),
        fontsize=9,
        color=color,
        **kwargs,
    )


def savefig(fig, path: str | Path, formats: Sequence[str] = ("png",)) -> list[Path]:
    """Write a figure in each requested format and close it."""
    import matplotlib.pyplot as plt

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for ext in formats:
        target = out_path.with_suffix(f".{ext}")
        fig.savefig(target)
        written.append(target)
    plt.close(fig)
    return written
