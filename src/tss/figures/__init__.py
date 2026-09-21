"""Layer: support
WP: WP8
Claims: C1, C2, C3, C4, C5
Status: partial
Requirements: -

Figure regeneration.

``tss figure <name>`` regenerates one figure from a stored result.  Nothing in
this package computes a result; a figure function loads what a study or a
campaign wrote and draws it.  That is the reproducibility promise of Section
4.2: every figure in the dissertation comes from a stored configuration and a
stored result, and rerunning the command reproduces it.

matplotlib is an optional dependency and is imported inside the drawing
functions, never at module scope.  ``import tss`` imports this package in order
to keep the command line working, and it must stay importable on a machine with
no plotting stack.

:func:`regenerate` returns a process exit status rather than raising, because
``tss.cli`` returns it straight to the shell:

    0   the figure was drawn and written
    2   the figure is declared but not yet written; the message says which
        study has to land first
    3   the figure exists but its input result is not on disk; the message
        names the file and the command that writes it
    4   no figure of that name is registered
    5   matplotlib is not installed

Those are distinguished rather than collapsed into a single failure because
they call for four different actions, and a script that reruns the figure set
after a campaign should be able to tell "not written yet" from "you have not
run the study".
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

from tss.figures import figures, style  # noqa: F401
from tss.figures.figures import (  # noqa: F401
    FIGURES,
    Figure,
    available_figures,
    default_figure_dir,
    register_figure,
)

__all__ = [
    "regenerate",
    "available_figures",
    "FIGURES",
    "Figure",
    "register_figure",
    "default_figure_dir",
    "style",
]


def regenerate(
    name: str,
    out: str | Path | None = None,
    *,
    stream=None,
    **kwargs: Any,
) -> int:
    """Regenerate one named figure.  Called by ``tss figure``.

    ``name`` of ``"list"`` prints the registry and returns 0, which is how a
    reader finds out what exists without reading this module.
    """
    out_stream = stream if stream is not None else sys.stdout

    if name == "list":
        for figure_name, record in available_figures().items():
            marker = "  " if record["status"] == "implemented" else "* "
            print(f"{marker}{figure_name:26s} {record['description']}", file=out_stream)
        print(
            "\n* declared but not yet written; the entry names the study that has "
            "to run first.",
            file=out_stream,
        )
        return 0

    figure = FIGURES.get(name)
    if figure is None:
        known = ", ".join(sorted(FIGURES)) or "nothing"
        print(f"no figure named {name!r}; registered: {known}", file=out_stream)
        return 4

    if figure.status != "implemented" or figure.draw is None:
        print(
            f"figure {name!r} is declared but not yet written.\n"
            f"  what it would show: {figure.description}\n"
            f"  it reads:           {', '.join(figure.inputs) or '(unstated)'}\n"
            f"  produced by:        {figure.produced_by or '(unstated)'}",
            file=out_stream,
        )
        return 2

    target = Path(out) if out is not None else default_figure_dir()
    try:
        written = figure.draw(target, **kwargs)
    except ModuleNotFoundError as exc:
        print(
            f"cannot draw {name!r}: {exc}.  matplotlib is an optional dependency; "
            "install it with `pip install 'tss-sim[figures]'`.",
            file=out_stream,
        )
        return 5
    except FileNotFoundError as exc:
        print(
            f"cannot draw {name!r}: the input {exc} is not on disk.\n"
            f"  it is written by: {figure.produced_by or '(unstated)'}",
            file=out_stream,
        )
        return 3

    for path in written:
        print(str(path), file=out_stream)
    return 0
