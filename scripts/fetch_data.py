#!/usr/bin/env python3
"""Fetch external data sets without installing the package.

A thin wrapper over :func:`tss.data.fetch`.  ``tss fetch`` is the supported
entry point; this exists so that a reader who has cloned the repository and not
run ``pip install -e .`` can still pull the inputs.

    python scripts/fetch_data.py list                 # what is available
    python scripts/fetch_data.py show spdf_isis_ne    # licence, citation, notes
    python scripts/fetch_data.py spdf_isis_ne         # fetch one source
    python scripts/fetch_data.py all --limit 1        # fetch everything, sparsely
    python scripts/fetch_data.py models               # probe local model installs

Nothing is downloaded by ``list``, ``show`` or ``models``.  Fetching honours the
cache unless ``--force`` is given, and every download writes a ``.meta.json``
sidecar recording the URL, the retrieval time and the SHA-256.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fetch_data.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "source",
        help="source name, or 'list', 'show', 'all', or 'models'",
    )
    p.add_argument("name", nargs="?", help="the source to describe, with 'show'")
    p.add_argument("--force", action="store_true", help="ignore the cache")
    p.add_argument("--limit", type=int, default=None,
                   help="cap how much is fetched, where the source honours it")
    p.add_argument("--json", action="store_true", help="machine-readable output")
    p.add_argument("--no-trial", action="store_true",
                   help="with 'models', skip the trial calls and only look for files")
    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    # Imported here so that --help works even if a dependency is missing.
    from tss.data import SOURCES, describe_sources, fetch
    from tss.data.base import FetchError, UnknownSource

    if args.source == "list":
        if args.json:
            print(json.dumps(describe_sources(), indent=2))
            return 0
        for name, src in sorted(SOURCES.items()):
            wps = ", ".join(src.work_packages) or "-"
            print(f"{name:22s} [{wps:20s}] {src.description}")
        return 0

    if args.source == "show":
        if not args.name:
            print("show needs a source name", file=sys.stderr)
            return 2
        try:
            src = SOURCES[args.name]
        except KeyError:
            print(f"no source named {args.name!r}; try 'list'", file=sys.stderr)
            return 2
        if args.json:
            print(json.dumps(src.describe(), indent=2))
            return 0
        print(f"{src.name}\n  {src.description}")
        print(f"  url        {src.url}")
        print(f"  licence    {src.licence}")
        print(f"  citation   {src.citation}")
        print(f"  cadence    {src.cadence}")
        print(f"  packages   {', '.join(src.work_packages) or '-'}")
        if src.acknowledgement:
            print(f"  {'acknowledge':<10s} {src.acknowledgement}")
        if src.notes:
            print(f"  notes      {src.notes}")
        return 0

    if args.source == "models":
        from tss.data.local_models import probe

        status = probe(trial=not args.no_trial)
        if args.json:
            print(json.dumps(status, indent=2))
            return 0
        for name, rec in status.items():
            if name.startswith("_"):
                continue
            mark = "yes" if rec["available"] else "no "
            trial = "trial" if rec["verified_by_trial"] else "looked"
            version = rec["version"] or "-"
            print(f"{name:10s} {mark}  {trial:6s} {version:12s} {rec['path']}")
            if not rec["available"]:
                print(f"           {rec['reason']}")
        return 0

    names = sorted(SOURCES) if args.source == "all" else [args.source]

    rc = 0
    for name in names:
        try:
            paths = fetch(name, force=args.force, limit=args.limit)
        except UnknownSource as exc:
            print(str(exc), file=sys.stderr)
            return 2
        except FetchError as exc:
            print(f"{name}: {exc}", file=sys.stderr)
            rc = 1
            continue
        for path in paths:
            print(path)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
