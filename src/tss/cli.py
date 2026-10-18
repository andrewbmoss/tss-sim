"""Command line entry point.

``tss`` is the single driver for the campaign.  Every figure and table in the
dissertation is produced by one of these subcommands from a stored
configuration, which is the reproducibility promise of Section 4.2 made
operational.

    tss tiers                       list registered tier implementations
    tss requirements                show the requirement set and check it
    tss config <campaign>           resolve and print a campaign configuration
    tss fetch <source>              download and cache an external data set
    tss study <name>                run a study and write results plus manifest
    tss run <campaign>              run a campaign through the full pipeline
    tss verify [--level V1..V4]     run the verification suite
    tss figure <name>               regenerate one figure from stored results
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Sequence
from pathlib import Path

__all__ = ["main", "build_parser"]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="tss", description=__doc__.splitlines()[0])
    p.add_argument("--version", action="store_true", help="print the version and exit")
    sub = p.add_subparsers(dest="command")

    s = sub.add_parser("tiers", help="list registered tier implementations")
    s.add_argument("--layer", choices=["l0", "l1", "l2", "l3", "l4", "l5"])

    s = sub.add_parser("requirements", help="show and validate the requirement set")
    s.add_argument("--set", dest="req_set", default="requirements/requirements")
    s.add_argument("--claim", help="show only requirements serving this claim, e.g. C3")
    s.add_argument("--json", action="store_true")

    s = sub.add_parser("config", help="resolve and print a campaign configuration")
    s.add_argument("campaign")
    s.add_argument("--set", dest="overrides", action="append", default=[],
                   metavar="a.b=value")

    s = sub.add_parser("fetch", help="download and cache an external data set")
    s.add_argument("source", help="source name, or 'list' to show what is available")
    s.add_argument("--force", action="store_true", help="ignore the cache")
    s.add_argument("--limit", type=int, default=None)

    s = sub.add_parser("study", help="run a study")
    s.add_argument("name")
    s.add_argument("--out", type=Path, default=None)
    s.add_argument("--set", dest="overrides", action="append", default=[])
    s.add_argument("--strict", action="store_true",
                   help="refuse to run from a dirty working tree")

    s = sub.add_parser("run", help="run a campaign through the full pipeline")
    s.add_argument("campaign")
    s.add_argument("--out", type=Path, default=None)
    s.add_argument("--set", dest="overrides", action="append", default=[])
    s.add_argument("--strict", action="store_true")

    s = sub.add_parser("verify", help="run the verification suite")
    s.add_argument("--level", choices=["V1", "V2", "V3", "V4", "all"], default="all")
    s.add_argument("--write", action="store_true", help="write the verification report")

    s = sub.add_parser("figure", help="regenerate a figure")
    s.add_argument("name")
    s.add_argument("--out", type=Path, default=None)

    return p


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.version:
        from tss import __version__
        print(__version__)
        return 0

    if args.command is None:
        build_parser().print_help()
        return 1

    if args.command == "tiers":
        from tss.registry import available
        for layer, names in available(args.layer).items():
            print(f"{layer}: {', '.join(names) or '(none registered)'}")
        return 0

    if args.command == "requirements":
        from tss.requirements import load_requirements
        reqs = load_requirements(args.req_set)
        problems = reqs.validate()
        selected = reqs.for_claim(args.claim) if args.claim else list(reqs.values())
        if args.json:
            print(json.dumps([r.__dict__ for r in selected], indent=2, default=str))
        else:
            print(f"requirement set {reqs.revision or '(no revision)'}: {len(reqs)} entries")
            for r in sorted(selected, key=lambda r: r.id):
                claims = ",".join(r.claims) or "-"
                print(f"  {r.id:8s} {r.title[:44]:44s} thr={r.threshold!s:>16s} "
                      f"goal={r.goal!s:>16s} {r.unit:8s} [{claims}]")
            cov = reqs.coverage()
            if cov["n_catalog"]:
                held = ", ".join(
                    f"{n} {cat}" for cat, n in cov["not_scored_by_category"].items() if n
                )
                print(f"\ncatalogue Rev {cov['catalog_revision']} "
                      f"({cov['catalog_workbook']}): {cov['n_catalog']} requirements, "
                      f"{cov['n_scored']} scored, {cov['n_not_scored']} not scored ({held})")
        if problems:
            print(f"\n{len(problems)} problem(s):", file=sys.stderr)
            for msg in problems:
                print(f"  {msg}", file=sys.stderr)
            return 2
        return 0

    if args.command == "config":
        from tss.config import load_campaign
        campaign = load_campaign(args.campaign, args.overrides)
        print(json.dumps(campaign.config, indent=2, default=str))
        return 0

    if args.command == "fetch":
        from tss.data import SOURCES, fetch
        if args.source == "list":
            for name, src in sorted(SOURCES.items()):
                print(f"{name:24s} {src.description}")
            return 0
        paths = fetch(args.source, force=args.force, limit=args.limit)
        for p in paths:
            print(p)
        return 0

    if args.command == "study":
        from tss.studies import run_study
        return run_study(args.name, out=args.out, overrides=args.overrides,
                         strict=args.strict)

    if args.command == "run":
        from tss.studies import run_campaign
        return run_campaign(args.campaign, out=args.out, overrides=args.overrides,
                            strict=args.strict)

    if args.command == "verify":
        from tss.verify import run_verification
        return run_verification(level=args.level, write=args.write)

    if args.command == "figure":
        from tss.figures import regenerate
        return regenerate(args.name, out=args.out)

    raise AssertionError(f"unhandled command {args.command!r}")


if __name__ == "__main__":
    raise SystemExit(main())
