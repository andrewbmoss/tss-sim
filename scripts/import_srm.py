"""Import the requirement catalogue from the SRM workbook.

    python scripts/import_srm.py PATH/TO/TSS_SRM_RevG.xlsx
    python scripts/import_srm.py PATH/TO/TSS_SRM_RevG.xlsx --check

The SRM workbook is the requirement set of record.  This script reads its three
requirement sheets (Measurement Reqs, Instrument & Mission, Software) and writes
``config/requirements/srm_catalog.yaml``: every requirement ID with its text
exactly as the workbook states it, the verification methods, and the SHA-256 of
the workbook it came from.  The catalogue is generated and is never edited by
hand.  Rerunning the script after the workbook changes and reading the
``git diff`` of the catalogue is how a new revision is reviewed, which is the
"diff review rather than a retype" that ``tss.requirements`` promises.

The catalogue is not what L5 scores against.  ``requirements.yaml`` carries the
numbers a simulation can score, with their sense, units and sources, and it
names the catalogue it belongs to.  :meth:`tss.requirements.RequirementSet.validate`
then checks that every catalogue ID is either scored there or listed under
``not_scored`` with a reason, so a requirement cannot drop out of the scored set
without anyone saying why.

The workbook is read with the standard library (``zipfile`` and
``xml.etree``) rather than openpyxl, so the script runs in the ``tss``
environment without a new dependency.  It reads the cached cell values that
Excel stores beside each formula, which is what ``data_only=True`` would give.

``--check`` writes nothing.  It exits 1 if the catalogue on disk differs from
what the workbook would generate, so a test or CI job can tell that the two
have drifted apart.

Exit status: 0 written or up to date, 1 out of date under ``--check``, 2 the
workbook could not be read or holds no requirements.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree

import yaml

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO / "config" / "requirements" / "srm_catalog.yaml"

_NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "rel": "http://schemas.openxmlformats.org/package/2006/relationships",
}

#: The sheets that hold requirements, and the meaning of each column after ID.
#: The header row of each sheet is checked against these names, so a workbook
#: whose columns have moved fails loudly instead of filing text under the
#: wrong key.
SHEETS: dict[str, tuple[tuple[str, str], ...]] = {
    "Measurement Reqs": (
        ("Requirement", "requirement"),
        ("Threshold", "threshold"),
        ("Goal", "goal"),
        ("Basis / source", "basis"),
        ("Verification", "verification"),
    ),
    "Instrument & Mission": (
        ("Requirement", "requirement"),
        ("Value", "value"),
        ("Flows from", "flows_from"),
        ("Current point design / status", "status"),
        ("Verification", "verification"),
    ),
    "Software (tss-sim)": (
        ("Requirement", "requirement"),
        ("Value", "value"),
        ("Flows from", "flows_from"),
        ("Verification", "verification"),
    ),
}

_ID = re.compile(r"^(SCI|INS|MIS|SW)-\d{2}[a-z]?$")
#: Verification methods named in the SRM: analysis, simulation, archived or
#: external data, bench, inspection, test, demonstration, in orbit.
_METHODS = set("ASRBITDO")
_METHOD_PREFIX = re.compile(r"^\s*([A-Z](?:\s*,\s*[A-Z])*)\s*\.")
_REVISION = re.compile(r"\bRev\s+([A-Z])\b")


# ---------------------------------------------------------------------------
# A minimal reader for the cached values of an .xlsx workbook
# ---------------------------------------------------------------------------


def _column_index(ref: str) -> int:
    """Zero-based column of a cell reference such as ``C14``."""
    letters = "".join(ch for ch in ref if ch.isalpha())
    index = 0
    for ch in letters:
        index = index * 26 + (ord(ch.upper()) - ord("A") + 1)
    return index - 1


def _text_of(node: ElementTree.Element) -> str:
    """Concatenated text of every ``<t>`` under a node, rich-text runs included."""
    return "".join(t.text or "" for t in node.iter(f"{{{_NS['m']}}}t"))


def read_workbook(path: Path) -> dict[str, list[list[str | None]]]:
    """Every sheet as a list of rows of cell strings, ``None`` for an empty cell."""
    with zipfile.ZipFile(path) as zf:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ElementTree.fromstring(zf.read("xl/sharedStrings.xml"))
            shared = [_text_of(si) for si in root.findall("m:si", _NS)]

        rels = ElementTree.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        target = {
            rel.get("Id"): rel.get("Target")
            for rel in rels.findall("rel:Relationship", _NS)
        }
        book = ElementTree.fromstring(zf.read("xl/workbook.xml"))

        sheets: dict[str, list[list[str | None]]] = {}
        for sheet in book.findall("m:sheets/m:sheet", _NS):
            name = sheet.get("name") or ""
            rid = sheet.get(f"{{{_NS['r']}}}id")
            part = str(target[rid]).lstrip("/")
            part = part if part.startswith("xl/") else f"xl/{part}"
            root = ElementTree.fromstring(zf.read(part))
            rows: list[list[str | None]] = []
            for row in root.findall("m:sheetData/m:row", _NS):
                cells: dict[int, str | None] = {}
                for cell in row.findall("m:c", _NS):
                    col = _column_index(cell.get("r") or "A1")
                    kind = cell.get("t")
                    if kind == "inlineStr":
                        node = cell.find("m:is", _NS)
                        value = _text_of(node) if node is not None else None
                    else:
                        v = cell.find("m:v", _NS)
                        raw = None if v is None else v.text
                        if raw is None:
                            value = None
                        elif kind == "s":
                            value = shared[int(raw)]
                        elif kind == "b":
                            value = "TRUE" if raw == "1" else "FALSE"
                        else:
                            value = raw
                    cells[col] = value
                width = max(cells) + 1 if cells else 0
                rows.append([cells.get(i) for i in range(width)])
            sheets[name] = rows
    return sheets


# ---------------------------------------------------------------------------
# From sheets to a catalogue
# ---------------------------------------------------------------------------


def _clean(text: str | None) -> str:
    """Collapse internal whitespace and strip; keep the wording exactly."""
    if text is None:
        return ""
    return re.sub(r"\s+", " ", str(text)).strip()


def verification_methods(text: str) -> list[str]:
    """The method letters that open a Verification cell, e.g. ``S, R.``."""
    match = _METHOD_PREFIX.match(text or "")
    if not match:
        return []
    letters = [ch.strip() for ch in match.group(1).split(",")]
    return [ch for ch in letters if ch in _METHODS]


def _header_row(rows: list[list[str | None]], sheet: str) -> int:
    for i, row in enumerate(rows):
        if row and _clean(row[0]) == "ID":
            return i
    raise ValueError(f"sheet {sheet!r} has no header row starting with 'ID'")


def catalogue_from_sheets(sheets: dict[str, list[list[str | None]]]) -> dict[str, dict]:
    """Requirement ID to its record, in workbook order."""
    out: dict[str, dict] = {}
    for sheet, columns in SHEETS.items():
        if sheet not in sheets:
            raise ValueError(f"the workbook has no sheet named {sheet!r}")
        rows = sheets[sheet]
        h = _header_row(rows, sheet)
        header = [_clean(c) for c in rows[h]]
        expected = ["ID", *[label for label, _ in columns]]
        if header[: len(expected)] != expected:
            raise ValueError(
                f"sheet {sheet!r}: header is {header[:len(expected)]}, expected {expected}"
            )
        for row in rows[h + 1:]:
            rid = _clean(row[0] if row else None)
            if not _ID.match(rid):
                continue
            if rid in out:
                raise ValueError(f"{rid} appears twice in the workbook")
            record: dict[str, object] = {"sheet": sheet, "family": rid.split("-")[0]}
            for offset, (_, key) in enumerate(columns, start=1):
                record[key] = _clean(row[offset] if offset < len(row) else None)
            record["verification_methods"] = verification_methods(str(record["verification"]))
            out[rid] = record
    if not out:
        raise ValueError("no requirement IDs found")
    return out


def _revision(sheets: dict[str, list[list[str | None]]], path: Path) -> str:
    """The revision letter, from the SRM sheet's subtitle or else the file name."""
    for row in sheets.get("SRM", [])[:4]:
        for cell in row:
            match = _REVISION.search(_clean(cell))
            if match:
                return match.group(1)
    match = re.search(r"Rev([A-Z])\b", path.stem)
    return match.group(1) if match else ""


def build_catalogue(path: Path) -> dict:
    """The catalogue document for one workbook."""
    sheets = read_workbook(path)
    requirements = catalogue_from_sheets(sheets)
    return {
        "generated_by": "scripts/import_srm.py",
        "source_workbook": path.name,
        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "revision": _revision(sheets, path),
        "counts": dict(Counter(r["family"] for r in requirements.values())),
        "requirements": requirements,
    }


_HEADER = """\
# GENERATED FILE - do not edit by hand.
#
# The requirement catalogue, transcribed from the SRM workbook by
#     python scripts/import_srm.py <workbook>
# Every ID in the workbook appears here with its text as the workbook states it.
# requirements.yaml names this file under `catalog:`, and `tss requirements`
# fails unless every ID here is either scored there or listed under
# `not_scored` with a reason.  To take in a new revision, rerun the script and
# review the diff of this file.
"""


def render(catalogue: dict) -> str:
    body = yaml.safe_dump(
        catalogue, sort_keys=False, allow_unicode=True, width=100, default_flow_style=False
    )
    return _HEADER + "\n" + body


def _summarise(old: dict | None, new: dict, stream) -> None:
    new_req = new["requirements"]
    counts = ", ".join(f"{k} {v}" for k, v in new["counts"].items())
    print(f"Rev {new['revision']} from {new['source_workbook']}: {len(new_req)} requirements ({counts})",
          file=stream)
    if not old:
        return
    old_req = old.get("requirements", {})
    added = sorted(set(new_req) - set(old_req))
    removed = sorted(set(old_req) - set(new_req))
    changed = sorted(r for r in set(new_req) & set(old_req) if new_req[r] != old_req[r])
    for label, ids in (("added", added), ("removed", removed), ("text changed", changed)):
        if ids:
            print(f"  {label}: {', '.join(ids)}", file=stream)
    if not (added or removed or changed):
        print("  no requirement changed since the catalogue on disk", file=stream)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("workbook", type=Path, help="the SRM .xlsx")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--check", action="store_true",
                        help="write nothing; exit 1 if the catalogue on disk is out of date")
    args = parser.parse_args(argv)

    try:
        catalogue = build_catalogue(args.workbook)
    except (OSError, KeyError, ValueError, zipfile.BadZipFile, ElementTree.ParseError) as exc:
        print(f"cannot read {args.workbook}: {exc}", file=sys.stderr)
        return 2

    text = render(catalogue)
    old_text = args.out.read_text(encoding="utf-8") if args.out.exists() else None
    old = yaml.safe_load(old_text) if old_text else None
    _summarise(old, catalogue, sys.stdout)

    if args.check:
        if old_text != text:
            print(f"{args.out} is out of date; rerun without --check", file=sys.stdout)
            return 1
        print(f"{args.out} is up to date", file=sys.stdout)
        return 0

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text, encoding="utf-8", newline="\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
