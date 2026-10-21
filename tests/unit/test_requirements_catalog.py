"""The Rev G catalogue, its import, and the rule that nothing drops out silently.

``config/requirements/srm_catalog.yaml`` is generated from the SRM workbook by
``scripts/import_srm.py``.  ``requirements.yaml`` scores part of it and lists
the rest under ``not_scored`` with a reason.  These tests hold three things:
the catalogue is the whole Rev G set; every catalogue ID is accounted for; and
the importer reads a workbook correctly, including the real one when it is on
this machine, so the catalogue cannot drift from the workbook unnoticed.
"""

from __future__ import annotations

import importlib.util
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import pytest

from tss.config import repo_root
from tss.requirements import (
    NOT_SCORED_CATEGORIES,
    NotScored,
    Requirement,
    RequirementSet,
    load_requirements,
)

REPO = repo_root()
REAL_WORKBOOK = REPO.parent / "Dissertation" / "TSS_SRM_RevG.xlsx"

#: The requirements of Table 3.1 of the proposal (Rev D).  Every one must be
#: scored; SCI-15 was missing from requirements.yaml until 2026-10-17.
TABLE_3_1 = (
    "INS-01a", "INS-01b", "INS-02", "INS-04", "INS-05", "INS-07", "INS-08", "INS-09",
    "INS-10", "INS-11", "INS-14", "INS-15", "INS-18", "INS-23", "INS-24", "INS-26",
    "SCI-05", "SCI-07", "SCI-09", "SCI-15", "SCI-23",
)


def _import_srm():
    spec = importlib.util.spec_from_file_location("import_srm", REPO / "scripts" / "import_srm.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# The catalogue and the scored set
# ---------------------------------------------------------------------------


def test_the_catalogue_is_the_full_rev_g_set():
    reqs = load_requirements()
    cov = reqs.coverage()
    assert cov["catalog_revision"] == "G"
    # 72 since the 2026-10-20 amendment (WP0 closure) added INS-28 and MIS-11.
    assert cov["n_catalog"] == 72
    families = [rid.split("-")[0] for rid in reqs.catalog]
    assert {f: families.count(f) for f in set(families)} == {
        "SCI": 28, "INS": 29, "MIS": 11, "SW": 4,
    }


def test_every_catalogue_requirement_is_scored_or_says_why():
    reqs = load_requirements()
    assert reqs.validate() == []
    cov = reqs.coverage()
    assert cov["unaccounted"] == []
    assert cov["n_scored"] + cov["n_not_scored"] == cov["n_catalog"]
    assert sum(cov["not_scored_by_category"].values()) == cov["n_not_scored"]


def test_every_table_3_1_requirement_is_scored():
    assert set(TABLE_3_1) <= set(load_requirements())


def test_every_not_scored_reason_is_in_a_known_category():
    for entry in load_requirements().not_scored.values():
        assert entry.category in NOT_SCORED_CATEGORIES
        assert entry.reason.strip()


def _set(*, catalog, scored=(), not_scored=()):
    reqs = [
        Requirement(id=rid, title=rid, sense="max", threshold=1.0, source="test")
        for rid in scored
    ]
    return RequirementSet(
        reqs, catalog={rid: {"requirement": rid} for rid in catalog}, not_scored=not_scored
    )


def test_a_requirement_that_drops_out_fails_validation():
    rs = _set(catalog=["A-01", "A-02"], scored=["A-01"])
    assert any("A-02" in p and "neither scored" in p for p in rs.validate())


def test_scored_and_not_scored_at_once_fails_validation():
    rs = _set(
        catalog=["A-01"], scored=["A-01"],
        not_scored=[NotScored("A-01", "analysis", "because")],
    )
    assert any("both scored" in p for p in rs.validate())


def test_a_scored_requirement_missing_from_the_catalogue_fails_validation():
    rs = _set(catalog=["A-01"], scored=["A-01", "A-09"])
    assert any("A-09" in p and "not in the catalogue" in p for p in rs.validate())


@pytest.mark.parametrize(
    "entry,fragment",
    [
        (NotScored("A-02", "nonsense", "x"), "category"),
        (NotScored("A-02", "analysis", "  "), "without a reason"),
        (NotScored("A-02", "flowed_down", "x"), "must name what it is scored through"),
        (NotScored("A-02", "flowed_down", "x", scored_through=("A-07",)), "not a scored"),
        (NotScored("A-02", "pending", "x"), "work package"),
    ],
)
def test_a_malformed_not_scored_entry_fails_validation(entry, fragment):
    rs = _set(catalog=["A-01", "A-02"], scored=["A-01"], not_scored=[entry])
    assert any(fragment in p for p in rs.validate()), rs.validate()


def test_a_set_with_no_catalogue_still_validates():
    """A requirement file that names no catalogue, such as a study's own, keeps
    working as it did before the catalogue existed."""
    rs = RequirementSet([Requirement(id="T-01", title="t", sense="max", threshold=1.0,
                                     source="test")])
    assert rs.validate() == []
    assert rs.coverage()["n_catalog"] == 0


# ---------------------------------------------------------------------------
# The report says what is not scored, and why
# ---------------------------------------------------------------------------


def test_the_report_lists_every_not_scored_requirement_with_its_reason():
    from tss.interfaces import DensityProfile, Ionogram, IonosphereState
    from tss.l5_assessment import StandardAssessment, render_json, render_text

    alt = np.linspace(2.0e5, 5.0e5, 50)
    truth = IonosphereState(alt_m=alt, ne=1.0e11 * np.exp(-(alt - 2.0e5) / 1.0e5))
    ion = Ionogram(freq_hz=np.array([1.0e6, 2.0e6]),
                   virtual_range_m=np.array([[1.0e5], [2.0e5]]),
                   snr_db=np.array([[20.0], [20.0]]), modes=("unresolved",))
    profile = DensityProfile(alt_m=alt, ne=truth.ne)
    assessment = StandardAssessment().score(profile, truth, ion, "catalogue")
    text = render_text(assessment)
    reqs = load_requirements()
    assert "NOT SCORED" in text
    for rid, entry in reqs.not_scored.items():
        assert rid in text
        assert entry.category in text
    assert '"not_scored"' in render_json(assessment)
    assert assessment.flags["coverage"]["n_catalog"] == 72


# ---------------------------------------------------------------------------
# The importer
# ---------------------------------------------------------------------------

_HEADERS = {
    "Measurement Reqs": ["ID", "Requirement", "Threshold", "Goal", "Basis / source",
                         "Verification"],
    "Instrument & Mission": ["ID", "Requirement", "Value", "Flows from",
                             "Current point design / status", "Verification"],
    "Software (tss-sim)": ["ID", "Requirement", "Value", "Flows from", "Verification"],
}


def _cell(ref: str, value: str, shared: list[str]) -> str:
    if value.startswith("@"):  # a shared string, to exercise that path
        shared.append(value[1:])
        return f'<c r="{ref}" t="s"><v>{len(shared) - 1}</v></c>'
    return f'<c r="{ref}" t="inlineStr"><is><t>{escape(value)}</t></is></c>'


def _write_xlsx(path: Path, sheets: dict[str, list[list[str]]]) -> None:
    shared: list[str] = []
    parts: dict[str, str] = {}
    rels, entries = [], []
    for i, (name, rows) in enumerate(sheets.items(), start=1):
        xml_rows = []
        for r, row in enumerate(rows, start=1):
            cells = "".join(
                _cell(f"{chr(ord('A') + c)}{r}", v, shared) for c, v in enumerate(row) if v
            )
            xml_rows.append(f'<row r="{r}">{cells}</row>')
        parts[f"xl/worksheets/sheet{i}.xml"] = (
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            f"<sheetData>{''.join(xml_rows)}</sheetData></worksheet>"
        )
        rels.append(f'<Relationship Id="rId{i}" Target="worksheets/sheet{i}.xml" '
                    'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"/>')
        entries.append(f'<sheet name="{escape(name)}" sheetId="{i}" r:id="rId{i}"/>')
    parts["xl/workbook.xml"] = (
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"<sheets>{''.join(entries)}</sheets></workbook>"
    )
    parts["xl/_rels/workbook.xml.rels"] = (
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f"{''.join(rels)}</Relationships>"
    )
    sst = "".join(f"<si><t>{escape(s)}</t></si>" for s in shared)
    parts["xl/sharedStrings.xml"] = (
        '<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f"{sst}</sst>"
    )
    with zipfile.ZipFile(path, "w") as zf:
        for name, body in parts.items():
            zf.writestr(name, body)


def _tiny_workbook(path: Path) -> None:
    _write_xlsx(path, {
        "SRM": [["Science Traceability Matrix"], ["Rev Q draft (2026-10-17)."]],
        "Measurement Reqs": [
            ["Measurement Requirements"], [], _HEADERS["Measurement Reqs"],
            ["SCI-01", "@Measure  the\nprofile.", "10%", "5%", "A source.", "S, R. In L5."],
            ["", "--- a section heading, not a requirement ---"],
        ],
        "Instrument & Mission": [
            _HEADERS["Instrument & Mission"],
            ["INS-01a", "Lowest frequency", "0.33 MHz", "SCI-03", "status", "A, R, S. Study."],
            ["MIS-01", "Orbit", "500-900 km", "SCI-02", "", "A. Mission level."],
        ],
        "Software (tss-sim)": [
            _HEADERS["Software (tss-sim)"],
            ["SW-01", "Score against truth", "rms", "SCI-05", "T, D. Unit tests."],
        ],
    })


def test_import_srm_reads_every_requirement_sheet(tmp_path):
    srm = _import_srm()
    book = tmp_path / "TSS_SRM_RevQ.xlsx"
    _tiny_workbook(book)
    cat = srm.build_catalogue(book)
    assert cat["revision"] == "Q"
    assert list(cat["requirements"]) == ["SCI-01", "INS-01a", "MIS-01", "SW-01"]
    sci = cat["requirements"]["SCI-01"]
    assert sci["requirement"] == "Measure the profile."      # whitespace collapsed
    assert (sci["threshold"], sci["goal"]) == ("10%", "5%")
    assert sci["verification_methods"] == ["S", "R"]
    assert cat["requirements"]["INS-01a"]["value"] == "0.33 MHz"
    assert cat["requirements"]["SW-01"]["verification_methods"] == ["T", "D"]
    assert cat["counts"] == {"SCI": 1, "INS": 1, "MIS": 1, "SW": 1}
    assert len(cat["source_sha256"]) == 64


def test_import_srm_refuses_a_sheet_whose_columns_moved(tmp_path):
    srm = _import_srm()
    book = tmp_path / "moved.xlsx"
    headers = dict(_HEADERS)
    _write_xlsx(book, {
        "Measurement Reqs": [["ID", "Requirement", "Goal", "Threshold", "Basis / source",
                              "Verification"]],
        "Instrument & Mission": [headers["Instrument & Mission"]],
        "Software (tss-sim)": [headers["Software (tss-sim)"]],
    })
    with pytest.raises(ValueError, match="header"):
        srm.build_catalogue(book)


def test_import_srm_check_mode_detects_drift(tmp_path):
    srm = _import_srm()
    book = tmp_path / "TSS_SRM_RevQ.xlsx"
    _tiny_workbook(book)
    out = tmp_path / "catalog.yaml"
    assert srm.main([str(book), "--out", str(out), "--check"]) == 1   # nothing on disk
    assert srm.main([str(book), "--out", str(out)]) == 0
    assert srm.main([str(book), "--out", str(out), "--check"]) == 0
    out.write_text(out.read_text(encoding="utf-8").replace("Lowest frequency", "Edited"),
                   encoding="utf-8")
    assert srm.main([str(book), "--out", str(out), "--check"]) == 1


@pytest.mark.skipif(not REAL_WORKBOOK.exists(), reason="the SRM workbook is not on this machine")
def test_the_catalogue_matches_the_workbook_on_this_machine():
    """Fails when TSS_SRM_RevG.xlsx has changed and the catalogue was not
    regenerated: rerun scripts/import_srm.py and review the diff."""
    srm = _import_srm()
    expected = srm.render(srm.build_catalogue(REAL_WORKBOOK))
    on_disk = (REPO / "config" / "requirements" / "srm_catalog.yaml").read_text(encoding="utf-8")
    assert on_disk == expected
