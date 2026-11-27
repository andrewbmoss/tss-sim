"""WP5 Phase 3 part 3: the diffuse-halo experiment script (``scripts/wp5_halo.py``)."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def script():
    path = Path(__file__).resolve().parents[2] / "scripts" / "wp5_halo.py"
    spec = importlib.util.spec_from_file_location("wp5_halo", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_cases_are_the_six_agreed_and_no_halo(script):
    """Andrew, 2026-11-27: 2 facet lengths at K = 10, 3, 1; 1 and 4 at K = 3;
    one K = 0 check."""
    assert script.CASES == {"no halo": None, "r2 K10": (2.0, 10.0), "r2 K3": (2.0, 3.0), "r2 K1": (2.0, 1.0),
                            "r1 K3": (1.0, 3.0), "r4 K3": (4.0, 3.0), "r2 K0": (2.0, 0.0)}


def test_one_realisation_gives_a_row_per_level_and_estimator(script):
    rows = script.trial(("r2 K10", 3), estimators=("pinned", "joint"))
    assert len(rows) == len(script.LEVELS) * 2
    assert {(r["level_db"], r["estimator"]) for r in rows} == {
        (lv, e) for lv in script.LEVELS for e in ("pinned", "joint")}
    clean = next(r for r in rows if r["level_db"] is None and r["estimator"] == "pinned")
    assert clean["detected"] and clean["formed"] and abs(clean["err_mps"]) < 25.0
    for r in rows:
        assert math.isfinite(r["err_mps"]) == r["formed"]


def test_the_summary_counts_detections_and_formed_drifts(script):
    rows = [{"case": "r2 K3", "seed": s, "level_db": 10.0, "estimator": "pinned", "detected": True,
             "formed": s < 3, "err_mps": float(s) if s < 3 else math.nan} for s in range(4)]
    out = script.summarise(rows)[("r2 K3", 10.0, "pinned")]
    assert out["n"] == 4 and out["detected"] == 4 and out["formed"] == 3
    assert out["rms"] == pytest.approx(math.sqrt(5 / 3)) and out["bias"] == pytest.approx(1.0)
