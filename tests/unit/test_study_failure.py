"""Code review F13: a study that raises unexpectedly still leaves a manifest.

Before the fix only ``NotImplementedError`` (exit 2) and ``StudyInputMissing``
(exit 3) were caught; anything else escaped :func:`tss.studies.run_study` and left
a run directory with partial outputs and no manifest saying what happened.
"""

from __future__ import annotations

import io
import json

import pytest

from tss.studies import STUDIES, StudyResult, register_study, run_study

_NAME = "_test_raises_unexpectedly"


@pytest.fixture
def failing_study():
    @register_study(_NAME, description="raises after writing one file", work_packages=["WP0"],
                    status="implemented")
    def _run(out_dir, *, config, manifest):
        (out_dir / "partial.csv").write_text("a,b\n1,2\n", encoding="utf-8")
        raise ZeroDivisionError("the denominator was empty")

    yield _NAME
    STUDIES.pop(_NAME, None)


def test_an_unexpected_exception_returns_exit_6(failing_study, tmp_path):
    out = io.StringIO()
    assert run_study(failing_study, out=tmp_path / "run", stream=out) == 6
    text = out.getvalue()
    assert "ZeroDivisionError" in text and "the denominator was empty" in text
    assert "manifest.json" in text


def test_the_failed_run_still_writes_a_manifest_with_the_error(failing_study, tmp_path):
    run_dir = tmp_path / "run"
    run_study(failing_study, out=run_dir, stream=io.StringIO())
    data = json.loads((run_dir / "manifest.json").read_text())
    assert data["error"]["type"] == "ZeroDivisionError"
    assert data["error"]["message"] == "the denominator was empty"
    assert "Traceback" in data["error"]["traceback"]
    assert any(p.endswith("partial.csv") for p in data["outputs"])


def test_a_failed_run_does_not_move_the_latest_pointer(failing_study, tmp_path):
    run_dir = tmp_path / "run"
    run_study(failing_study, out=run_dir, stream=io.StringIO())
    assert not (run_dir / "latest.json").exists()
    assert not (run_dir / "summary.md").exists()


def test_a_successful_run_records_no_error(tmp_path):
    @register_study("_test_succeeds", description="ok", work_packages=["WP0"], status="implemented")
    def _run(out_dir, *, config, manifest):
        return StudyResult(summary="fine")

    try:
        assert run_study("_test_succeeds", out=tmp_path / "ok", stream=io.StringIO()) == 0
        data = json.loads((tmp_path / "ok" / "manifest.json").read_text())
        assert data["error"] is None
    finally:
        STUDIES.pop("_test_succeeds", None)
