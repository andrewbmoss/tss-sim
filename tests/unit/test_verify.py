"""Unit tests for the verification harness.

The harness is tested, not the physics.  What matters is that it reports
honestly: a check whose subject is not written must come back ``not_run`` with
the work package named, a check that could run but for a missing MATLAB must
come back ``skipped``, and neither may ever come back ``pass``.  A harness that
turned an unwritten algorithm into a silent pass would make every verification
report worthless, and it would do so invisibly.

``run_verification`` is also tested for the property that its exit status alone
is not evidence: a run in which nothing ran and nothing failed returns zero, so
the summary line has to say how many checks actually ran.
"""

from __future__ import annotations

import io
import json

import pytest

from tss.verify import (
    LEVELS,
    VERDICTS,
    CheckResult,
    NotRun,
    Skipped,
    registered,
    render_markdown,
    run_all,
    run_level,
    run_verification,
    summarise,
    write_report,
)


def test_all_four_levels_have_registered_checks():
    """A level with no checks is a level that cannot report anything."""
    names = registered()
    assert set(names) == set(LEVELS)
    for level in LEVELS:
        assert names[level], f"{level} has no registered checks"


def test_v1_checks_actually_run_and_pass():
    """The closed-form cases are real checks, not placeholders.

    If none of V1 passed, the harness would be reporting an empty suite as a
    clean one.
    """
    results = run_level("V1")
    verdicts = [r.verdict for r in results]
    assert verdicts.count("pass") >= 8
    assert "fail" not in verdicts


def test_no_check_reports_an_unknown_verdict():
    for records in run_all().values():
        for record in records:
            assert record.verdict in VERDICTS


def test_unwritten_physics_reports_not_run_and_never_pass():
    """The property the whole harness exists for.

    Checks registered against WP5 and WP2 have nothing to test yet.  They must
    say so, name the work package, and must not pass.
    """
    by_name = {r.name: r for r in run_level("V1") + run_level("V2")}
    for name in ("c3_estimator_recovers_truth_velocity", "leakage_cancellation_meets_ins24"):
        record = by_name[name]
        assert record.verdict == "not_run"
        assert "WP5" in record.reason
    assert by_name["analytic_vs_raytrace_group_delay"].verdict == "not_run"
    assert "WP2" in by_name["analytic_vs_raytrace_group_delay"].reason


def test_every_not_run_and_skipped_check_carries_a_reason():
    """A verdict without a reason gives the reader nothing to act on."""
    for records in run_all().values():
        for record in records:
            if record.verdict in ("not_run", "skipped"):
                assert record.reason, f"{record.level}/{record.name} has no reason"


def test_missing_external_dependency_is_skipped_not_not_run():
    """The two are different and the difference is the point.

    A skipped check might pass on another machine.  A not-run check will not
    pass anywhere until its work package lands.  Collapsing them would let the
    codebase's real coverage hide behind an absent MATLAB install.
    """
    by_name = {r.name: r for r in run_level("V2")}
    assert by_name["python_matlab_transcription"].verdict == "skipped"
    assert "MATLAB" in by_name["python_matlab_transcription"].reason


def test_v3_registers_all_five_replications():
    """Barry, Poole, Zawdie, Franklin and MacLean, Chartier."""
    names = registered("V3")["V3"]
    for stem in ("barry", "poole", "zawdie", "franklin", "chartier"):
        assert any(stem in name for name in names), f"no V3 check mentions {stem}"


def test_v4_records_the_permanent_limitation_as_a_check():
    """The archive cannot validate the FMCW half, and the report says so.

    Registering it as a check rather than leaving it in the document means the
    limitation appears in every verification report, where a reader of the
    report will actually see it.
    """
    by_name = {r.name: r for r in run_level("V4")}
    record = by_name["fmcw_half_not_reachable_from_archive"]
    assert record.verdict == "not_run"
    assert "pulse sounders" in record.reason


def test_a_check_that_raises_unexpectedly_is_a_failure_not_a_skip():
    """A crashing check has found something and must not be excused."""
    from tss.verify.levels import _Registration, _run_one

    def _boom():
        raise ZeroDivisionError("division by zero")

    record = _run_one(_Registration(level="V1", name="boom", fn=_boom))
    assert record.verdict == "fail"
    assert "ZeroDivisionError" in record.detail


def test_not_run_and_skipped_are_distinguished_by_exception_type():
    from tss.verify.levels import _Registration, _run_one

    not_run = _run_one(
        _Registration(level="V1", name="a", fn=lambda: (_ for _ in ()).throw(NotRun("WP9: no bench")))
    )
    skipped = _run_one(
        _Registration(level="V1", name="b", fn=lambda: (_ for _ in ()).throw(Skipped("no MATLAB")))
    )
    assert not_run.verdict == "not_run"
    assert skipped.verdict == "skipped"


def test_a_stubs_own_message_becomes_the_reason():
    """A check written against a stub reports the work package automatically.

    It starts passing or failing honestly the moment the algorithm lands, with
    no change to the check.
    """
    from tss.verify.levels import _Registration, _run_one

    def _calls_a_stub():
        from tss.l3_signal.novel.leakage_cancel import residual_floor_db

        residual_floor_db(None)

    record = _run_one(_Registration(level="V1", name="stubby", fn=_calls_a_stub))
    assert record.verdict == "not_run"
    assert "WP5" in record.reason


def test_run_verification_returns_zero_when_nothing_failed():
    """And prints how many checks actually ran, because the status alone is not evidence."""
    buffer = io.StringIO()
    status = run_verification("all", write=False, stream=buffer)
    output = buffer.getvalue()
    assert status == 0
    assert "actually ran" in output
    assert "not run" in output


def test_run_verification_rejects_an_unknown_level():
    with pytest.raises(ValueError, match="unknown verification level"):
        run_verification("V9", stream=io.StringIO())


def test_summary_keeps_zero_counts():
    counts = summarise({level: [] for level in LEVELS})
    assert set(counts) == set(VERDICTS)
    assert all(value == 0 for value in counts.values())


def test_report_says_when_it_is_weak():
    """``docs/verification.md`` requires a report with many skips to say so.

    Built from a synthetic result set in which almost nothing ran, so the test
    does not depend on how much of the real suite happens to be written today.
    """
    results = {
        "V1": [CheckResult("V1", "ok", "pass")],
        "V2": [CheckResult("V2", "a", "not_run", reason="WP2")],
        "V3": [CheckResult("V3", "b", "not_run", reason="WP8")],
        "V4": [CheckResult("V4", "c", "skipped", reason="no cache")],
    }
    markdown = render_markdown(results)
    assert "weak report" in markdown
    assert "Only 1 of 4 checks actually ran" in markdown


def test_report_writes_json_and_markdown():
    """Both renderings land under the given directory and the JSON parses."""
    import tempfile
    from pathlib import Path

    results = run_all()
    with tempfile.TemporaryDirectory() as directory:
        json_path, md_path = write_report(results, directory)
        assert json_path.exists() and md_path.exists()
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        assert set(payload["levels"]) == set(LEVELS)
        assert payload["total_checks"] == sum(len(v) for v in results.values())
        assert "Verification report" in Path(md_path).read_text(encoding="utf-8")


def test_report_lists_every_registered_check():
    """Nothing is omitted for brevity.  The report is the full registry."""
    results = run_all()
    markdown = render_markdown(results)
    for _level, names in registered().items():
        for name in names:
            assert f"`{name}`" in markdown
