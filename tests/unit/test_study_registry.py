"""The study registry tells the truth about what runs (code review F6, 2026-11-08).

``transcription`` was registered ``implemented`` but could never produce a
result: it raised ``StudyInputMissing`` whether MATLAB was present or not, and
sent the reader to MATLAB files that are not in the repository.  A study that
cannot run is a ``stub`` naming the work package that owns it, and running it
returns exit status 2 with that work package in the message.
"""

from __future__ import annotations

import io

import pytest

from tss.studies import STUDIES, run_study


def test_the_transcription_study_is_a_stub_owned_by_wp8():
    study = STUDIES["transcription"]
    assert study.status == "stub"
    assert study.blocked_on == "WP8"


def test_every_stub_names_the_work_package_it_waits_on():
    stubs = [s for s in STUDIES.values() if s.status == "stub"]
    assert stubs
    for study in stubs:
        assert study.blocked_on, study.name  # its own or another (trade_space waits on WP4)


@pytest.mark.parametrize("name", sorted(n for n, s in STUDIES.items() if s.status == "stub"))
def test_running_a_stub_returns_2_and_names_its_work_package(name, tmp_path):
    stream = io.StringIO()
    status = run_study(name, out=tmp_path / name, stream=stream)
    assert status == 2, stream.getvalue()
    assert STUDIES[name].blocked_on in stream.getvalue()
    assert not (tmp_path / name).exists(), "a stub must not leave a results directory behind"
