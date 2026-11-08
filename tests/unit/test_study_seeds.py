"""Code review F9: every seed a study draws from is named and reaches its manifest.

``Manifest.seed`` is the campaign seed from ``config["seed"]``.  A study that
seeds its own generators from a literal buried in a function leaves a manifest
that cannot reproduce it, so such seeds are module constants and the study
records them in ``Manifest.seeds``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

import tss.studies as studies_pkg
from tss.provenance import Manifest
from tss.studies import leakage_budget, leakage_receive_chain, wp4_dispersion

_LITERAL_SEED = re.compile(r"default_rng\(\s*\d+\s*\)")


def test_the_manifest_writes_the_seeds_it_was_given(tmp_path):
    m = Manifest(study="x", seeds={"leakage": 0, "receive": [1, 2]})
    data = json.loads(m.write(tmp_path).read_text())
    assert data["seeds"] == {"leakage": 0, "receive": [1, 2]}


def test_a_manifest_starts_with_no_seeds():
    assert Manifest(study="x").seeds == {}


def test_no_study_seeds_a_generator_from_a_buried_literal():
    """``default_rng(<int>)`` is allowed only as the ``rng is None`` fallback of a
    helper whose study passes its own generator."""
    offenders = []
    for path in sorted(Path(studies_pkg.__file__).parent.glob("*.py")):
        for no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if _LITERAL_SEED.search(line) and "if rng is None" not in line:
                offenders.append(f"{path.name}:{no}: {line.strip()}")
    assert not offenders, "\n".join(offenders)


class _Stop(Exception):
    pass


def _stop(*args, **kwargs):
    raise _Stop


@pytest.mark.parametrize(
    ("module", "expected"),
    [
        (leakage_receive_chain, {"leakage": leakage_receive_chain.LEAKAGE_SEED}),
        (leakage_budget, {"receive": list(leakage_budget.SEEDS),
                          "start_up": leakage_budget.START_UP_SEED,
                          "validation_sample": leakage_budget.VALIDATION_SAMPLE_SEED}),
        (wp4_dispersion, {"l1_echoes": wp4_dispersion.L1_ECHO_SEED,
                          "receive": [int(s) for s in wp4_dispersion.DEFAULT_SEEDS]}),
    ],
    ids=lambda v: getattr(v, "__name__", "").rsplit(".", 1)[-1] or None,
)
def test_the_study_records_its_seeds_before_it_computes(module, expected, monkeypatch, tmp_path):
    monkeypatch.setattr(module, "compute", _stop)
    m = Manifest(study=module.__name__)
    with pytest.raises(_Stop):
        module.run(tmp_path, config={}, manifest=m)
    assert m.seeds == expected


def test_the_seeds_keep_their_historical_values():
    """Naming the seeds must not move the published numbers."""
    assert leakage_receive_chain.LEAKAGE_SEED == 0
    assert leakage_budget.VALIDATION_SAMPLE_SEED == 12345
    assert leakage_budget.START_UP_SEED == 0
    assert wp4_dispersion.L1_ECHO_SEED == 0
