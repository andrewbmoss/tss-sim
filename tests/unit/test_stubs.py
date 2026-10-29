"""The stubs raise, and they say which work package will fill them in.

A stub that raises a bare :class:`NotImplementedError` tells a reader that
something is missing.  A stub that names its work package tells them whether it
is missing because nobody has got to it, because it is a goal-level objective,
or because it is blocked on a decision that has not been made.  That difference
is the reason these tests exist: they hold the message, not only the exception.

Both stubs are importable and, where they are registered, buildable.  That is so
a campaign naming one can be validated before it is run, which is what
``tss/registry.py`` promises in its module docstring.

Covered here:

``tss.l0_ionosphere.perturbations``   WP10, goal level.  Equatorial plasma
                                      bubble depletions.
``tss.l2_instrument.leakage``         WP3, critical path, C4.
                                      Transmit-to-receive leakage and its
                                      de-chirped residual.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from tss.interfaces import IonosphereState, SweepSchedule  # noqa: E402
from tss.l0_ionosphere import perturbations  # noqa: E402
from tss.l2_instrument import leakage  # noqa: E402


def _state() -> IonosphereState:
    alt = np.arange(300.0e3, 900.0e3, 10.0e3)
    return IonosphereState(alt_m=alt, ne=np.full(alt.shape, 1.0e11))


def _sweep() -> SweepSchedule:
    return SweepSchedule(
        f_start_hz=3.3e5, f_stop_hz=2.0e7, k_f_hz_per_s=5.6e6, sweep_period_s=3.3
    )


# ---------------------------------------------------------------------------
# WP10: plasma bubble depletions
# ---------------------------------------------------------------------------
def test_perturbations_bubble_parameters_names_wp10():
    with pytest.raises(NotImplementedError) as excinfo:
        perturbations.BubbleParameters()
    assert "WP10" in str(excinfo.value)


def test_perturbations_depletion_factor_names_wp10():
    rng = np.random.default_rng(0)
    with pytest.raises(NotImplementedError) as excinfo:
        perturbations.bubble_depletion_factor(
            0.0, 0.0, np.array([300.0e3, 400.0e3]), "2026-01-01T00:00:00", None, rng
        )
    assert "WP10" in str(excinfo.value)


def test_perturbations_perturb_state_names_wp10():
    rng = np.random.default_rng(0)
    with pytest.raises(NotImplementedError) as excinfo:
        perturbations.perturb_state(_state(), None, rng)
    assert "WP10" in str(excinfo.value)


def test_perturbations_keeps_the_measured_statistics_it_is_parameterised_on():
    """The defaults are measurements, so they belong to the stub and not to the
    implementation that will eventually consume them."""
    assert perturbations.DEFAULT_ZONAL_ENVELOPE_M == pytest.approx(300.0e3)
    assert perturbations.DEFAULT_FIELD_ALIGNED_SCALE_M == pytest.approx(550.0e3)
    assert tuple(perturbations.DEFAULT_DEPLETION_DECADES) == (2.0, 3.0)
    assert tuple(perturbations.DEFAULT_EASTWARD_DRIFT_MPS) == (75.0, 250.0)


# ---------------------------------------------------------------------------
# WP3: transmit-to-receive leakage.  The model (Phase 1), the receive chain
# (Phase 3) and the residual and its score (Phase 4) are built and tested in
# test_leakage.py, test_receive_chain.py and test_leakage_scoring.py.
# ---------------------------------------------------------------------------
def test_leakage_keeps_the_requirement_thresholds_it_will_be_scored_against():
    assert leakage.BASELINE_ISOLATION_DB == pytest.approx(40.0)
    assert leakage.BASELINE_COUPLING_DELAY_S == pytest.approx(3.0e-9)


# ---------------------------------------------------------------------------
# A run made without them has to say so
# ---------------------------------------------------------------------------
def test_the_l2_baseline_tier_records_that_leakage_is_not_modelled():
    """The gap is recorded on every record the tier produces, not only in a
    document, because a reader of a result does not have the document."""
    from tss.l2_instrument import BaselineInstrument

    record = BaselineInstrument().describe()
    assert record["leakage_modelled"] is False
