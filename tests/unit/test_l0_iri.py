"""The l0/iri tier: IRI-2020 with the implementation a configuration option.

The default and the dispatch are checked without MATLAB or PyIRI (the PHaRLAP
route is built with ``probe=False``, which constructs it without starting
MATLAB).  The PyIRI route is checked against the ``iri2020`` tier in
tests/unit/test_l0_pyiri.py.  The ``pharlap``-marked test at the end runs the
baseline campaign's L0 through MATLAB and PHaRLAP and checks the solar level
IRI itself reports; it is deselected by default and run with ``-m pharlap``.
"""

from __future__ import annotations

import datetime as dt
import inspect
import json

import numpy as np
import pytest

from tss.config import load_campaign
from tss.l0_ionosphere.backends import (
    IRI_DEFAULT_MODEL,
    IRI_MODELS,
    IRIBackend,
    MatlabIRIBackend,
)
from tss.l0_ionosphere.base import iri_f107_from_r12, iri_r12_from_f107
from tss.pipeline import build_tiers
from tss.registry import build, tier_class

ALT = np.arange(1.0e5, 2.0e6 + 1.0, 1.0e4)


def _pharlap(**kw):
    return build("l0", "iri", **{"f107": 120.0, "probe": False, **kw})


def test_the_tier_is_registered_and_defaults_to_pharlap():
    assert tier_class("l0", "iri") is IRIBackend
    assert IRI_DEFAULT_MODEL == "pharlap"
    assert IRI_MODELS == {"pharlap": "matlab_iri", "pyiri": "iri2020"}
    default = inspect.signature(IRIBackend).parameters["model"].default
    assert default == "pharlap"


def test_without_a_model_the_tier_is_pharlap():
    tier = _pharlap()
    assert tier.model == "pharlap"
    assert isinstance(tier.impl, MatlabIRIBackend)
    assert "model=pharlap" in tier.label and "PHaRLAP" in tier.provenance


def test_the_model_name_is_case_and_space_tolerant():
    assert _pharlap(model=" PHaRLAP ").model == "pharlap"


@pytest.mark.parametrize("bad", ["matlab_iri", "iri2016", "nequick", ""])
def test_an_unknown_model_is_refused(bad):
    with pytest.raises(ValueError, match="model must be one of"):
        _pharlap(model=bad)


def test_a_solar_index_is_required_and_only_one():
    with pytest.raises(ValueError, match="one solar index"):
        build("l0", "iri", probe=False)
    with pytest.raises(ValueError, match="one solar index"):
        build("l0", "iri", f107=120.0, r12=70.0, probe=False)


def test_iri_relation_round_trips_and_matches_the_mex_constant():
    r12 = np.array([1.0, 10.0, 71.1, 150.0, 200.0])
    assert np.allclose(iri_r12_from_f107(iri_f107_from_r12(r12)), r12)
    assert iri_f107_from_r12(0.0) == pytest.approx(63.75)
    assert float(iri_r12_from_f107(120.0)) == pytest.approx(71.09, abs=0.01)


def test_f107_reaches_pharlap_as_the_r12_iri_maps_back_to_it():
    tier = _pharlap()
    assert tier.impl.r12 == pytest.approx(float(iri_r12_from_f107(120.0)))
    assert iri_f107_from_r12(tier.impl.r12) == pytest.approx(120.0)
    assert tier.f107 == 120.0 and tier.solar_index_stated == "f107"


def test_r12_can_be_stated_instead():
    tier = _pharlap(f107=None, r12=50.0)
    assert tier.impl.r12 == 50.0
    assert tier.f107 == pytest.approx(float(iri_f107_from_r12(50.0)))


@pytest.mark.parametrize("f107", [60.0, 260.0])
def test_a_solar_level_outside_the_pharlap_range_is_refused(f107):
    with pytest.raises(ValueError, match="R12 in 1..200"):
        _pharlap(f107=f107)


def test_options_reach_the_chosen_implementation():
    tier = _pharlap(topside="IRI-2001", timeout=30.0, pharlap_root="X:/pharlap")
    assert tier.impl.topside == "IRI-2001"
    assert tier.impl.timeout == 30.0 and tier.impl.pharlap_root == "X:/pharlap"
    # Attributes this class lacks are read from the implementation.
    assert tier.pharlap_root == "X:/pharlap" and tier.extra is None


def test_an_option_the_implementation_does_not_take_is_an_error():
    with pytest.raises(TypeError, match="MatlabIRIBackend"):
        _pharlap(foF2_coeff="URSI")


def test_describe_records_the_choice_in_json():
    record = dict(_pharlap().describe())
    json.dumps(record)
    assert record["name"] == "iri" and record["model"] == "pharlap"
    assert record["implementation"] == "l0/matlab_iri"
    assert record["implementation_settings"]["r12"] == pytest.approx(71.09, abs=0.01)


def test_profiles_and_extras_are_the_implementations(monkeypatch):
    tier = _pharlap()
    ne = np.full((1, ALT.size), 1.0e11)
    monkeypatch.setattr(tier.impl, "profiles", lambda *a, **k: ne)
    monkeypatch.setattr(tier.impl, "_extra_for_point", lambda i: {"nmf2": 2.0e11, "hmf2_km": 310.0})
    state = tier.state_at(10.0, 20.0, "2024-03-01T12:00:00Z", ALT)
    assert np.array_equal(state.ne, ne[0])
    assert state.extra["nmf2"] == 2.0e11 and state.extra["hmf2_km"] == 310.0
    assert state.source == tier.label


def test_the_baseline_campaign_uses_pharlap_and_states_its_solar_level():
    l0 = load_campaign("baseline").config["tiers"]["l0"]
    assert l0["name"] == "iri" and l0["model"] == "pharlap" and l0["f107"] == 120.0


def test_a_campaign_that_omits_the_model_gets_pharlap():
    tiers = build_tiers({"l0": {"name": "iri", "f107": 120.0, "probe": False}})
    assert tiers["l0"].model == "pharlap"


@pytest.mark.pharlap
def test_pharlap_route_runs_at_the_stated_solar_level():
    """Through MATLAB and PHaRLAP: IRI's own reported F10.7 is the one stated."""
    tier = build("l0", "iri", f107=120.0)
    direct = build("l0", "matlab_iri", r12=float(iri_r12_from_f107(120.0)))
    when = [dt.datetime(2024, 3, 1, 12)]
    ne = tier.profiles([42.6], [-71.5], when, ALT)
    assert np.allclose(ne, direct.profiles([42.6], [-71.5], when, ALT))
    extra = tier._extra_for_point(0)
    assert extra["f107"] == pytest.approx(120.0, abs=0.05)
    assert 1.0e11 < extra["nmf2"] < 5.0e12 and 200.0 < extra["hmf2_km"] < 450.0
