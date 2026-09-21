"""Unit tests for the L0 ionosphere layer.

Every test here runs on a bare machine.  Nothing imports iri2016, PyIRI, ppigrf
or MATLAB, nothing opens a socket and nothing shells out.  That is the property
the ``chapman`` tier exists to have: it is named as the default in
``config/default.yaml``, so if it needed an optional dependency then the package
default would be unrunnable on a fresh checkout, which is the sort of thing that
is discovered at the worst moment.

The tests check four things.  That the default tier returns a well-formed
:class:`~tss.interfaces.IonosphereState` with every slot L1 needs populated.
That the plasma frequency relation round-trips, because everything about the
band edges is that relation.  That the registry holds the tier names a campaign
configuration is allowed to name.  And that ``state_at`` does not modify its
inputs, which is contract rule 1 and is the one rule that fails silently.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from tss.constants import density_from_plasma_frequency, plasma_frequency  # noqa: E402
from tss.interfaces import IonosphereState, L0Ionosphere  # noqa: E402
from tss.l0_ionosphere import (  # noqa: E402
    AnalyticBackend,
    ChapmanIonosphere,
    chapman_density,
    chapman_scale_height_m,
    dipole_field_ned_t,
    state_from_profile,
)
from tss.registry import UnknownTier, available, build  # noqa: E402

EPOCH = "2015-06-21T12:00:00"


def _grid():
    return np.arange(90.0e3, 2000.0e3 + 1.0, 5.0e3)


# ---------------------------------------------------------------------------
# The default tier on a machine with nothing installed
# ---------------------------------------------------------------------------
def test_chapman_returns_a_well_formed_state():
    tier = ChapmanIonosphere()
    alt = _grid()
    state = tier.state_at(41.7, -111.8, EPOCH, alt)

    assert isinstance(state, IonosphereState)
    assert state.alt_m.shape == alt.shape
    assert state.ne.shape == alt.shape
    assert np.all(np.isfinite(state.ne))
    assert np.all(state.ne > 0.0)

    # L1 absorption cannot run without nu, so a tier that leaves it None is not
    # usable for a scored run.  The architecture document says so; this holds it.
    assert state.nu is not None and state.nu.shape == alt.shape
    assert np.all(np.isfinite(state.nu)) and np.all(state.nu > 0.0)
    assert state.te is not None and np.all(state.te > 0.0)

    # The magnetoionic splitting at L1 and the X-mode cutoff of the band-edge
    # study both need the field, and it has to come out as (3, N) north-east-down.
    assert state.b_field_t is not None
    assert state.b_field_t.shape == (3, alt.size)
    assert np.all(np.isfinite(state.b_field_t))

    assert state.lat_deg == pytest.approx(41.7)
    assert state.lon_deg == pytest.approx(-111.8)
    assert state.epoch == EPOCH
    assert state.source


def test_chapman_peaks_at_hmf2_and_decays_upward():
    tier = ChapmanIonosphere(nmf2=5.0e11, hmf2_m=300.0e3, plasma_scale_height_m=100.0e3)
    alt = _grid()
    state = tier.state_at(0.0, 0.0, EPOCH, alt)
    peak_index = int(np.argmax(state.ne))
    assert alt[peak_index] == pytest.approx(300.0e3, abs=5.0e3)
    assert math.isclose(float(np.max(state.ne)), 5.0e11, rel_tol=1e-6)
    above = state.ne[peak_index:]
    assert np.all(np.diff(above) < 0.0)


def test_chapman_uses_the_plasma_scale_height_not_the_chapman_one():
    """Far above the peak the density decays with ``Hp``, not with ``Hp/2``.

    This is the single easiest thing in the module to get wrong, and getting it
    wrong is a factor of two in topside density that looks entirely plausible.
    """
    hp = 100.0e3
    assert chapman_scale_height_m(hp) == pytest.approx(hp / 2.0)
    hmf2 = 300.0e3
    far = np.array([hmf2 + 6.0 * hp, hmf2 + 7.0 * hp])
    ne = chapman_density(far, 5.0e11, hmf2, hp)
    decay_scale = (far[1] - far[0]) / math.log(ne[0] / ne[1])
    assert decay_scale == pytest.approx(hp, rel=1e-3)


def test_chapman_needs_no_optional_dependency():
    """The default tier must run with numpy alone.

    ppigrf is optional, so the geomagnetic field has to fall back.  Asking for
    the dipole explicitly is the bare-machine case, and ``describe()`` has to say
    which route answered, because a result that used the dipole and a result that
    used IGRF are not the same result.
    """
    tier = ChapmanIonosphere(geomagnetic_model="dipole")
    state = tier.state_at(0.0, 0.0, EPOCH, _grid())
    assert state.b_field_t is not None
    record = tier.describe()
    assert record["geomagnetic_provenance"]["model"] == "centred dipole"
    assert record["name"] == "chapman"
    assert record["layer"] == "l0"
    assert record["plasma_scale_height_m"] == pytest.approx(100.0e3)


def test_chapman_satisfies_the_l0_protocol():
    assert isinstance(ChapmanIonosphere(), L0Ionosphere)


# ---------------------------------------------------------------------------
# The plasma relation
# ---------------------------------------------------------------------------
def test_plasma_frequency_round_trips():
    for ne in (1e9, 1e10, 1e11, 1e12):
        assert math.isclose(
            float(density_from_plasma_frequency(plasma_frequency(ne))), ne, rel_tol=1e-12
        )


def test_plasma_frequency_round_trips_through_a_state():
    """foF2 taken from the state's own peak must match the reported fof2_hz."""
    tier = ChapmanIonosphere(nmf2=1.0e12)
    state = tier.state_at(0.0, 0.0, EPOCH, _grid())
    fof2 = float(plasma_frequency(np.max(state.ne)))
    assert math.isclose(fof2, float(state.extra["fof2_hz"]), rel_tol=1e-9)
    assert math.isclose(
        float(density_from_plasma_frequency(fof2)), float(np.max(state.ne)), rel_tol=1e-12
    )


# ---------------------------------------------------------------------------
# The registry
# ---------------------------------------------------------------------------
def test_registry_holds_the_expected_l0_tier_names():
    names = available("l0")["l0"]
    for expected in ("analytic", "chapman", "iri2020", "iri_wrapper", "matlab_iri"):
        assert expected in names, f"{expected} is not registered; got {names}"


def test_registry_builds_the_default_tier_by_name():
    tier = build("l0", "chapman", plasma_scale_height_m=150.0e3)
    assert tier.name == "chapman"
    assert tier.plasma_scale_height_m == pytest.approx(150.0e3)


def test_optional_tiers_are_registered_even_where_their_model_is_absent():
    """A campaign naming a tier must validate before it runs.

    The registry holds the class whether or not the model behind it is installed.
    Building it is what fails, and it fails with a message about the model rather
    than with a KeyError about the name.
    """
    names = available("l0")["l0"]
    assert "matlab_iri" in names
    # The failure must be about the missing model, not about the tier name, so
    # anything except UnknownTier is the correct outcome here.
    try:
        build("l0", "iri_wrapper", module="a_module_that_does_not_exist")
    except UnknownTier:  # pragma: no cover - would be a registry bug
        raise AssertionError("registration should not depend on the model being installed") from None
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Contract rule 1: inputs are never mutated
# ---------------------------------------------------------------------------
def test_state_at_does_not_mutate_its_inputs():
    tier = ChapmanIonosphere()
    alt = _grid()
    before = alt.copy()
    state = tier.state_at(10.0, 20.0, EPOCH, alt)
    assert np.array_equal(alt, before)

    # And the returned state must not alias the caller's array, or a later write
    # by the caller would silently rewrite the truth a run was scored against.
    assert state.alt_m is not alt
    state.alt_m[0] += 1.0
    assert np.array_equal(alt, before)


def test_analytic_backend_does_not_mutate_its_inputs():
    tier = AnalyticBackend()
    alt = _grid()
    before = alt.copy()
    state = tier.state_at(-20.0, 130.0, EPOCH, alt)
    assert np.array_equal(alt, before)
    assert state.alt_m is not alt
    assert state.ne.shape == alt.shape
    assert np.all(state.ne > 0.0)


def test_state_from_profile_does_not_mutate_the_profile_or_the_grid():
    class _Profile:
        def __init__(self):
            self.calls = 0

        def ne(self, h):
            self.calls += 1
            return 5.0e11 * np.exp(-(np.asarray(h, dtype=float) - 300.0e3) / 200.0e3)

        def te(self, h):
            return np.full(np.shape(np.asarray(h, dtype=float)), 1500.0)

    profile = _Profile()
    alt = _grid()
    before = alt.copy()
    state = state_from_profile(profile, alt, lat_deg=1.0, lon_deg=2.0, epoch=EPOCH)
    assert np.array_equal(alt, before)
    assert state.alt_m is not alt
    assert state.nu is not None and np.all(state.nu > 0.0)
    assert "hmf2_km" in state.extra


# ---------------------------------------------------------------------------
# The geomagnetic fallback
# ---------------------------------------------------------------------------
def test_dipole_field_is_stronger_at_the_pole_than_at_the_equator():
    equator = dipole_field_ned_t(0.0, 0.0, np.array([300.0e3]))
    pole = dipole_field_ned_t(80.65, -72.68, np.array([300.0e3]))
    assert np.linalg.norm(pole) > 1.9 * np.linalg.norm(equator)


def test_dipole_field_falls_off_as_the_inverse_cube():
    from tss.constants import R_EARTH_M

    alt = np.array([0.0, R_EARTH_M])
    b = dipole_field_ned_t(0.0, 0.0, alt)
    magnitude = np.linalg.norm(b, axis=0)
    assert magnitude[0] / magnitude[1] == pytest.approx(8.0, rel=1e-9)


# ---------------------------------------------------------------------------
# Grid validation
# ---------------------------------------------------------------------------
def test_a_bad_altitude_grid_raises_rather_than_returning_nonsense():
    """A malformed grid is a configuration error, not a physical impossibility.

    Contract rule 3 says physical impossibility is a result and not an error.  A
    descending altitude grid is neither; it is a bug in the caller, and raising
    is right.
    """
    tier = ChapmanIonosphere()
    with pytest.raises(ValueError):
        tier.state_at(0.0, 0.0, EPOCH, np.array([500.0e3, 400.0e3, 300.0e3]))
    with pytest.raises(ValueError):
        tier.state_at(0.0, 0.0, EPOCH, np.array([500.0e3]))
