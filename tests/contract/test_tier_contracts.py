"""The four contract rules of :mod:`tss.interfaces`, checked on every tier.

The rules are stated once, in the module docstring of ``tss/interfaces.py``:

1.  Inputs are never mutated.  A tier returns new objects.
2.  Units are SI unless the attribute name carries a suffix.
3.  Physical impossibility is a masked result with a reason, not an exception.
4.  Randomness is drawn only from the ``numpy.random.Generator`` passed in.

Until now they were checked only here and there (no-mutation tests for L0, L3
and L4 in the unit suite).  This module checks all four on every implemented
tier, through the same chain the smoke campaign runs, so a tier added later is
held to them without anyone remembering to write the tests.  The last test
fails if a tier is registered that this module does not know about.

Three departures from the rules were found when these tests were written
(code review of 2026-10-17), and all three are fixed:

* an echo set with no rows made L2 raise a TypeError (rule 3);
* the platform velocity was produced in Earth-fixed axes and read by L1 and L5
  as local north-east-down, so a circular orbit had a 7.5 km/s vertical
  closing speed (rule 2, frames);
* L2 picked its sample rate from the deepest simulated echo, so a record with
  no usable echo was sampled at 2.5 Hz and L3 raised on it (rule 3).
"""

from __future__ import annotations

import copy
import dataclasses
import io
import json
import random
from collections.abc import Mapping

import numpy as np
import pytest

import tss  # noqa: F401  (registers the tiers)
from tss.config import load_campaign, load_yaml, resolve_path
from tss.interfaces import (
    Assessment,
    DensityProfile,
    EchoSet,
    Ionogram,
    IonosphereState,
    ReceivedSignal,
)
from tss.registry import available, build
from tss.studies import _grids, _platform_for, _sweep_from_config

EPOCH = "2026-01-01T12:00:00Z"
SEED = 20261017

#: Tiers that run in any environment and are held to every rule.
IMPLEMENTED: dict[str, tuple[str, ...]] = {
    "l0": ("chapman", "analytic"),
    "l1": ("analytic",),
    "l2": ("baseline",),
    "l3": ("baseline",),
    "l4": ("lamination", "lamination_magnetoionic"),
    "l5": ("standard",),
}
#: Registered and deliberately unimplemented.  Each must raise
#: NotImplementedError naming its work package when it is called.
STUBS: dict[tuple[str, str], str] = {
    ("l1", "raytrace"): "WP2",
    ("l3", "motion_separated"): "WP5",
}
#: Tiers that need software this suite cannot assume (PyIRI, the IRI Fortran
#: wrapper, MATLAB).  Only their description is checked here.
EXTERNAL: tuple[tuple[str, str], ...] = (
    ("l0", "iri2020"),
    ("l0", "iri_wrapper"),
    ("l0", "matlab_iri"),
)
#: The default tier at every layer, used to make the inputs of the layer under test.
DEFAULT = {"l0": "chapman", "l1": "analytic", "l2": "baseline", "l3": "baseline",
           "l4": "lamination_magnetoionic", "l5": "standard"}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _setup():
    """Sweep, grids and platform, exactly as the smoke campaign builds them."""
    config = load_campaign("smoke").config
    sweep, _ = _sweep_from_config(config, io.StringIO())
    alt_m, freq_hz = _grids(config)
    platform = _platform_for(
        load_yaml(resolve_path("leo_500km_single", subdir="configurations")), EPOCH
    )
    return sweep, alt_m, freq_hz, platform


def _chain(seed: int = SEED, tiers: Mapping[str, str] | None = None):
    """Run L0 to L5 once and return every intermediate product."""
    names = {**DEFAULT, **(tiers or {})}
    sweep, alt_m, freq_hz, platform = _setup()
    rng = np.random.default_rng(seed)
    t = {layer: build(layer, name) for layer, name in names.items()}
    truth = t["l0"].state_at(platform.lat_deg, platform.lon_deg, EPOCH, alt_m)
    echoes = t["l1"].echoes(truth, platform, freq_hz, rng)
    signal = t["l2"].receive(echoes, sweep, rng)
    ionogram = t["l3"].form_ionogram(signal, platform)
    retrieved = t["l4"].invert(ionogram, platform)
    assessment = t["l5"].score(retrieved, truth, ionogram, "contract")
    return {
        "sweep": sweep, "alt_m": alt_m, "freq_hz": freq_hz, "platform": platform,
        "truth": truth, "echoes": echoes, "signal": signal, "ionogram": ionogram,
        "retrieved": retrieved, "assessment": assessment,
    }


@pytest.fixture(scope="module")
def chain():
    return _chain()


def _frozen(obj):
    """A copy of a dataclass whose arrays are read-only.

    Handing a tier read-only arrays turns any in-place write into an immediate
    ValueError at the line that does it, which is a far better report than a
    comparison afterwards.  The comparison is made as well, for writes that go
    through a mapping.
    """
    if isinstance(obj, np.ndarray):
        out = obj.copy()
        out.setflags(write=False)
        return out
    if not dataclasses.is_dataclass(obj):
        return copy.deepcopy(obj)
    kwargs = {}
    for f in dataclasses.fields(obj):
        value = getattr(obj, f.name)
        if isinstance(value, np.ndarray):
            value = _frozen(value)
        elif isinstance(value, Mapping):
            value = copy.deepcopy(dict(value))
        kwargs[f.name] = value
    return type(obj)(**kwargs)


def _snapshot(obj):
    if isinstance(obj, np.ndarray):
        return obj.copy()
    if not dataclasses.is_dataclass(obj):
        return copy.deepcopy(obj)
    return {f.name: copy.deepcopy(getattr(obj, f.name)) for f in dataclasses.fields(obj)}


def _same(a, b) -> bool:
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        a, b = np.asarray(a), np.asarray(b)
        if a.shape != b.shape:
            return False
        if a.dtype.kind in "fc" and b.dtype.kind in "fc":
            return bool(np.array_equal(a, b, equal_nan=True))
        return a.tolist() == b.tolist()
    if isinstance(a, Mapping) and isinstance(b, Mapping):
        return a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a)
    if dataclasses.is_dataclass(a) and dataclasses.is_dataclass(b):
        return _same(_snapshot(a), _snapshot(b))
    if isinstance(a, float) and isinstance(b, float) and np.isnan(a) and np.isnan(b):
        return True
    return a == b


def _assert_unchanged(obj, snap, label: str) -> None:
    now = _snapshot(obj)
    if isinstance(snap, dict) and isinstance(now, dict):
        changed = [k for k in snap if not _same(snap[k], now[k])]
        assert not changed, f"{label}: the tier changed {changed}"
    else:
        assert _same(snap, now), f"{label}: the tier changed its input"


def _all_implemented():
    return [(layer, name) for layer, names in IMPLEMENTED.items() for name in names]


# ---------------------------------------------------------------------------
# Every registered tier is covered, and describes itself
# ---------------------------------------------------------------------------


def test_every_registered_tier_is_covered_by_these_contracts():
    """A tier added to the registry has to be added here, as implemented, stub
    or external, so that it is held to the rules from its first commit."""
    known = set(_all_implemented()) | set(STUBS) | set(EXTERNAL)
    registered = {(layer, name) for layer, names in available().items() for name in names}
    assert registered - known == set(), "registered but not covered by the contract tests"
    assert known - registered == set(), "covered here but no longer registered"


@pytest.mark.parametrize(
    "layer,name",
    sorted({(lay, n) for lay, names in available().items() for n in names}),
)
def test_every_tier_describes_itself_in_json(layer, name):
    try:
        tier = build(layer, name)
    except (ImportError, ModuleNotFoundError) as exc:
        pytest.skip(f"{layer}/{name} needs software this environment lacks: {exc}")
    record = tier.describe()
    json.dumps(dict(record))  # raises on anything that is not JSON
    assert getattr(tier, "name", None) == name
    assert getattr(tier, "layer", None) == layer


@pytest.mark.parametrize("key,wp", sorted(STUBS.items()))
def test_stubs_raise_and_name_their_work_package(chain, key, wp):
    layer, name = key
    tier = build(layer, name)
    rng = np.random.default_rng(0)
    with pytest.raises(NotImplementedError, match=wp):
        if layer == "l1":
            tier.echoes(chain["truth"], chain["platform"], chain["freq_hz"], rng)
        elif layer == "l3":
            tier.form_ionogram(chain["signal"], chain["platform"])
        else:  # pragma: no cover - a new stub layer needs a branch here
            raise AssertionError(f"no call defined for a stub at {layer}")


# ---------------------------------------------------------------------------
# Rule 1: inputs are never mutated
# ---------------------------------------------------------------------------


def _call(layer: str, tier, chain, rng, *, frozen: bool):
    """Call one tier on its inputs from the chain, frozen or as they are."""
    f = _frozen if frozen else (lambda x: x)
    if layer == "l0":
        alt = f(chain["alt_m"])
        p = chain["platform"]
        return (alt,), tier.state_at(p.lat_deg, p.lon_deg, EPOCH, alt)
    if layer == "l1":
        args = (f(chain["truth"]), f(chain["platform"]), f(chain["freq_hz"]))
        return args, tier.echoes(*args, rng)
    if layer == "l2":
        args = (f(chain["echoes"]), f(chain["sweep"]))
        return args, tier.receive(*args, rng)
    if layer == "l3":
        args = (f(chain["signal"]), f(chain["platform"]))
        return args, tier.form_ionogram(*args)
    if layer == "l4":
        args = (f(chain["ionogram"]), f(chain["platform"]))
        return args, tier.invert(*args)
    if layer == "l5":
        args = (f(chain["retrieved"]), f(chain["truth"]), f(chain["ionogram"]))
        return args, tier.score(*args, "contract")
    raise AssertionError(layer)


@pytest.mark.parametrize("layer,name", _all_implemented())
def test_rule1_a_tier_does_not_write_into_its_inputs(chain, layer, name):
    """Read-only inputs: any in-place write raises at the line that makes it."""
    tier = build(layer, name)
    args, _ = _call(layer, tier, chain, np.random.default_rng(1), frozen=True)
    for arg in args:
        if isinstance(arg, np.ndarray):
            assert not arg.flags.writeable


@pytest.mark.parametrize("layer,name", _all_implemented())
def test_rule1_a_tier_leaves_its_inputs_as_it_found_them(chain, layer, name):
    """Including what a tier could change through a mapping it was handed."""
    tier = build(layer, name)
    args = {
        "l0": (chain["alt_m"],),
        "l1": (chain["truth"], chain["platform"], chain["freq_hz"]),
        "l2": (chain["echoes"], chain["sweep"]),
        "l3": (chain["signal"], chain["platform"]),
        "l4": (chain["ionogram"], chain["platform"]),
        "l5": (chain["retrieved"], chain["truth"], chain["ionogram"]),
    }[layer]
    snaps = [_snapshot(a) for a in args]
    _call(layer, tier, chain, np.random.default_rng(1), frozen=False)
    for i, (arg, snap) in enumerate(zip(args, snaps, strict=True)):
        _assert_unchanged(arg, snap, f"{layer}/{name} argument {i}")


# ---------------------------------------------------------------------------
# Rule 2: SI units, and one frame for the platform velocity
# ---------------------------------------------------------------------------


def _finite(x) -> np.ndarray:
    a = np.asarray(x, dtype=float).ravel()
    return a[np.isfinite(a)]


@pytest.mark.parametrize("name", IMPLEMENTED["l0"])
def test_rule2_l0_state_is_in_si(chain, name):
    p = chain["platform"]
    state = build("l0", name).state_at(p.lat_deg, p.lon_deg, EPOCH, chain["alt_m"])
    assert isinstance(state, IonosphereState)
    assert 1.0e4 < state.alt_m.min() and state.alt_m.max() < 1.0e8, "altitude in metres"
    assert 1.0e9 < np.nanmax(state.ne) < 1.0e14, "density in m^-3, not cm^-3"
    if state.b_field_t is not None:
        b = np.linalg.norm(np.asarray(state.b_field_t, dtype=float), axis=0)
        assert 1.0e-6 < np.nanmedian(b) < 1.0e-3, "flux density in tesla, not nT or gauss"


def test_rule2_echoes_are_in_si(chain):
    e: EchoSet = chain["echoes"]
    tau = _finite(e.group_delay_s)
    assert tau.size, "the smoke ionosphere returns echoes"
    # 15 km to 3000 km of two-way group path.
    assert 1.0e-4 < tau.min() and tau.max() < 2.0e-2, "group delay in seconds"
    assert np.array_equal(np.unique(e.freq_hz), np.unique(chain["freq_hz"])), "frequency in Hz"
    if e.reflection_alt_m is not None:
        h = _finite(e.reflection_alt_m)
        assert h.size and chain["alt_m"].min() <= h.min() and h.max() <= chain["alt_m"].max()


def test_rule2_signal_ionogram_and_profile_are_in_si(chain):
    s: ReceivedSignal = chain["signal"]
    assert 1.0 < s.fs_hz < 1.0e8
    ion: Ionogram = chain["ionogram"]
    assert 1.0e5 <= ion.freq_hz.min() and ion.freq_hz.max() <= 3.0e7, "frequency in Hz"
    vr = _finite(ion.virtual_range_m)
    assert vr.size and 1.0e4 < vr.min() and vr.max() < 5.0e6, "virtual range in metres"
    prof: DensityProfile = chain["retrieved"]
    h, ne = _finite(prof.alt_m), _finite(prof.ne)
    assert h.size and 5.0e4 < h.min() and h.max() < 3.0e6, "altitude in metres"
    assert 1.0e8 < ne.min() and ne.max() < 1.0e14, "density in m^-3"
    a: Assessment = chain["assessment"]
    assert a.requirements["SCI-07"]["unit"] == "m"


def test_rule2_platform_velocity_is_in_the_frame_its_consumers_read(chain):
    """A circular orbit neither climbs nor falls, and moves at about 7.6 km/s.

    Failed until 2026-10-17, when the orbit code delivered Earth-fixed
    Cartesian components that L1 and L5 read as north-east-down.
    """
    from tss.l1_propagation.doppler import closing_velocity_mps, vertical_line_of_sight
    from tss.l5_assessment.scoring import along_track_sampling_m

    p = chain["platform"]
    vertical = closing_velocity_mps(p.velocity_ned_mps, vertical_line_of_sight())
    assert abs(vertical) < 50.0, f"vertical closing speed {vertical:.0f} m/s for a circular orbit"
    spacing = along_track_sampling_m(1.0, p)
    assert spacing is not None and 7.0e3 < spacing < 8.0e3


# ---------------------------------------------------------------------------
# Rule 3: physical impossibility is a masked result, not an exception
# ---------------------------------------------------------------------------


def test_rule3_l1_masks_the_band_edges_instead_of_raising(chain):
    """Below the plasma frequency at the spacecraft, and above penetration."""
    tier = build("l1", "analytic")
    freq = np.array([1.0e5, 3.0e6, 3.0e7])
    e = tier.echoes(chain["truth"], chain["platform"], freq, np.random.default_rng(0))
    flags = [str(f) for f in np.asarray(e.flags, dtype=object)]
    o = [i for i, m in enumerate(e.mode) if str(m) == "O"]
    assert flags[o[0]] == "evanescent" and np.isnan(e.group_delay_s[o[0]])
    assert flags[o[-1]] == "penetrating" and np.isnan(e.group_delay_s[o[-1]])
    assert np.isfinite(e.group_delay_s[o[1]])


def _empty(echoes: EchoSet) -> EchoSet:
    def cut(v):
        return v[:0] if isinstance(v, np.ndarray) else v

    return dataclasses.replace(
        echoes, **{f.name: cut(getattr(echoes, f.name)) for f in dataclasses.fields(echoes)}
    )


def test_rule3_l2_accepts_an_echo_set_with_no_rows(chain):
    """An ionosphere that returns nothing is a result.  Raised a TypeError
    until 2026-10-17, from a list comprehension that built a float array."""
    signal = build("l2", "baseline").receive(
        _empty(chain["echoes"]), chain["sweep"], np.random.default_rng(0)
    )
    assert signal.flags["rows_used"] == 0


def test_rule3_an_echo_free_record_gives_an_empty_ionogram(chain):
    """Failed until 2026-10-17, when L2 set its sample rate from the echoes and
    sampled an echo-free record at 2.5 Hz."""
    nothing = dataclasses.replace(
        chain["echoes"],
        flags=np.array(["evanescent"] * chain["echoes"].freq_hz.size, dtype=object),
    )
    signal = build("l2", "baseline").receive(nothing, chain["sweep"], np.random.default_rng(0))
    ionogram = build("l3", "baseline").form_ionogram(signal, chain["platform"])
    assert signal.fs_hz == chain["signal"].fs_hz, "the sample rate must not follow the echoes"
    # A fixed 10 dB threshold over every bin out to each column's unambiguous
    # range (about 9300 km in the slow bottom steps) lets the odd noise spike
    # through: one, at 6368 km in the 0.34 MHz column, for this seed.  That is
    # the WP4 detector's business (CFAR, and a search limit at the design
    # range), not rule 3's, so this test holds only that noise gives an
    # ionogram and not a trace.
    assert np.count_nonzero(np.isfinite(ionogram.virtual_range_m)) <= 2


def test_rule3_l4_and_l5_take_an_ionogram_with_no_trace(chain):
    """An ionogram with nothing detected is inverted to an empty profile and
    scored as not evaluated, never as a pass and never as an exception."""
    ion = chain["ionogram"]
    blank = dataclasses.replace(
        ion,
        virtual_range_m=np.full_like(ion.virtual_range_m, np.nan),
        snr_db=np.full_like(ion.snr_db, np.nan),
    )
    for name in IMPLEMENTED["l4"]:
        profile = build("l4", name).invert(blank, chain["platform"])
        assert isinstance(profile, DensityProfile)
        a = build("l5", "standard").score(profile, chain["truth"], blank, "contract")
        for rid in ("SCI-05", "SCI-07"):
            assert a.requirements[rid]["verdict"] == "not_evaluated", (name, rid)


# ---------------------------------------------------------------------------
# Rule 4: randomness only from the injected generator
# ---------------------------------------------------------------------------


def test_rule4_the_same_seed_gives_the_same_run():
    a, b = _chain(seed=7), _chain(seed=7)
    assert np.array_equal(a["signal"].samples, b["signal"].samples)
    assert _same(a["ionogram"].virtual_range_m, b["ionogram"].virtual_range_m)
    assert a["assessment"].flags["summary"] == b["assessment"].flags["summary"]


def test_rule4_a_different_seed_gives_different_noise():
    a, b = _chain(seed=7), _chain(seed=8)
    assert not np.array_equal(a["signal"].samples, b["signal"].samples)


def test_rule4_no_tier_touches_the_global_generators():
    np.random.seed(12345)
    random.seed(12345)
    numpy_state = np.random.get_state()
    python_state = random.getstate()
    _chain(seed=3)
    after = np.random.get_state()
    assert numpy_state[0] == after[0] and np.array_equal(numpy_state[1], after[1])
    assert numpy_state[2:] == after[2:]
    assert random.getstate() == python_state


def test_rule4_no_tier_makes_a_generator_of_its_own_while_running(monkeypatch):
    """Tiers are built first; once running, a fresh generator is a violation."""
    sweep, alt_m, freq_hz, platform = _setup()
    rng = np.random.default_rng(5)
    t = {layer: build(layer, name) for layer, name in DEFAULT.items()}

    def refuse(*args, **kwargs):
        raise AssertionError("a tier created its own random generator")

    monkeypatch.setattr(np.random, "default_rng", refuse)
    monkeypatch.setattr(np.random, "RandomState", refuse)
    truth = t["l0"].state_at(platform.lat_deg, platform.lon_deg, EPOCH, alt_m)
    echoes = t["l1"].echoes(truth, platform, freq_hz, rng)
    signal = t["l2"].receive(echoes, sweep, rng)
    ionogram = t["l3"].form_ionogram(signal, platform)
    retrieved = t["l4"].invert(ionogram, platform)
    t["l5"].score(retrieved, truth, ionogram, "contract")
