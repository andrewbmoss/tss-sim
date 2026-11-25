"""WP5 Phase 3, first item: ``sheet_echo`` wrapped in the facet tier (D4's second half).

Exit for this item (``Plans\\WP5_plan.md`` P3): a noise-free sequence reproduces
the ``sheet_echo`` phase to 0.01 rad and the Doppler to 1e-3, and
``field_aligned_verify`` is unchanged.  ``sheet_echo`` is flat, so the exit is
met on the straight track; the curved sheet rests on the ray engine's own
checks (``test_l1_sphere``, ``test_l1_facet_3d``).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss import sweep_law
from tss.constants import C_LIGHT
from tss.interfaces import IonosphereState, PlatformState
from tss.l1_propagation import plate_echo as pe
from tss.l1_propagation import sheet_echo as se
from tss.l2_instrument.baseline import BaselineInstrument
from tss.pipeline import Pipeline
from tss.platform.orbits import CircularOrbit, platform_state_at
from tss.registry import build

F = 5.0e6
PROF = pe.ExponentialTopside()
SHEET = se.FieldAlignedSheet(dip_rad=math.radians(30.0), rotation_rad=math.radians(10.0),
                             width_m=1.3e3, length_m=20e3, curvature_radius_m=2.0e6,
                             velocity_n_mps=20.0, velocity_u_mps=5.0)
H, HD = 250e3, math.radians(40.0)
SPEC = {"id": "s1", "kind": "sheet", "dip_deg": 30.0, "rotation_deg": 10.0, "width_m": 1.3e3,
        "length_m": 20e3, "curvature_radius_m": 2.0e6, "velocity_n_mps": 20.0, "velocity_u_mps": 5.0,
        "height_m": H, "heading_deg": 40.0}


def _state():
    alt = np.linspace(50e3, 700e3, 1301)
    return IonosphereState(alt_m=alt, ne=(1.5e6**2 * np.exp((600e3 - alt) / 200e3)) / 80.6163)


def test_echo_track_reproduces_both_sheet_models_exactly():
    ts = np.linspace(-1.0, 1.0, 7)
    a = se.echo_stratified(ts, F, SHEET, PROF, height_m=H, heading_rad=HD)
    b = se.echo_track(ts, F, SHEET, PROF, height_m=H, heading_rad=HD)
    assert np.array_equal(a.field, b.field) and np.array_equal(a.group_path_m, b.group_path_m)
    c = se.echo_free_space(ts, F, SHEET, height_m=H, heading_rad=HD)
    d = se.echo_track(ts, F, SHEET, None, height_m=H, heading_rad=HD)
    assert np.array_equal(c.field, d.field) and np.array_equal(c.group_path_m, d.group_path_m)
    # The field's phase sits a quarter cycle behind the shortest ray's, the
    # stationary phase of a two-dimensional minimum, within a quarter cycle:
    # the facet tier centres its cycle there so the phase delay never jumps.
    k0 = 2 * np.pi * F / C_LIGHT
    res = np.angle(b.field * np.exp(2j * k0 * b.phase_path_m) * np.exp(0.5j * np.pi))
    assert np.max(np.abs(res)) < np.pi / 2


def test_straight_sheet_rows_are_the_sheet_model():
    tier = build("l1", "facet", facets=[SPEC])
    ts = np.array([-0.4, 0.0, 0.25])
    ref = se.echo_stratified(ts, F, SHEET, PROF, height_m=H, heading_rad=HD)
    h = 0.5e-3
    for t, e_ref, g_ref in zip(ts, ref.field, ref.group_path_m, strict=True):
        e = tier.echoes(_state(), PlatformState(alt_m=600e3, lat_deg=0, lon_deg=0, t_s=float(t)),
                        np.array([F]), np.random.default_rng(0))
        assert str(e.flags[0]) == "facet:s1;sheet"
        assert e.group_delay_s[0] == pytest.approx(2 * g_ref / C_LIGHT, rel=1e-12)
        dphi = np.angle(np.exp(1j * (2 * np.pi * F * e.phase_delay_s[0] + np.angle(e_ref))))
        assert abs(dphi) < 1e-6
        pair = se.echo_stratified(np.array([t - h, t + h]), F, SHEET, PROF, height_m=H, heading_rad=HD)
        want = np.angle(pair.field[1] * np.conj(pair.field[0])) / (2 * np.pi * 2 * h)
        assert e.doppler_hz[0] == pytest.approx(want, rel=1e-6, abs=1e-6)


def _carrier_phase(rec, echo, sweep):
    """As in test_run_sequence: 2 pi f tau_p at the record centre and its geometry time."""
    x = rec.samples[0]
    fs = rec.fs_hz
    tt = np.arange(x.size) / fs
    tc = 0.5 * x.size / fs
    i = 0
    tau = float(echo.group_delay_s[i])
    fd = float(echo.doppler_hz[i])
    k = sweep.k_f_hz_per_s
    live = np.abs(x) > 0
    z = np.sum(x[live] * np.exp(-2j * np.pi * (k * tau - fd) * (tt[live] - tc)))
    beat = 2 * np.pi * (sweep.f_start_hz * tau + k * tau * tc - k * tau**2 / 2 - F * tau)
    return float(np.angle(z)) - beat, tc - tau / 2  # the geometry time: reflected half a delay earlier


# For the exit the platform crosses the sheet's normal plane square on, so the
# Doppler at the centre is the drift's and the phase steps of 45 Hz sweeps do
# not alias (along the normal's azimuth the platform's own Doppler is about
# 100 Hz, which the drift mode samples only modulo 45 Hz; the echo set carries
# it exactly, test_straight_sheet_rows_are_the_sheet_model).
_N_AZ = math.degrees(math.atan2(math.sin(SHEET.rotation_rad),
                                math.cos(SHEET.rotation_rad) * math.sin(SHEET.dip_rad)))
HD_X = math.radians(_N_AZ + 90.0)


def test_exit_noise_free_sequence_reproduces_the_sheet_phase_and_doppler():
    spec = dict(SPEC, heading_deg=math.degrees(HD_X))
    tiers = {"l0": build("l0", "chapman"), "l1": build("l1", "facet", facets=[spec]),
             "l2": BaselineInstrument(quantise=False, thermal_noise=False)}
    sweep = sweep_law.drift_mode(F)
    n, i0 = 41, 20
    t = sweep_law.drift_mode_times(n, t0_s=-i0 / 45 - 1 / 90)
    plats = [PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0) for _ in t]
    res = Pipeline(tiers, seed=5).run_sequence("sheet", sweep=sweep, t_s=t, platforms=plats,
                                               alt_m=np.linspace(100e3, 600e3, 101),
                                               epoch="2026-01-01T12:00:00Z")
    errs, phs, tcs = [], [], []
    for rec, echo, tn in zip(res.sequence.records, res.echoes, t, strict=True):
        phase, tc = _carrier_phase(rec, echo, sweep)
        ref = se.echo_stratified(np.array([tn + tc]), F, SHEET, PROF, height_m=H, heading_rad=HD_X)
        errs.append(abs(np.angle(np.exp(1j * (phase + np.angle(ref.field[0]))))))
        phs.append(phase)
        tcs.append(tn + tc)
    assert max(errs) < 0.01
    # The Doppler from the L2 phase steps about the middle sweep, against the
    # sheet model's own phase over the same interval, sampled finely.
    ph = np.unwrap(phs)
    T = sweep.sweep_period_s
    est = -(ph[i0 + 1] - ph[i0 - 1]) / (2 * np.pi * 2 * T)
    fine = se.echo_stratified(np.linspace(tcs[i0 - 1], tcs[i0 + 1], 41), F, SHEET, PROF, height_m=H,
                              heading_rad=HD_X)
    arg = np.unwrap(np.angle(fine.field))
    want = (arg[-1] - arg[0]) / (2 * np.pi * (tcs[i0 + 1] - tcs[i0 - 1]))
    assert est == pytest.approx(want, rel=1e-3)


ORBIT = CircularOrbit(600e3, 60.0, u0_deg=20.0)


def _plat(t):
    p = platform_state_at(ORBIT, t)
    return PlatformState(alt_m=p.alt_m, lat_deg=p.lat_deg, lon_deg=p.lon_deg,
                         velocity_ned_mps=p.velocity_ned_mps, t_s=t, epoch="2026-01-01T12:00:00Z")


def test_orbit_horizontal_sheet_is_the_3d_plate():
    """With the field overridden horizontal (dip 0, declination 0) and no
    rotation, the sheet is a horizontal rectangle along magnetic (= true) north:
    the 3-D plate of the same size, laid along north, gives the same echo."""
    p0 = _plat(0.0)
    common = {"lat_deg": p0.lat_deg, "lon_deg": p0.lon_deg, "alt_m": 300e3}
    sheet = dict(common, id="s", kind="sheet", dip_deg=0.0, declination_deg=0.0, rotation_deg=0.0,
                 width_m=1.5e3, length_m=2.0e3, curvature_radius_m=math.inf, tukey_fraction=0.0,
                 velocity_n_mps=25.0)
    plate = dict(common, id="p", length_m=2.0e3, width_m=1.5e3, azimuth_deg=0.0, drift_mps=25.0)
    ts_ = build("l1", "facet", geometry="orbit", ref_alt_m=600e3, facets=[sheet])
    tp = build("l1", "facet", geometry="orbit", ref_alt_m=600e3, plate_model="3d", facets=[plate])
    st = _state()
    for t in (-0.3, 0.0, 0.2):
        a = ts_.echoes(st, _plat(t), np.array([F]), np.random.default_rng(0))
        b = tp.echoes(st, _plat(t), np.array([F]), np.random.default_rng(0))
        assert abs(np.angle(np.exp(2j * np.pi * F * (a.phase_delay_s[0] - b.phase_delay_s[0])))) < 0.01
        assert a.doppler_hz[0] == pytest.approx(b.doppler_hz[0], rel=1e-3, abs=0.02)


def test_orbit_sheet_orientation_override_and_normal_perpendicular_to_b():
    p0 = _plat(0.0)
    spec = {"id": "s", "kind": "sheet", "lat_deg": p0.lat_deg, "lon_deg": p0.lon_deg, "alt_m": 300e3,
            "dip_deg": 45.0, "declination_deg": 10.0, "rotation_deg": 20.0, "length_m": 10e3}
    tier = build("l1", "facet", geometry="orbit", ref_alt_m=600e3, facets=[spec])
    g = tier.sheet_geometry("s", _state())
    assert g["field_model"] == "override"
    assert abs(float(np.dot(g["b_local"], g["normal_local"]))) < 1e-12
    # b: 45 degrees down, horizontal part 10 degrees east of north.
    assert g["b_local"] == pytest.approx([math.cos(math.radians(45)) * math.sin(math.radians(10)),
                                          math.cos(math.radians(45)) * math.cos(math.radians(10)),
                                          -math.sin(math.radians(45))], abs=1e-12)


def test_orbit_sheet_field_from_igrf_south_points_up():
    pytest.importorskip("ppigrf")
    spec = {"id": "s", "kind": "sheet", "lat_deg": -35.0, "lon_deg": 140.0, "alt_m": 300e3}
    tier = build("l1", "facet", geometry="orbit", ref_alt_m=600e3, facets=[spec])
    g = tier.sheet_geometry("s", IonosphereState(alt_m=np.array([0.0, 1.0]), ne=np.ones(2),
                                                 epoch="2026-01-01T12:00:00Z"))
    assert g["field_model"] == "IGRF via ppigrf"
    assert g["dip_deg"] < -40.0  # southern hemisphere: the field points up
    assert g["b_local"][2] > 0.0
    assert tier.describe()["facets"][0]["field_model"] == "IGRF via ppigrf"


def test_orbit_sheet_needs_an_epoch_for_igrf():
    p0 = _plat(0.0)
    spec = {"id": "s", "kind": "sheet", "lat_deg": p0.lat_deg, "lon_deg": p0.lon_deg, "alt_m": 300e3}
    tier = build("l1", "facet", geometry="orbit", ref_alt_m=600e3, facets=[spec])
    with pytest.raises(ValueError):
        tier.sheet_geometry("s", IonosphereState(alt_m=np.array([0.0, 1.0]), ne=np.ones(2)))
