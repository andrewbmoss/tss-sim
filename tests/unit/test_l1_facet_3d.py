"""WP5 Phase 3: the facet tier in three dimensions (Andrew, 2026-11-25, decisions (a), (b)).

Every scatterer's ray is solved in its own plane through the Earth's centre
(``tss.l1_propagation.sphere``, held to Bouguer's law in ``test_l1_sphere``).
Here the three-dimensional surfaces are held to the verified flat sheet model
in the flat limit, a tilted plate to the sheet it equals, and the 3-D plate on
the orbit to P2's 2-D plate up to the constant the cross-track integral adds.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss.constants import C_LIGHT
from tss.interfaces import IonosphereState, PlatformState
from tss.l1_propagation import plate_echo as pe
from tss.l1_propagation import sheet_echo as se
from tss.l1_propagation.facet import pair_phase_delay_difference, plate_surface
from tss.l1_propagation.sphere import SphericalMedium, ray_paths
from tss.platform.orbits import RE, CircularOrbit, platform_state_at
from tss.registry import build

F = 5.0e6
K0 = 2 * np.pi * F / C_LIGHT
PROF = pe.ExponentialTopside()
FLAT_R0 = 1.0e15
ORBIT = CircularOrbit(600e3, 60.0, u0_deg=20.0)


def _engine(pts_local, w, plat_local, base_height, f=F):
    med = SphericalMedium(PROF, base_height)
    ph = ray_paths(plat_local, pts_local, f, med, FLAT_R0)[0]
    return np.sum(w * np.exp(-2j * 2 * np.pi * f / C_LIGHT * ph))


def _mag_to_local(v):
    """Magnetic (north, east, up) to local (east, north, up) with no declination."""
    v = np.asarray(v, dtype=float)
    return v[..., [1, 0, 2]]


def test_engine_in_the_flat_limit_reproduces_the_flat_sheet_model():
    """The 3-D ray engine, on a sphere of 1e15 m, against sheet_echo.echo_stratified
    with a dipping, rotated, curved, moving sheet and an oblique track."""
    sheet = se.FieldAlignedSheet(dip_rad=math.radians(30.0), rotation_rad=math.radians(10.0),
                                 width_m=1.3e3, length_m=20e3, curvature_radius_m=2.0e6,
                                 velocity_n_mps=20.0, velocity_u_mps=5.0)
    h, hd = 250e3, math.radians(40.0)
    ts = np.linspace(-1.0, 1.0, 9)
    ref = se.echo_stratified(ts, F, sheet, PROF, height_m=h, heading_rad=hd)
    pts, w = se.sheet_points(sheet, F)
    p0 = se.platform_start_stratified(sheet, h, PROF, F)
    vel, th = se.sheet_velocity(sheet), se.heading_vector(hd)
    for t, e_ref in zip(ts, ref.field, strict=True):
        plat = p0 + (7600.0 * th - vel) * t  # the medium frame, as sheet_echo
        e = _engine(_mag_to_local(pts), w, _mag_to_local(plat), p0[2])
        assert abs(np.angle(e * np.conj(e_ref))) < 0.01
        assert abs(e) == pytest.approx(abs(e_ref), rel=0.01, abs=1e-4)


@pytest.mark.parametrize("tilt_deg", [0.0, 1.0])
def test_a_tilted_plate_is_the_sheet_it_equals(tilt_deg):
    """A uniform rectangle tilted along north is a sheet with that dip, no
    rotation, no curvature and no taper: two independent point layouts."""
    a = math.radians(tilt_deg)
    plate = pe.Plate(length_m=2.0e3, tilt_rad=a)
    pts, w, _ = plate_surface(plate, width_m=1.5e3, azimuth_rad=0.0, cross_tilt_rad=0.0, f_hz=F)
    sheet = se.FieldAlignedSheet(dip_rad=a, width_m=1.5e3, length_m=2.0e3, tukey_fraction=0.0)
    spts, sw = se.sheet_points(sheet, F)
    h = 300e3
    # Along-track positions inside the main lobe, about the specular point: in
    # the sidelobes the field is the edges' and the two layouts' end weights
    # (quarter against half a wavelength apart) differ there by design.
    ys = h * math.tan(a)
    for y in (ys - 1.5e3, ys, ys + 1.0e3):  # along-track position of the platform, north
        plat = np.array([400.0, y, h])
        e1 = _engine(pts, w, plat, h)
        e2 = _engine(_mag_to_local(spts), sw, plat, h)
        assert abs(np.angle(e1 * np.conj(e2))) < 0.01
        assert abs(e1) == pytest.approx(abs(e2), rel=0.02)


def test_plate_surface_centre_and_normal():
    plate = pe.Plate(length_m=2.0e3, tilt_rad=math.radians(2.0))
    pts, w, normal = plate_surface(plate, width_m=1.0e3, azimuth_rad=math.radians(90.0),
                                   cross_tilt_rad=math.radians(1.0), f_hz=F)
    assert pts[0] == pytest.approx([0.0, 0.0, 0.0])  # the centre, weight zero
    assert w[0] == 0.0 and w.sum() == pytest.approx(1.0)
    assert np.max(np.abs((pts[1:] - pts[0]) @ normal)) < 1e-6  # all in the plane
    # Tilt along the azimuth (east here) leans the normal east; cross tilt
    # leans it toward the azimuth's right (south).
    assert normal[0] == pytest.approx(math.sin(math.radians(2.0)), rel=2e-3)
    assert normal[1] < 0.0


def _plat(t):
    p = platform_state_at(ORBIT, t)
    return PlatformState(alt_m=p.alt_m, lat_deg=p.lat_deg, lon_deg=p.lon_deg,
                         velocity_ned_mps=p.velocity_ned_mps, t_s=t)


def _state():
    alt = np.linspace(50e3, 700e3, 1301)
    return IonosphereState(alt_m=alt, ne=(1.5e6**2 * np.exp((600e3 - alt) / 200e3)) / 80.6163)


def _tier(model, *, medium="plasma", drift=25.0, length=2.0e3, **kw):
    p0 = _plat(0.0)
    facet = {"id": "o1", "lat_deg": p0.lat_deg, "lon_deg": p0.lon_deg, "alt_m": 300e3, "drift_mps": drift,
             "length_m": length}
    facet.update(kw)
    return build("l1", "facet", geometry="orbit", medium=medium, profile="exponential", ref_alt_m=600e3,
                 plate_model=model, facets=[facet])


def test_orbit_3d_plate_matches_the_2d_plate_up_to_a_constant():
    """Same group delay, same Doppler; the phase differs by a constant (the
    cross-track integral), so its spread over the lobe is far under 0.01 rad."""
    t2, t3 = _tier("2d"), _tier("3d", width_m=2.0e3)
    st = _state()
    ts = np.arange(-20, 21) / 45.0
    e2 = [t2.echoes(st, _plat(t), np.array([F]), np.random.default_rng(0)) for t in ts]
    e3 = [t3.echoes(st, _plat(t), np.array([F]), np.random.default_rng(0)) for t in ts]
    g2 = np.array([e.group_delay_s[0] for e in e2])
    g3 = np.array([e.group_delay_s[0] for e in e3])
    assert np.max(np.abs(g2 - g3)) * C_LIGHT / 2 < 0.05
    d2 = np.array([e.doppler_hz[0] for e in e2])
    d3 = np.array([e.doppler_hz[0] for e in e3])
    i0 = 20
    assert d3[i0] == pytest.approx(d2[i0], rel=1e-3)
    dphi = np.unwrap(2 * np.pi * F * np.array([b.phase_delay_s[0] - a.phase_delay_s[0]
                                               for a, b in zip(e2, e3, strict=True)]))
    lobe = np.array([e.amplitude_db[0] for e in e2]) > max(e.amplitude_db[0] for e in e2) - 20
    assert np.ptp(dphi[lobe]) < 0.01


def test_orbit_3d_free_space_range_includes_the_cross_track_offset():
    """A facet 3 km to the side of the track: the 3-D group path is the
    spherical chord to its centre, cross-track included."""
    p0 = _plat(0.0)
    v = np.asarray(p0.velocity_ned_mps, dtype=float)
    hd = math.atan2(v[1], v[0])
    # 3 km to the right of the track at the facet's height.
    dpsi = 3e3 / (RE + 300e3)
    lat = math.degrees(math.radians(p0.lat_deg) - dpsi * math.sin(hd))
    lon = p0.lon_deg + math.degrees(dpsi * math.cos(hd) / math.cos(math.radians(p0.lat_deg)))
    tier = build("l1", "facet", geometry="orbit", medium="free", ref_alt_m=600e3, plate_model="3d",
                 facets=[{"id": "x", "lat_deg": lat, "lon_deg": lon, "alt_m": 300e3}])
    for t in (-0.5, 0.0, 0.7):
        p = _plat(t)
        e = tier.echoes(_state(), p, np.array([F]), np.random.default_rng(0))
        a1, o1 = math.radians(p.lat_deg), math.radians(p.lon_deg)
        a2, o2 = math.radians(lat), math.radians(lon)
        u1 = np.array([math.cos(a1) * math.cos(o1), math.cos(a1) * math.sin(o1), math.sin(a1)])
        u2 = np.array([math.cos(a2) * math.cos(o2), math.cos(a2) * math.sin(o2), math.sin(a2)])
        chord = np.linalg.norm((RE + p.alt_m) * u1 - (RE + 300e3) * u2)
        assert e.group_delay_s[0] == pytest.approx(2 * chord / C_LIGHT, rel=1e-9)
    assert tier.max_cross_track_m == pytest.approx(3e3, rel=0.01)


def test_orbit_3d_overflight_doppler_is_the_drift_doppler():
    tier = _tier("3d")
    e = tier.echoes(_state(), _plat(0.0), np.array([F]), np.random.default_rng(0))
    mu = math.sqrt(1 - (1.5e6 / F) ** 2)
    from tss.l1_propagation.doppler import doppler_shift_hz

    assert e.doppler_hz[0] == pytest.approx(float(doppler_shift_hz(25.0, F, mu_sat=mu)), rel=2e-3)


def test_straight_3d_plate_is_the_flat_limit():
    """geometry='straight', plate_model='3d': the same rays on a sphere of 1e15 m.
    Its Doppler at overflight is the drift Doppler and its group delay the 2-D one."""
    t2 = build("l1", "facet", facets=[{"id": "f", "depth_m": 300e3, "drift_mps": 25.0}])
    t3 = build("l1", "facet", plate_model="3d", facets=[{"id": "f", "depth_m": 300e3, "drift_mps": 25.0}])
    st = _state()
    for t in (0.0, 0.3):
        p = PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0, t_s=t)
        a = t2.echoes(st, p, np.array([F]), np.random.default_rng(0))
        b = t3.echoes(st, p, np.array([F]), np.random.default_rng(0))
        assert b.group_delay_s[0] == pytest.approx(a.group_delay_s[0], rel=1e-9)
        assert b.doppler_hz[0] == pytest.approx(a.doppler_hz[0], rel=2e-3, abs=0.01)


def test_pair_difference_survives_a_wrap():
    """The Doppler pair is taken from the conjugate product, so a residual phase
    that crosses +/-pi between the two times does not jump by a cycle."""
    pc = np.array([1000.0, 1000.0 + 1e-3])
    res = np.array([np.pi - 0.01, -np.pi + 0.01])  # crosses the cut
    field = np.exp(1j * res) * np.exp(-2j * K0 * pc)
    d = pair_phase_delay_difference(field[0], field[1], pc[0], pc[1], F)
    want = 2 * ((pc[1] - pc[0]) - 0.02 / (2 * K0)) / C_LIGHT
    assert d == pytest.approx(want, rel=1e-9)


def test_plate_model_is_checked():
    with pytest.raises(ValueError):
        build("l1", "facet", plate_model="4d")
