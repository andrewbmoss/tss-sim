"""Tests for IRTAM map evaluation (tss.l0_ionosphere.irtam_maps).

The reference values below were computed on 2026-09-23 with PyIRTAM 0.0.7 /
PyIRI 0.1.7 (``IRTAM_diurnal_functions``, ``IRTAM_read_files``,
``IRTAM_set_gl_G`` and the matrix product in ``IRTAM_gamma``) on the two real
coefficient files in ``tests/fixtures/irtam/``, at fixed points with fixed
modip.  They pin this implementation to PyIRTAM without needing PyIRTAM
installed.  No test here needs the network or ppigrf.
"""

from __future__ import annotations

from datetime import timedelta
from pathlib import Path

import numpy as np
import pytest

FIX = Path(__file__).resolve().parents[2] / "tests" / "fixtures" / "irtam"
ARCHIVE = FIX / "IRTAM_foF2_COEFFS_20121017_233000.ASC"
SERVLET = FIX / "servlet_foF2_20240301T1200.txt"

LAT = np.array([42.6, 0.0, -11.95, 13.6, -45.0, 70.0])
LON = np.array([288.5, 0.0, 283.13, 144.9, 170.0, 20.0])
MODIP = np.array([53.6, -26.0, -1.4, 12.0, -60.0, 72.0])

PYIRTAM_REF = {
    ARCHIVE: [6.1839073142, 6.9275991261, 10.6198300331, 10.8602744513, 6.3728182616, 2.5504336976],
    SERVLET: [6.7941956247, 13.2701879314, 9.5837543271, 14.0568668056, 6.2438092808, 10.106207473],
}


def _coeffs(path):
    from tss.data.irtam import parse_coefficients

    return parse_coefficients(path.read_bytes())


@pytest.mark.parametrize("path", [ARCHIVE, SERVLET])
def test_evaluate_reproduces_pyirtam(path) -> None:
    from tss.l0_ionosphere.irtam_maps import evaluate

    got = evaluate(_coeffs(path), LAT, LON, modip_deg=MODIP)
    np.testing.assert_allclose(got, PYIRTAM_REF[path], rtol=0, atol=1e-9)


def test_spatial_functions_count_order_and_constant_term() -> None:
    from tss.l0_ionosphere.irtam_maps import N_SPATIAL, spatial_functions

    g = spatial_functions(LAT, LON, MODIP)
    assert g.shape == (N_SPATIAL, LAT.size)
    np.testing.assert_array_equal(g[0], 1.0)
    np.testing.assert_allclose(g[1], np.sin(np.radians(MODIP)))
    # First longitude-order-1 pair: cos(lat) cos(lon), cos(lat) sin(lon).
    np.testing.assert_allclose(g[12], np.cos(np.radians(LAT)) * np.cos(np.radians(LON)))
    np.testing.assert_allclose(g[13], np.cos(np.radians(LAT)) * np.sin(np.radians(LON)))


def test_temporal_functions_ccir_and_irtam_forms() -> None:
    from tss.l0_ionosphere.irtam_maps import temporal_functions

    d13 = temporal_functions(12.0)
    assert d13.shape == (13,)
    # At 12 UT, h = 0: 1, then sin(nh) = 0 and cos(nh) = 1.
    np.testing.assert_allclose(d13, [1] + [0, 1] * 6, atol=1e-15)
    d14 = temporal_functions(12.0, 720.0)
    assert d14.shape == (14,) and d14[1] == 720.0
    np.testing.assert_array_equal(np.delete(d14, 1), d13)


def test_evaluate_defaults_to_the_time_of_validity_and_refuses_distant_times() -> None:
    from tss.l0_ionosphere.irtam_maps import evaluate

    c = _coeffs(SERVLET)
    at_tov = evaluate(c, LAT, LON, modip_deg=MODIP, time=c.valid_time)
    np.testing.assert_array_equal(at_tov, evaluate(c, LAT, LON, modip_deg=MODIP))
    evaluate(c, LAT, LON, modip_deg=MODIP, time=c.valid_time + timedelta(minutes=7))
    with pytest.raises(ValueError, match="TOV"):
        evaluate(c, LAT, LON, modip_deg=MODIP, time=c.valid_time + timedelta(minutes=30))
    far = evaluate(c, LAT, LON, modip_deg=MODIP, time=c.valid_time + timedelta(minutes=30),
                   allow_extrapolation=True)
    assert np.all(np.isfinite(far))


def test_evaluate_keeps_the_input_shape() -> None:
    from tss.l0_ionosphere.irtam_maps import evaluate, global_grid

    lat, lon = global_grid(30.0, 60.0)
    modip = np.clip(lat, -80, 80)
    v = evaluate(_coeffs(ARCHIVE), lat, lon, modip_deg=modip)
    assert v.shape == lat.shape == (7, 6)


def test_b0_and_b1_are_clipped_like_iri(monkeypatch) -> None:
    from dataclasses import replace

    from tss.l0_ionosphere.irtam_maps import B0_LIMITS_KM, B1_LIMITS, evaluate

    c = _coeffs(SERVLET)  # a foF2 file, relabelled to exercise the clip
    for name, (lo, hi) in (("B0", B0_LIMITS_KM), ("B1", B1_LIMITS)):
        v = evaluate(replace(c, characteristic=name), LAT, LON, modip_deg=MODIP)
        assert v.min() >= lo and v.max() <= hi
        raw = evaluate(replace(c, characteristic=name), LAT, LON, modip_deg=MODIP, clip=False)
        np.testing.assert_allclose(raw, PYIRTAM_REF[SERVLET], atol=1e-9)


def test_irtam_at_uses_the_nearest_epoch_and_records_provenance() -> None:
    from tss.l0_ionosphere import irtam_maps as im

    c = _coeffs(SERVLET)
    asked: list[tuple] = []

    def fake_fetch(ch, t, route="auto"):
        asked.append((ch, t, route))
        return c, SERVLET

    maps = im.irtam_at(
        c.valid_time + timedelta(minutes=5), LAT, LON,
        characteristics=("foF2",), modip_prefer="dipole", fetch=fake_fetch,
    )
    assert asked[0][1] == c.valid_time
    assert maps.epochs["foF2"] == c.valid_time
    assert maps.modip_source["model"] == "centred dipole"
    assert maps["foF2"].shape == LAT.shape
    np.testing.assert_allclose(maps.NmF2, 1.2404e10 * maps["foF2"] ** 2)
    assert maps.units == {"foF2": "MHz"}


def test_dipole_is_never_the_silent_default() -> None:
    import inspect

    from tss.l0_ionosphere import geomag, irtam_maps

    assert inspect.signature(irtam_maps.irtam_at).parameters["modip_prefer"].default == "igrf"
    assert inspect.signature(geomag.modip_deg).parameters["prefer"].default == "igrf"
    assert "modip_deg" in inspect.signature(irtam_maps.evaluate).parameters


def test_modip_dipole_route_is_rawer_modip_of_the_dipole_inclination() -> None:
    from tss.l0_ionosphere.geomag import inclination_deg, modip_deg

    lat = np.array([-60.0, 0.0, 45.0])
    lon = np.array([10.0, 100.0, 250.0])
    inc, _ = inclination_deg(lat, lon, prefer="dipole")
    mod, rec = modip_deg(lat, lon, prefer="dipole")
    expect = np.degrees(np.arctan(np.radians(inc) / np.sqrt(np.cos(np.radians(lat)))))
    np.testing.assert_allclose(mod, expect)
    assert rec["model"] == "centred dipole" and rec["quantity"] == "modip (Rawer)"
    assert np.all(np.sign(inc) == np.sign(mod))


def test_modip_igrf_route_raises_rather_than_falling_back(monkeypatch) -> None:
    import builtins

    from tss.l0_ionosphere.geomag import modip_deg

    real_import = builtins.__import__

    def no_ppigrf(name, *a, **k):
        if name == "ppigrf":
            raise ImportError("ppigrf deliberately hidden")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_ppigrf)
    with pytest.raises(ImportError):
        modip_deg(np.array([0.0]), np.array([0.0]))
    _m, rec = modip_deg(np.array([0.0]), np.array([0.0]), prefer="auto")
    assert rec["model"] == "centred dipole" and "reason" in rec


def test_physical_anchor_millstone_hill_at_noon_ut() -> None:
    """IRTAM at an assimilated station should sit close to that station's own measurement.

    Coefficients: the real 2024-03-01 12:00 UT foF2 file in the fixtures.
    Measurement: Millstone Hill (MHJ45) DIDBase records captured on 2026-09-23,
    11:52:30 to 12:07:30 UT, foF2 6.575, 6.588 and 6.675 MHz with confidence
    scores 80, 70 and 90.  Modip 53.6 degrees is IGRF-14 at 300 km, mid-March
    2024.  A wrong basis ordering or coefficient layout misses by megahertz,
    not by a tenth.
    """
    from tss.l0_ionosphere.irtam_maps import evaluate

    measured = np.mean([6.575, 6.588, 6.675])
    got = evaluate(_coeffs(SERVLET), np.array([42.6]), np.array([288.5]), modip_deg=np.array([53.6]))[0]
    assert abs(got - measured) < 0.3, got


@pytest.mark.parametrize("epoch,expect", [
    ("2024-03-15T00:00:00+00:00", "2024-03-15T00:00:00"),
    ("2024-03-15T00:00:00Z", "2024-03-15T00:00:00"),
    ("2024-03-15T06:00:00+06:00", "2024-03-15T00:00:00"),
    ("2024-03-15", "2024-03-15T00:00:00"),
])
def test_epoch_parser_returns_naive_utc(epoch, expect) -> None:
    """Found live on 2026-09-23: ppigrf raises on an aware datetime, and
    irtam_at passes an aware mid-month epoch."""
    from tss.l0_ionosphere.geomag import _parse_epoch_to_datetime

    dt = _parse_epoch_to_datetime(epoch)
    assert dt.tzinfo is None and dt.isoformat() == expect


def test_modip_through_ppigrf_with_an_aware_epoch() -> None:
    """The IGRF route end to end, where ppigrf is installed."""
    pytest.importorskip("ppigrf")
    from tss.l0_ionosphere.geomag import modip_deg

    m, rec = modip_deg(np.array([42.6, 0.0, 90.0]), np.array([288.5, 0.0, 0.0]),
                       epoch="2024-03-15T00:00:00+00:00")
    assert rec["model"] == "IGRF via ppigrf" and rec["version"]
    assert abs(m[0] - 53.615) < 0.01
    assert np.all(np.isfinite(m)), "the poles must not come back NaN"
