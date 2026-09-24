"""Layer: L0
WP: WP1
Claims: C2, C5
Status: implemented
Requirements: INS-08

Geomagnetic field for the ``b_field_t`` slot of
:class:`~tss.interfaces.IonosphereState`.

L0 owns the field because two things downstream need it and neither should be
inventing its own.  The magnetoionic splitting in L1 needs the field to place
the X trace, and the X-mode cutoff used by the band-edge study needs it to say
where the bottom of the band actually is.  A tier that leaves ``b_field_t`` as
``None`` is usable for a link-budget sweep and is not usable for a scored run.

Two routes, in order of preference.

**ppigrf.**  A pure-Python IGRF evaluation.  It is an optional dependency and
is imported lazily inside :func:`field_ned_t`, so ``import tss`` does not need
it.  When it is importable and answers, it is used, and ``describe()`` records
``"ppigrf"`` together with the package version.

**A tilted centred dipole.**  Used when ppigrf is absent or raises.  The dipole
axis is placed through the geomagnetic north pole at 80.65 N, 72.68 W, which is
the IGRF-13 epoch-2020 centred-dipole pole, and the equatorial surface field is
taken as 3.12e-5 T.  The field is

    B = B0 (a/r)^3 [ n - 3 (n . r) r ]

with ``n`` the unit vector toward the geomagnetic north pole and ``r`` the unit
position vector, then resolved into local north, east and down.

What the dipole costs is worth stating plainly.  Against IGRF it is wrong by up
to about 20 per cent in magnitude and by several degrees in direction, and it is
worst over the South Atlantic where the real field departs furthest from a
dipole.  For the O and X trace separation that error moves the separation by the
same fraction, which is tens of kilohertz at HF and is visible on an ionogram.
So the dipole is a fallback that keeps a run possible, not a fallback that keeps
a run publishable, and ``describe()`` says which one produced a given result
precisely so that a reviewer can tell the two apart.

The field is returned in local north-east-down components, tesla, as an array of
shape ``(3, N)`` parallel to the altitude grid, which is what
:class:`~tss.interfaces.IonosphereState` documents.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from tss.constants import R_EARTH_M

__all__ = [
    "B0_EQUATORIAL_T",
    "DIPOLE_POLE_LAT_DEG",
    "DIPOLE_POLE_LON_DEG",
    "dipole_field_ned_t",
    "field_ned_t",
    "geomagnetic_latitude_deg",
    "inclination_deg",
    "modip_deg",
]

#: Equatorial surface field of the centred dipole, tesla.
B0_EQUATORIAL_T = 3.12e-5

#: Centred-dipole north pole, IGRF-13 at epoch 2020.
DIPOLE_POLE_LAT_DEG = 80.65
DIPOLE_POLE_LON_DEG = -72.68


def _unit_from_lat_lon(lat_deg: float, lon_deg: float) -> np.ndarray:
    """Unit vector in geocentric Cartesian components from a spherical pair."""
    lat = math.radians(float(lat_deg))
    lon = math.radians(float(lon_deg))
    return np.array(
        [math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)],
        dtype=float,
    )


def geomagnetic_latitude_deg(lat_deg: float, lon_deg: float) -> float:
    """Latitude in the centred-dipole frame, degrees.

    Reported in ``extra`` by the tiers that use the dipole, because the topside
    is organised by magnetic latitude and a diagnostic keyed on geographic
    latitude alone hides that organisation.
    """
    axis = _unit_from_lat_lon(DIPOLE_POLE_LAT_DEG, DIPOLE_POLE_LON_DEG)
    point = _unit_from_lat_lon(lat_deg, lon_deg)
    return math.degrees(math.asin(float(np.clip(np.dot(axis, point), -1.0, 1.0))))


def dipole_field_ned_t(
    lat_deg: float,
    lon_deg: float,
    alt_m: np.ndarray | float,
) -> np.ndarray:
    """Tilted centred-dipole field in local north-east-down components, tesla.

    Returns an array of shape ``(3, N)`` for an altitude grid of ``N`` samples.
    The altitude is treated as a height above a spherical Earth of radius
    :data:`tss.constants.R_EARTH_M`.  At the altitudes this package works at, the
    difference between that and a geodetic treatment is smaller than the error
    the dipole approximation already carries, so it is not worth carrying an
    ellipsoid here.
    """
    alt = np.atleast_1d(np.asarray(alt_m, dtype=float))
    r = R_EARTH_M + alt

    axis = _unit_from_lat_lon(DIPOLE_POLE_LAT_DEG, DIPOLE_POLE_LON_DEG)
    r_hat = _unit_from_lat_lon(lat_deg, lon_deg)

    lat = math.radians(float(lat_deg))
    lon = math.radians(float(lon_deg))
    e_north = np.array(
        [-math.sin(lat) * math.cos(lon), -math.sin(lat) * math.sin(lon), math.cos(lat)],
        dtype=float,
    )
    e_east = np.array([-math.sin(lon), math.cos(lon), 0.0], dtype=float)
    e_down = -r_hat

    scale = B0_EQUATORIAL_T * (R_EARTH_M / r) ** 3          # (N,)
    direction = axis - 3.0 * float(np.dot(axis, r_hat)) * r_hat   # (3,)

    b_cart = np.outer(direction, scale)                     # (3, N)
    return np.vstack(
        [
            e_north @ b_cart,
            e_east @ b_cart,
            e_down @ b_cart,
        ]
    )


def _ppigrf_field_ned_t(
    lat_deg: float,
    lon_deg: float,
    alt_m: np.ndarray,
    epoch: str,
) -> np.ndarray:
    """IGRF through ppigrf, in north-east-down tesla.

    ppigrf returns east, north and up components in nanotesla, so the sign of
    the third component flips and all three are scaled.  The import is inside
    this function so that ``import tss`` does not require the package.
    """
    import datetime as _dt

    import ppigrf  # type: ignore

    when = _parse_epoch_to_datetime(epoch)
    if when is None:
        when = _dt.datetime(2020, 1, 1, 12, 0, 0)

    alt_km = np.asarray(alt_m, dtype=float) / 1.0e3
    lon = np.full(alt_km.shape, float(lon_deg))
    lat = np.full(alt_km.shape, float(lat_deg))
    be, bn, bu = ppigrf.igrf(lon, lat, alt_km, when)
    be = np.asarray(be, dtype=float).reshape(-1)[: alt_km.size]
    bn = np.asarray(bn, dtype=float).reshape(-1)[: alt_km.size]
    bu = np.asarray(bu, dtype=float).reshape(-1)[: alt_km.size]
    return np.vstack([bn, be, -bu]) * 1.0e-9


def _parse_epoch_to_datetime(epoch: str):
    """ISO 8601 to :class:`datetime.datetime`, or ``None`` if it will not parse.

    A field model that cannot read the epoch should fall back to a default date
    rather than fail the run, because the secular variation over a few years
    moves the field by far less than the dipole fallback already does.

    The result is always naive UTC.  ppigrf compares the date against a
    timezone-naive table and raises on an aware datetime, so an epoch carrying
    an offset (``+00:00``, as ``datetime.isoformat`` writes for an aware time)
    is converted to UTC and its offset dropped.
    """
    import datetime as _dt

    text = str(epoch or "").strip()
    if not text:
        return None
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    for fmt in (None, "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
        try:
            if fmt is None:
                dt = _dt.datetime.fromisoformat(text)
            else:
                dt = _dt.datetime.strptime(text, fmt)
        except (TypeError, ValueError):
            continue
        if dt.tzinfo is not None:
            dt = dt.astimezone(_dt.timezone.utc).replace(tzinfo=None)
        return dt
    return None


def field_ned_t(
    lat_deg: float,
    lon_deg: float,
    alt_m: np.ndarray | float,
    epoch: str = "",
    *,
    prefer: str = "auto",
) -> tuple[np.ndarray, dict[str, Any]]:
    """Geomagnetic field in north-east-down tesla, with a provenance record.

    Parameters
    ----------
    lat_deg, lon_deg
        Geodetic location.
    alt_m
        Altitude grid, metres.  Scalars are promoted to a one-element grid.
    epoch
        ISO 8601 UTC timestamp.  Used only by the IGRF route.
    prefer
        ``"auto"`` tries ppigrf and falls back to the dipole.  ``"igrf"``
        demands ppigrf and raises if it cannot be used, which is what a scored
        campaign should ask for.  ``"dipole"`` skips ppigrf entirely, which is
        what a test on a bare machine should ask for.

    Returns
    -------
    (b_ned, provenance)
        ``b_ned`` has shape ``(3, N)``.  ``provenance`` is JSON-serialisable and
        carries ``model``, ``version`` and, when the dipole was used because
        ppigrf failed rather than because it was not asked for, ``reason``.
    """
    alt = np.atleast_1d(np.asarray(alt_m, dtype=float))

    if prefer not in ("auto", "igrf", "dipole"):
        raise ValueError(f"unknown prefer {prefer!r}, expected auto, igrf or dipole")

    if prefer == "dipole":
        return dipole_field_ned_t(lat_deg, lon_deg, alt), {
            "model": "centred dipole",
            "version": "",
            "pole_lat_deg": DIPOLE_POLE_LAT_DEG,
            "pole_lon_deg": DIPOLE_POLE_LON_DEG,
            "b0_equatorial_t": B0_EQUATORIAL_T,
        }

    try:
        b = _ppigrf_field_ned_t(lat_deg, lon_deg, alt, epoch)
        import ppigrf  # type: ignore

        return b, {
            "model": "IGRF via ppigrf",
            "version": str(getattr(ppigrf, "__version__", "") or _package_version("ppigrf")),
        }
    except Exception as exc:                       # noqa: BLE001 - see below
        # Any failure of the optional route is the same thing to the caller:
        # the field has to come from somewhere else.  The reason is carried into
        # the manifest rather than swallowed, so a run that silently used the
        # dipole is not possible.
        if prefer == "igrf":
            raise
        return dipole_field_ned_t(lat_deg, lon_deg, alt), {
            "model": "centred dipole",
            "version": "",
            "pole_lat_deg": DIPOLE_POLE_LAT_DEG,
            "pole_lon_deg": DIPOLE_POLE_LON_DEG,
            "b0_equatorial_t": B0_EQUATORIAL_T,
            "reason": f"ppigrf unavailable or failed: {type(exc).__name__}: {exc}",
        }


# ---------------------------------------------------------------------------
# Vectorised inclination and modified dip, for global maps
# ---------------------------------------------------------------------------


def inclination_deg(
    lat_deg: np.ndarray | float,
    lon_deg: np.ndarray | float,
    *,
    alt_km: float = 300.0,
    epoch: str = "",
    prefer: str = "igrf",
) -> tuple[np.ndarray, dict[str, Any]]:
    """Magnetic inclination (dip) in degrees, positive downward, at many points.

    :func:`field_ned_t` answers for one location on an altitude grid, which is
    what a profile needs.  A global map needs the opposite shape, many
    locations at one altitude, and looping :func:`field_ned_t` over a few
    thousand grid points would call ppigrf a few thousand times.  This function
    makes one ppigrf call for all of them.

    ``prefer`` defaults to ``"igrf"`` here, not ``"auto"``, because the one
    consumer so far is the IRTAM and CCIR map evaluation, where modip enters
    every spatial basis function and the dipole's error of several degrees in
    modip moves the equatorial anomaly crests.  Ask for ``"dipole"`` explicitly
    to run on a machine without ppigrf, and expect maps that are qualitatively
    right and quantitatively not.
    """
    lat = np.asarray(lat_deg, dtype=float)
    lon = np.asarray(lon_deg, dtype=float)
    lat, lon = np.broadcast_arrays(lat, lon)
    shape = lat.shape
    flat_lat, flat_lon = lat.reshape(-1), lon.reshape(-1)

    if prefer not in ("auto", "igrf", "dipole"):
        raise ValueError(f"unknown prefer {prefer!r}, expected auto, igrf or dipole")

    def _dipole() -> np.ndarray:
        alt_m = np.array([alt_km * 1.0e3])
        out = np.empty(flat_lat.size)
        for i, (la, lo) in enumerate(zip(flat_lat, flat_lon, strict=True)):
            bn, be, bd = dipole_field_ned_t(float(la), float(lo), alt_m)[:, 0]
            out[i] = np.degrees(np.arctan2(bd, np.hypot(bn, be)))
        return out

    dipole_record = {
        "model": "centred dipole",
        "version": "",
        "alt_km": float(alt_km),
        "pole_lat_deg": DIPOLE_POLE_LAT_DEG,
        "pole_lon_deg": DIPOLE_POLE_LON_DEG,
    }
    if prefer == "dipole":
        return _dipole().reshape(shape), dipole_record

    try:
        import datetime as _dt

        import ppigrf  # type: ignore

        when = _parse_epoch_to_datetime(epoch) or _dt.datetime(2020, 1, 1, 12, 0, 0)
        h = np.full(flat_lat.shape, float(alt_km))
        # ppigrf divides by sin(colatitude) and returns NaN exactly at a pole.
        # A thousandth of a degree is about 100 m, far below anything a map at
        # degree resolution can resolve.
        safe_lat = np.clip(flat_lat, -89.999, 89.999)
        be, bn, bu = ppigrf.igrf(flat_lon, safe_lat, h, when)
        be = np.asarray(be, dtype=float).reshape(-1)[: flat_lat.size]
        bn = np.asarray(bn, dtype=float).reshape(-1)[: flat_lat.size]
        bd = -np.asarray(bu, dtype=float).reshape(-1)[: flat_lat.size]
        inc = np.degrees(np.arctan2(bd, np.hypot(bn, be)))
        record = {
            "model": "IGRF via ppigrf",
            "version": _package_version("ppigrf"),
            "alt_km": float(alt_km),
            "epoch": when.isoformat(),
        }
        return inc.reshape(shape), record
    except Exception as exc:  # noqa: BLE001 - same policy as field_ned_t
        if prefer == "igrf":
            raise
        rec = dict(dipole_record)
        rec["reason"] = f"ppigrf unavailable or failed: {type(exc).__name__}: {exc}"
        return _dipole().reshape(shape), rec


def modip_deg(
    lat_deg: np.ndarray | float,
    lon_deg: np.ndarray | float,
    *,
    alt_km: float = 300.0,
    epoch: str = "",
    prefer: str = "igrf",
) -> tuple[np.ndarray, dict[str, Any]]:
    """Rawer's modified dip latitude, degrees, at many points.

    ``tan(modip) = I / sqrt(cos(lat))`` with the inclination ``I`` in radians.
    This is the coordinate the CCIR, URSI and IRTAM foF2 and hmF2 maps are
    expanded in.  The defaults, IGRF at 300 km, are the ones IRI and PyIRI use;
    PyIRI also evaluates the field at mid-month, which a caller reproduces by
    passing that date as ``epoch``.
    """
    inc, record = inclination_deg(lat_deg, lon_deg, alt_km=alt_km, epoch=epoch, prefer=prefer)
    lat = np.radians(np.broadcast_to(np.asarray(lat_deg, dtype=float), inc.shape))
    coslat = np.clip(np.cos(lat), 1.0e-12, None)
    mod = np.degrees(np.arctan(np.radians(inc) / np.sqrt(coslat)))
    return mod, dict(record, quantity="modip (Rawer)")


def _package_version(name: str) -> str:
    try:
        import importlib.metadata as md

        return md.version(name)
    except Exception:  # noqa: BLE001
        return ""
