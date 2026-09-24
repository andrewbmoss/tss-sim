"""Layer: L0
WP: WP1
Claims: C1
Status: implemented (map evaluation); not yet an L0 tier
Requirements: -

Evaluate IRTAM coefficients into maps of foF2, hmF2, B0 and B1.

What it does.  :mod:`tss.data.irtam` fetches and parses the coefficient files
the Lowell GIRO Data Center publishes every 15 minutes.  This module turns one
of them into values at any set of geographic points, which is what an
event-driven L0 truth ionosphere needs: the F2 peak and the bottomside shape
that actually existed at a given UT, rather than IRI's monthly median.

The expansion.  IRTAM keeps the CCIR foF2 map formalism of Jones and Gallet
(1965) and Jones and Graham (1966), with one extra temporal term.  A value at
UT ``t`` and position ``x`` is

    value(t, x) = D(t) . U . G(x)

with ``U`` the ``(14, 76)`` coefficient matrix from
:attr:`tss.data.irtam.IrtamCoefficients.matrix`, ``D`` the 14 temporal
functions and ``G`` the 76 spatial functions.

* Temporal.  With ``h = 15 t - 180`` degrees and ``t`` in hours UT, the CCIR
  functions are ``1, sin h, cos h, sin 2h, cos 2h, ..., sin 6h, cos 6h``.
  IRTAM inserts a linear-trend term as the second function, equal to the
  minutes from the coefficient set's time of validity (TOV) plus 720.  At the
  TOV itself the trend function is 720.
* Spatial.  ``sin(modip)^i cos(lat)^j cos(j lon)`` and ``... sin(j lon)``, for
  longitude order ``j`` from 0 to 8 and ``i`` from 0 to ``Q_j - 1`` with
  ``Q = (12, 12, 9, 5, 2, 1, 1, 1, 1)``, and no sine term for ``j = 0``: 76
  functions.  ``modip`` is Rawer's modified dip latitude from
  :func:`tss.l0_ionosphere.geomag.modip_deg`.

Where this came from and how it was checked.  The basis, its ordering, the
trend term and the clipping of B0 and B1 follow PyIRTAM 0.0.7 and PyIRI 0.1.7
(Forsythe et al. 2023, 2024, Space Weather), which implement IRI's CCIR
formalism in Python and IRTAM on top of it.  The tests reproduce PyIRTAM's
matrix product to rounding error on the real coefficient files in
``tests/fixtures/irtam/``.  That shows this code computes what PyIRTAM
computes.  It does not by itself show that PyIRTAM matches Lowell's own
evaluation; the validation against GIRO measurements in
``docs/irtam-validation.md`` is the evidence for that.

Modip.  PyIRI evaluates IGRF-13 at 300 km, mid-month, treating geographic
latitude as geocentric.  :func:`~tss.l0_ionosphere.geomag.modip_deg` uses
IGRF-14 through ppigrf with a geodetic conversion.  The two differ by up to
about half a degree of modip; the effect on the maps is measured in the
validation note.  The dipole fallback is refused by default, because its modip
error reaches tens of degrees in the South Atlantic anomaly.

Time.  Coefficients are published every 15 minutes, each with its own TOV.
:func:`irtam_at` evaluates the file whose epoch is nearest the requested time,
at the requested time, so the offset is at most 7.5 minutes and the trend term
carries it.  Evaluating one file far from its TOV is refused unless asked for,
since nothing documents how far the trend term can be trusted.

What this is not yet.  It is not an L0 tier.  An L0 tier needs a full profile:
IRTAM gives the peak and the bottomside shape, and the topside above hmF2 has
to come from somewhere else (IRI's topside option, or NeQuick).  Building that
tier is the next step, and it is where the topside-model sensitivity the
proposal worries about enters.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from datetime import date, datetime, timezone
from typing import Any

import numpy as np

__all__ = [
    "Q_F0F2",
    "N_SPATIAL",
    "N_TEMPORAL",
    "B0_LIMITS_KM",
    "B1_LIMITS",
    "spatial_functions",
    "temporal_functions",
    "evaluate",
    "IrtamMaps",
    "irtam_at",
    "global_grid",
    "foF2_to_NmF2",
]

Q_F0F2: tuple[int, ...] = (12, 12, 9, 5, 2, 1, 1, 1, 1)
"""Highest power of sin(modip), plus one, for each longitude order 0 to 8."""

N_SPATIAL = 76
N_TEMPORAL = 14
N_TEMPORAL_CCIR = 13

B0_LIMITS_KM = (1.0, 350.0)
B1_LIMITS = (1.0, 6.0)
"""Clipping applied to B0 and B1, as PyIRTAM and IRI apply it, to keep a
poorly constrained expansion from producing unphysical bottomside shapes."""

MAX_OFFSET_MIN = 7.5
"""Largest distance from a coefficient set's TOV that :func:`evaluate` accepts
without ``allow_extrapolation``: half the 15-minute publication cadence."""

_UNITS = {"foF2": "MHz", "hmF2": "km", "B0": "km", "B1": ""}


# ---------------------------------------------------------------------------
# Basis functions
# ---------------------------------------------------------------------------


def spatial_functions(lat_deg: Any, lon_deg: Any, modip_deg: Any) -> np.ndarray:
    """The 76 Jones and Gallet spatial functions, shape ``(76, N)``.

    Points are flattened in C order.  The ordering matches PyIRI's
    ``set_global_functions``: longitude order outermost, then the power of
    sin(modip), then cosine before sine.
    """
    lat = np.radians(np.asarray(lat_deg, dtype=float).reshape(-1))
    lon = np.radians(np.asarray(lon_deg, dtype=float).reshape(-1))
    mod = np.radians(np.asarray(modip_deg, dtype=float).reshape(-1))
    if not (lat.size == lon.size == mod.size):
        raise ValueError("lat, lon and modip must have the same number of points")

    smod = np.sin(mod)
    clat = np.cos(lat)
    g = np.empty((N_SPATIAL, lat.size))
    k = 0
    for j, q in enumerate(Q_F0F2):
        cl = clat**j
        cj, sj = np.cos(j * lon), np.sin(j * lon)
        for i in range(q):
            base = smod**i * cl
            g[k] = base * cj
            k += 1
            if j != 0:
                g[k] = base * sj
                k += 1
    assert k == N_SPATIAL
    return g


def temporal_functions(ut_hours: float, trend_minutes: float | None = None) -> np.ndarray:
    """The temporal functions at one UT.

    With ``trend_minutes`` None, the 13 CCIR functions.  Otherwise the 14 IRTAM
    functions, with ``trend_minutes`` (minutes from the TOV plus 720) inserted
    second.
    """
    h = np.radians(15.0 * float(ut_hours) - 180.0)
    d = [1.0]
    for n in range(1, 7):
        d += [np.sin(n * h), np.cos(n * h)]
    d = np.array(d[:N_TEMPORAL_CCIR])
    if trend_minutes is None:
        return d
    return np.concatenate([d[:1], [float(trend_minutes)], d[1:]])


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------


def _as_utc(value: str | date | datetime) -> datetime:
    if isinstance(value, datetime):
        return value.astimezone(timezone.utc) if value.tzinfo else value.replace(tzinfo=timezone.utc)
    if isinstance(value, date):
        return datetime(value.year, value.month, value.day, tzinfo=timezone.utc)
    dt = datetime.fromisoformat(str(value).strip().replace("Z", "+00:00"))
    return dt.astimezone(timezone.utc) if dt.tzinfo else dt.replace(tzinfo=timezone.utc)


def _ut_hours(t: datetime) -> float:
    return t.hour + t.minute / 60.0 + (t.second + t.microsecond / 1e6) / 3600.0


def evaluate(
    coeffs: Any,
    lat_deg: Any,
    lon_deg: Any,
    *,
    modip_deg: Any,
    time: str | datetime | None = None,
    clip: bool = True,
    allow_extrapolation: bool = False,
) -> np.ndarray:
    """Evaluate one coefficient set at points, returning an array shaped like ``lat_deg``.

    ``coeffs`` is an :class:`~tss.data.irtam.IrtamCoefficients`.  ``time``
    defaults to the set's time of validity.  ``modip_deg`` must be supplied,
    so that the caller decides where the field model comes from and a map is
    never evaluated on a silently substituted dipole.

    B0 and B1 are clipped to :data:`B0_LIMITS_KM` and :data:`B1_LIMITS` when
    ``clip`` is true, as PyIRTAM and IRI clip them.
    """
    if (coeffs.n_temporal, coeffs.n_spatial) != (N_TEMPORAL, N_SPATIAL):
        raise ValueError(
            f"expected a {N_TEMPORAL} x {N_SPATIAL} IRTAM basis, got "
            f"{coeffs.n_temporal} x {coeffs.n_spatial}"
        )
    tov = coeffs.valid_time
    t = tov if time is None else _as_utc(time)
    offset_min = (t - tov).total_seconds() / 60.0
    if abs(offset_min) > MAX_OFFSET_MIN and not allow_extrapolation:
        raise ValueError(
            f"requested time is {offset_min:+.1f} min from the coefficient set's TOV "
            f"({tov:%Y-%m-%d %H:%M} UT); use the set nearest the time, or pass "
            "allow_extrapolation=True"
        )

    shape = np.shape(np.asarray(lat_deg))
    g = spatial_functions(lat_deg, lon_deg, modip_deg)
    d = temporal_functions(_ut_hours(t), offset_min + 720.0)
    values = (d @ coeffs.matrix) @ g

    if clip and coeffs.characteristic == "B0":
        values = np.clip(values, *B0_LIMITS_KM)
    elif clip and coeffs.characteristic == "B1":
        values = np.clip(values, *B1_LIMITS)
    return values.reshape(shape)


def foF2_to_NmF2(foF2_mhz: Any) -> np.ndarray:
    """Peak electron density in m^-3 from foF2 in MHz: ``Nm = 1.24e10 foF2^2``."""
    f = np.asarray(foF2_mhz, dtype=float)
    return 1.2404e10 * f * f


@dataclass(frozen=True)
class IrtamMaps:
    """IRTAM values at a set of points and one time.

    Every array has the shape of the ``lat`` passed in.  ``epochs`` names the
    coefficient set used for each characteristic, ``modip_source`` records the
    field model, and ``sources`` the cached coefficient files, for the
    manifest.
    """

    time: datetime
    lat_deg: np.ndarray
    lon_deg: np.ndarray
    modip_deg: np.ndarray
    values: dict[str, np.ndarray]
    epochs: dict[str, datetime]
    modip_source: dict[str, Any]
    sources: dict[str, str] = field(default_factory=dict)

    def __getitem__(self, key: str) -> np.ndarray:
        return self.values[key]

    @property
    def units(self) -> dict[str, str]:
        return {k: _UNITS[k] for k in self.values}

    @property
    def NmF2(self) -> np.ndarray:
        return foF2_to_NmF2(self.values["foF2"])


def irtam_at(
    time: str | datetime,
    lat_deg: Any,
    lon_deg: Any,
    *,
    characteristics: Sequence[str] = ("foF2", "hmF2", "B0", "B1"),
    modip_prefer: str = "igrf",
    route: str = "auto",
    clip: bool = True,
    fetch: Any = None,
) -> IrtamMaps:
    """IRTAM foF2, hmF2, B0 and B1 at points and a time, fetching as needed.

    Uses the coefficient set nearest ``time`` for each characteristic, through
    :func:`tss.data.irtam.fetch_coefficients` and its cache, and modip from
    IGRF at 300 km at mid-month, as IRI and PyIRI do.

    ``fetch`` replaces the fetch function; it takes ``(characteristic, time,
    route=...)`` and returns ``(IrtamCoefficients, path)``.  Tests use it to
    stay offline.
    """
    from tss.data import irtam as irtam_data
    from tss.l0_ionosphere.geomag import modip_deg

    t = _as_utc(time)
    lat = np.asarray(lat_deg, dtype=float)
    lon = np.asarray(lon_deg, dtype=float)
    lat, lon = np.broadcast_arrays(lat, lon)
    mid_month = datetime(t.year, t.month, 15, tzinfo=timezone.utc)
    modip, modip_source = modip_deg(lat, lon, epoch=mid_month.isoformat(), prefer=modip_prefer)

    get = fetch or irtam_data.fetch_coefficients
    values: dict[str, np.ndarray] = {}
    epochs: dict[str, datetime] = {}
    sources: dict[str, str] = {}
    for ch in characteristics:
        coeffs, path = get(ch, irtam_data.round_to_cadence(t), route=route)
        values[ch] = evaluate(coeffs, lat, lon, modip_deg=modip, time=t, clip=clip)
        epochs[ch] = coeffs.valid_time
        sources[ch] = str(path)
    return IrtamMaps(
        time=t, lat_deg=lat, lon_deg=lon, modip_deg=modip, values=values,
        epochs=epochs, modip_source=modip_source, sources=sources,
    )


def global_grid(dlat_deg: float = 2.5, dlon_deg: float = 5.0) -> tuple[np.ndarray, np.ndarray]:
    """A regular latitude-longitude grid for global maps, as two 2-D arrays.

    Latitudes run from -90 to 90 inclusive, longitudes from 0 up to but not
    including 360.  The default spacing is finer than IRTAM's own 46 x 45
    assimilation grid, which is the right direction for plotting and the wrong
    direction for claiming resolution the model does not have.
    """
    lats = np.arange(-90.0, 90.0 + 1e-9, dlat_deg)
    lons = np.arange(0.0, 360.0 - 1e-9, dlon_deg)
    lon2, lat2 = np.meshgrid(lons, lats)
    return lat2, lon2
