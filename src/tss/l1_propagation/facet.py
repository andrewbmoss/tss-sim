"""Layer: L1
WP: WP5
Claims: C3
Status: partial
Requirements: SCI-15, INS-09, INS-15

The facet tier (``l1/facet``): the verified physical-optics plate and sheet
models as an L1 tier, one echo set per sweep (WP5 Phases 2 and 3, decision D4).

WHAT IT DOES
------------
For each sweep of a drift-mode sequence the pipeline calls :meth:`FacetTier.echoes`
with the platform state of that sweep, whose ``t_s`` is the time from the start
of the sequence.  For every facet and every sounding frequency the tier returns
one EchoSet row:

* ``group_delay_s``: twice the group path of the reference ray over ``c`` (a
  plate's centre; a sheet's shortest ray, as :mod:`sheet_echo` takes it);
* ``phase_delay_s``: the phase delay of the summed field, ``2 P_eff / c`` with
  ``P_eff = P_c - arg(E exp(+j 2 k0 P_c)) / (2 k0)``, ``P_c`` the reference
  ray's phase path; ``2 pi f tau_p`` is then ``-arg(E)`` to rounding, in L2's sign;
* ``doppler_hz``: ``-f d tau_p / dt`` over ``doppler_dt_s``, the difference
  taken from the conjugate product of the two fields so that it cannot jump by
  a cycle (closing positive, as :mod:`tss.l1_propagation.doppler`);
* ``amplitude_db``: the free-space spreading of the analytic tier over twice
  the group path, plus ``reflection_db``, plus ``20 log10 |E|``.  ``|E|`` is
  near 1 at the specular peak of a facet inside its Fresnel zone, so
  ``reflection_db`` sets the echo strength there.  A 3-D plate's ``|E|`` is not
  the 2-D plate's: the cross-track integral changes it;
* ``flags``: ``facet:<id>`` for a plate and ``facet:<id>;sheet`` for a sheet,
  the name L2 keys fading and sequence state on.

A facet the frequency cannot reach (below the platform's plasma frequency, or
at or past the reflection level) is a masked row, ``...;beyond_reflection``,
with NaN delays (contract rule 3).  The plate and sheet models are called
unchanged (D4).

FOUR WAYS TO COMPUTE A PLATE, ONE FOR A SHEET ON EACH TRACK
-----------------------------------------------------------
``geometry`` is ``"straight"`` (the flat models' own straight level track at
``v_s_mps``; the platform state supplies only the time) or ``"orbit"`` (the
platform's position and velocity from each platform state).  ``plate_model``
is ``"2d"`` (the plate model's line of points, as in Phase 2) or ``"3d"`` (a
rectangle of ``length_m`` by ``width_m``, every point's ray solved in its own
plane through the Earth's centre by :mod:`tss.l1_propagation.sphere`).  The
default, ``"auto"``, is 2-D on the straight track (M3 is scored there, like
for like with ``c3_verify``; decision (c)) and 3-D on the orbit (decision (b)).
The other two combinations are kept to split what the plate model does from
what the sphere does: straight 3-D is the same ray engine on a sphere of
``flat_radius_m`` (1e15 m), and orbit 2-D is Phase 2's flattened plate.

A sheet (``kind: sheet``) on the straight track is :func:`sheet_echo.echo_track`,
which reproduces the verified flat model exactly; on the orbit it is the 3-D
ray engine with the sheet laid out in the magnetic frame at its reference
point: dip and declination from IGRF (``field_model``, provenance in
:meth:`describe`) unless given, the field line's curvature from the dipole
formula of :mod:`sheet_echo` unless given.

STRAIGHT GEOMETRY (verified)
----------------------------
A plate is placed by ``depth_m`` (or ``mu_facet`` at ``mu_at_hz``) and ``x_m``
relative to the platform at ``t = 0`` (``y_m`` to the right, 3-D only); a sheet
by ``height_m`` below the platform and the platform's ``heading_deg`` from
magnetic north.  ``medium="plasma"`` with ``profile="exponential"`` is the
reference topside of Technical Volume Chapter 8 (``f_n_sat_hz`` 1.5 MHz,
``scale_height_m`` 200 km), the c3_verify profile.

ORBIT GEOMETRY AND THE L0 PROFILE (Phases 2 and 3)
--------------------------------------------------
A facet is placed by its sub-point (``lat_deg``, ``lon_deg``, spherical, as
:mod:`tss.platform.orbits`) and altitude (``alt_m``).  A 3-D plate lies along
``azimuth_deg`` (from north; by default the platform's heading at the first
call, latched and reported).  The 3-D rays are held to Bouguer's law
integrated in radius (``tests/unit/test_l1_sphere.py``), to the flat sheet
model in the flat limit and to the 2-D plate up to a constant
(``tests/unit/test_l1_facet_3d.py``).  What the sphere does to route one is
real: near overflight the range's second derivative is ``v^2 (r_f / r_p) / D``,
not the flat ``v^2 / D``, 4.3 % less at 300 km below a 600 km orbit
(``r_f / r_p = 0.957``; ``tests/unit/test_l1_facet_orbit.py``).  It moves
route one's slope ratio, the compact-facet test, by that factor; the drift
route one forms is, noise-free and to first order, the Doppler at the power
peak and does not move (measured, ``tests/unit/test_l3_drift.py``; 0.06 to
0.28 m/s with an asymmetric envelope; the Phase 2 note that it would be 4 %
wrong was corrected on 2026-11-25).

The tier keeps state between calls (the latched plate azimuth, the sheet's
field direction, the surfaces and the Newton warm start): build one per scene.

Drift is the strata moving as a whole: the profile rises radially at
``drift_mps`` and the plate with it (a sheet at the vertical part of its
velocity); ``drift_x_mps`` moves a plate along its azimuth.  The 2-D orbit
plate (Phase 2) is two-dimensional: it ignores the platform's cross-track
offset, whose largest value is in :attr:`FacetTier.max_cross_track_m`; the 3-D
plate does not.

``profile="l0"`` takes the electron density from the L0 state the pipeline
passes (``f_N^2 = 80.616 N_e``, interpolated in ``log N_e``), with
``ref_alt_m`` as the level depths are counted from.

THE LAYER ECHO
--------------
With ``layer_echo=True`` each frequency also gets the vertical echo of the
layer below, from the L0 state the pipeline passes, flagged ``layer:O``: the
field-free group path ``int dz / mu`` and phase path ``int mu dz`` from the
platform to the reflection level (by the substitution ``z = z_r - s^2``, which
takes the square-root singularity out), its spreading over twice the group
path, and the Doppler of the platform's vertical velocity.  No absorption, no
magnetoionic split (it is not the ``analytic`` tier's echo), no layer drift.
A frequency below the platform's plasma frequency or above the layer's
critical frequency gives a masked row.

THE DIFFUSE HALO (PHASE 3, DECISION (E))
----------------------------------------
A plate may carry ``halo: {k, ratio, n_points, seed}`` (:class:`HaloSpec`):
``n_points`` (200) weak point scatterers, equal in all directions, drawn
uniformly along the plate's line over ``ratio`` (2) plate lengths about its
centre, each with a random phase, moving with the medium (option A, Andrew,
2026-11-27).  Each gives its own row, flagged ``halo:<id>``, with its own
group and phase delay and Doppler, so L2 sums them as it sums echoes and the
facet's specular fading does not touch them.  ``K`` is the facet's power over
the halo's mean power **at overflight** (the platform over the facet's centre,
at its height): every scatterer has the amplitude ``a`` with
``n_points a^2 = |E_facet(overflight)|^2 / K``, spreading apart, so the halo's
power stays about constant through a pass while the facet's rises and falls.
``K = inf`` is no halo.  ``K = 0`` leaves the facet's row out and gives the
halo the power the facet would have had.  The field of such a patch
decorrelates as ``sinc(pi B tau)``, ``B = 2 v L / (lambda P')``
(:func:`tss.l2_instrument.fading.halo_decorrelation_s`).  The scatterers lie
on a line along the track, as the 2-D plate does: there is no cross-track
extent.  Not on a sheet, nor on the orbit's 2-D plate.

WHAT IS LEFT OUT
----------------
Everything the plate and sheet models leave out (the magnetic field's effect on
the wave, absorption, reflectivity variation across the surface); the halo's
scatterers do not interact with the facet or with each other.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from tss.constants import C_LIGHT
from tss.interfaces import EchoSet, IonosphereState, PlatformState
from tss.l1_propagation import sheet_echo as se
from tss.l1_propagation.plate_echo import (
    ExponentialTopside,
    Plate,
    echo_track,
    plate_points,
    point_paths,
)
from tss.l1_propagation.sphere import LocalFrame, SphericalMedium, ray_paths
from tss.platform.orbits import RE
from tss.registry import register

__all__ = [
    "FacetSpec",
    "SheetSpec",
    "FacetTier",
    "DEFAULT_FACET",
    "TabulatedTopside",
    "FlattenedProfile",
    "FN2_PER_NE",
    "plate_surface",
    "pair_phase_delay_difference",
]

#: The facet used when none is configured: a 2 km uniform plate at nadir,
#: 400 km below the platform, still.
DEFAULT_FACET: Mapping[str, Any] = {"id": "f1", "depth_m": 400e3}

#: ``f_N^2 / N_e``, Hz^2 m^3: ``e^2 / (4 pi^2 eps0 m_e)``.
FN2_PER_NE = 80.6163


class TabulatedTopside:
    """Field-free index from a tabulated density, depth counted from ``ref_alt_m``.

    ``mu^2 = 1 - f_N^2 / f^2`` with ``f_N^2 = 80.616 N_e`` and ``N_e``
    interpolated in ``log N_e`` on the state's altitude grid.
    """

    def __init__(self, alt_m: np.ndarray, ne: np.ndarray, ref_alt_m: float) -> None:
        alt = np.asarray(alt_m, dtype=float)
        n = np.clip(np.asarray(ne, dtype=float), 1.0, None)
        order = np.argsort(alt)
        self.alt_m = alt[order]
        self.log_ne = np.log(n[order])
        self.ref_alt_m = float(ref_alt_m)

    def f_n2(self, depth_m: np.ndarray | float) -> np.ndarray:
        d = np.asarray(depth_m, dtype=float)
        alt = self.ref_alt_m - d
        return FN2_PER_NE * np.exp(np.interp(alt.ravel(), self.alt_m, self.log_ne).reshape(d.shape))

    def mu_squared(self, f_hz: float, depth_m: np.ndarray | float) -> np.ndarray:
        return 1.0 - self.f_n2(depth_m) / float(f_hz) ** 2


class _FreeSpace:
    def mu_squared(self, f_hz: float, depth_m: np.ndarray | float) -> np.ndarray:
        return np.ones_like(np.asarray(depth_m, dtype=float))


class FlattenedProfile:
    """Earth-flattened index ``mu'^2(d) = mu^2(r) (r/R)^2`` (Phase 2's 2-D orbit plate).

    ``base.mu_squared(f, depth_true)`` is the index at true depth
    ``R - r`` below the reference radius ``R``.
    """

    def __init__(self, base: Any, r_ref_m: float) -> None:
        self.base = base
        self.r_ref_m = float(r_ref_m)

    def radius(self, d_flat: np.ndarray | float) -> np.ndarray:
        return self.r_ref_m * np.exp(-np.asarray(d_flat, dtype=float) / self.r_ref_m)

    def flat_depth(self, r_m: float) -> float:
        return float(self.r_ref_m * math.log(self.r_ref_m / float(r_m)))

    def mu_squared(self, f_hz: float, depth_m: np.ndarray | float) -> np.ndarray:
        r = self.radius(depth_m)
        return self.base.mu_squared(f_hz, self.r_ref_m - r) * (r / self.r_ref_m) ** 2


def pair_phase_delay_difference(e0: complex, e1: complex, pc0: float, pc1: float, f_hz: float) -> float:
    """``tau_p(1) - tau_p(0)`` from two fields and their reference phase paths.

    ``tau_p = 2 (P_c - res / (2 k0)) / c`` with ``res = arg(E exp(+j 2 k0 P_c))``;
    the difference of the residuals is taken as the angle of their conjugate
    product, so it is the short way round even when each residual sits near
    the cut at +/- pi.
    """
    k0 = 2.0 * math.pi * float(f_hz) / C_LIGHT
    dpc = float(pc1) - float(pc0)
    d_res = float(np.angle(complex(e1) * np.conj(complex(e0)) * np.exp(2j * k0 * dpc)))
    return 2.0 * (dpc - d_res / (2.0 * k0)) / C_LIGHT


def _plate_normal_and_along(plate: Plate, azimuth_rad: float,
                            cross_tilt_rad: float) -> tuple[np.ndarray, np.ndarray]:
    """A 3-D plate's unit normal and its along-track tangent, local (east, north, up)."""
    a_hat = np.array([math.sin(azimuth_rad), math.cos(azimuth_rad), 0.0])
    c_hat = np.array([math.cos(azimuth_rad), -math.sin(azimuth_rad), 0.0])
    up = np.array([0.0, 0.0, 1.0])
    n = math.tan(float(plate.tilt_rad)) * a_hat + math.tan(float(cross_tilt_rad)) * c_hat + up
    n /= np.linalg.norm(n)
    t_c = np.cross(a_hat, n)
    t_c /= np.linalg.norm(t_c)
    return n, np.cross(n, t_c)


def plate_surface(
    plate: Plate, *, width_m: float, azimuth_rad: float, cross_tilt_rad: float, f_hz: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A plate as a rectangle of points in a local (east, north, up) frame.

    Centred on the origin, ``plate.length_m`` along ``azimuth_rad`` (from north)
    and ``width_m`` across it, tilted by ``plate.tilt_rad`` toward the azimuth
    and by ``cross_tilt_rad`` toward its right.  The along-track samples and
    taper are :func:`plate_points`'s; across, a quarter wavelength, uniform.
    Returns the points with the centre first (weight zero, the reference ray),
    the weights (sum one) and the unit normal.
    """
    lam = C_LIGHT / float(f_hz)
    s, wa = plate_points(plate, f_hz)
    nc = max(16, int(np.ceil(float(width_m) / (lam / 4.0))) | 1)
    y = np.linspace(-float(width_m) / 2.0, float(width_m) / 2.0, nc)
    wc = np.full(nc, 1.0 / nc)
    a_hat = np.array([math.sin(azimuth_rad), math.cos(azimuth_rad), 0.0])
    n, t_a = _plate_normal_and_along(plate, azimuth_rad, cross_tilt_rad)
    t_c = np.cross(a_hat, n)
    t_c /= np.linalg.norm(t_c)
    pts = s[:, None, None] * t_a[None, None, :] + y[None, :, None] * t_c[None, None, :]
    w = np.outer(wa, wc).ravel()
    pts = np.vstack([np.zeros((1, 3)), pts.reshape(-1, 3)])
    return pts, np.concatenate([[0.0], w / w.sum()]), n


@dataclass(frozen=True)
class HaloSpec:
    """The diffuse halo around a plate (module docstring): ``K``, size, count and draw."""

    k: float
    ratio: float = 2.0
    n_points: int = 200
    seed: int = 0

    @classmethod
    def from_mapping(cls, m: Mapping[str, Any]) -> HaloSpec:
        d = dict(m)
        unknown = set(d) - {"k", "ratio", "n_points", "seed"}
        if unknown:
            raise ValueError(f"unknown halo keys {sorted(unknown)}")
        if "k" not in d:
            raise ValueError("a halo needs k, the facet's power over the halo's")
        spec = cls(k=float(d["k"]), ratio=float(d.get("ratio", 2.0)), n_points=int(d.get("n_points", 200)),
                   seed=int(d.get("seed", 0)))
        if not spec.k >= 0.0:
            raise ValueError("halo k must not be negative")
        if not (spec.ratio > 0.0 and math.isfinite(spec.ratio)):
            raise ValueError("halo ratio must be positive")
        if spec.n_points < 1:
            raise ValueError("halo n_points must be at least 1")
        return spec

    @property
    def present(self) -> bool:
        """False for ``K = inf``, which is no halo."""
        return math.isfinite(self.k)

    def draw(self, facet_length_m: float) -> tuple[np.ndarray, np.ndarray]:
        """Positions along the plate's line about its centre, m, and phases, rad."""
        rng = np.random.default_rng(self.seed)
        length = self.ratio * float(facet_length_m)
        s = rng.uniform(-0.5 * length, 0.5 * length, self.n_points)
        return s, rng.uniform(-math.pi, math.pi, self.n_points)


@dataclass(frozen=True)
class FacetSpec:
    """One plate: its id, its plate, its reflection (dB), its 3-D extent and, on an orbit, its place."""

    facet_id: str
    plate: Plate
    reflection_db: float = 0.0
    lat_deg: float | None = None
    lon_deg: float | None = None
    alt_m: float | None = None
    width_m: float | None = None
    cross_tilt_rad: float = 0.0
    y_m: float = 0.0
    azimuth_deg: float | None = None
    kind: str = "plate"
    halo: HaloSpec | None = None

    @classmethod
    def from_mapping(cls, m: Mapping[str, Any], profile: Any, *, orbit: bool) -> FacetSpec:
        """From a configuration mapping.

        Keys: ``id``; on the straight track ``depth_m``, or ``mu_facet`` with
        ``mu_at_hz`` (inside the plasma), and ``x_m`` (``y_m`` for a 3-D plate);
        on an orbit ``lat_deg``, ``lon_deg`` and ``alt_m`` (``azimuth_deg`` for
        a 3-D plate); and ``tilt_deg``, ``cross_tilt_deg``, ``length_m``,
        ``width_m`` (3-D; default the length), ``drift_mps``, ``drift_x_mps``,
        ``taper``, ``reflection_db``, ``halo`` (:class:`HaloSpec`).
        """
        d = dict(m)
        lat = lon = alt = None
        if orbit:
            if not all(k in d for k in ("lat_deg", "lon_deg", "alt_m")):
                raise ValueError("on an orbit a facet needs lat_deg, lon_deg and alt_m")
            lat, lon, alt = float(d["lat_deg"]), float(d["lon_deg"]), float(d["alt_m"])
            depth = 0.0  # set by the tier
        elif "depth_m" in d:
            depth = float(d["depth_m"])
        elif "mu_facet" in d:
            if profile is None or not hasattr(profile, "depth_for_mu"):
                raise ValueError("mu_facet places a facet only inside the exponential plasma")
            depth = float(profile.depth_for_mu(float(d["mu_at_hz"]), float(d["mu_facet"])))
        else:
            raise ValueError("a facet needs depth_m, or mu_facet and mu_at_hz")
        plate = Plate(
            x_m=float(d.get("x_m", 0.0)), depth_m=depth,
            tilt_rad=math.radians(float(d.get("tilt_deg", 0.0))),
            length_m=float(d.get("length_m", 2.0e3)), drift_mps=float(d.get("drift_mps", 0.0)),
            taper=str(d.get("taper", "uniform")), drift_x_mps=float(d.get("drift_x_mps", 0.0)),
        )
        width = d.get("width_m")
        az = d.get("azimuth_deg")
        return cls(facet_id=str(d.get("id", "f1")), plate=plate,
                   reflection_db=float(d.get("reflection_db", 0.0)), lat_deg=lat, lon_deg=lon, alt_m=alt,
                   width_m=None if width is None else float(width),
                   cross_tilt_rad=math.radians(float(d.get("cross_tilt_deg", 0.0))),
                   y_m=float(d.get("y_m", 0.0)), azimuth_deg=None if az is None else float(az),
                   halo=None if d.get("halo") is None else HaloSpec.from_mapping(d["halo"]))

    def describe(self) -> dict[str, Any]:
        p = self.plate
        halo = {} if self.halo is None else {"halo": {
            "k": self.halo.k, "ratio": self.halo.ratio, "length_m": self.halo.ratio * p.length_m,
            "n_points": self.halo.n_points, "seed": self.halo.seed}}
        return {**halo, "id": self.facet_id, "kind": "plate", "x_m": p.x_m, "depth_m": p.depth_m,
                "tilt_deg": math.degrees(p.tilt_rad), "length_m": p.length_m,
                "drift_mps": p.drift_mps, "drift_x_mps": p.drift_x_mps, "taper": p.taper,
                "reflection_db": self.reflection_db, "lat_deg": self.lat_deg,
                "lon_deg": self.lon_deg, "alt_m": self.alt_m, "width_m": self.width_m,
                "cross_tilt_deg": math.degrees(self.cross_tilt_rad), "y_m": self.y_m,
                "azimuth_deg": self.azimuth_deg}


@dataclass(frozen=True)
class SheetSpec:
    """One field-aligned sheet (``kind: sheet``); the :class:`sheet_echo.FieldAlignedSheet` keys.

    Straight track: ``height_m`` below the platform, the platform's
    ``heading_deg`` from magnetic north, ``dip_deg`` (default 0).  Orbit:
    ``lat_deg``, ``lon_deg``, ``alt_m``; ``dip_deg`` and ``declination_deg``
    each override IGRF.  ``curvature_radius_m`` is a number (``inf`` for a
    straight field line) or ``"dipole"`` (the orbit's default).
    """

    facet_id: str
    rotation_deg: float = 0.0
    width_m: float = 1.3e3
    length_m: float = 60e3
    curvature: Any = "dipole"
    tukey_fraction: float = 0.3
    velocity_n_mps: float = 0.0
    velocity_u_mps: float = 0.0
    reflection_db: float = 0.0
    dip_deg: float | None = None
    declination_deg: float | None = None
    height_m: float | None = None
    heading_deg: float = 0.0
    lat_deg: float | None = None
    lon_deg: float | None = None
    alt_m: float | None = None
    kind: str = "sheet"

    @classmethod
    def from_mapping(cls, m: Mapping[str, Any], *, orbit: bool) -> SheetSpec:
        d = dict(m)
        if d.get("halo") is not None:
            raise ValueError("a halo is modelled around a plate, not around a sheet")
        curv: Any = d.get("curvature_radius_m", "dipole" if orbit else math.inf)
        if isinstance(curv, str):
            if curv.lower() in ("inf", "flat", "straight"):
                curv = math.inf
            elif curv.lower() == "dipole":
                curv = "dipole"
            else:
                raise ValueError(f"curvature_radius_m: a number, inf or 'dipole', not {curv!r}")
        else:
            curv = float(curv)
        lat = lon = alt = None
        height = None
        if orbit:
            if not all(k in d for k in ("lat_deg", "lon_deg", "alt_m")):
                raise ValueError("on an orbit a sheet needs lat_deg, lon_deg and alt_m")
            lat, lon, alt = float(d["lat_deg"]), float(d["lon_deg"]), float(d["alt_m"])
        else:
            if "height_m" not in d:
                raise ValueError("on the straight track a sheet needs height_m")
            height = float(d["height_m"])
            if curv == "dipole":
                raise ValueError("on the straight track give curvature_radius_m as a number")
        dip = d.get("dip_deg")
        dec = d.get("declination_deg")
        return cls(
            facet_id=str(d.get("id", "s1")), rotation_deg=float(d.get("rotation_deg", 0.0)),
            width_m=float(d.get("width_m", 1.3e3)), length_m=float(d.get("length_m", 60e3)),
            curvature=curv, tukey_fraction=float(d.get("tukey_fraction", 0.3)),
            velocity_n_mps=float(d.get("velocity_n_mps", 0.0)),
            velocity_u_mps=float(d.get("velocity_u_mps", 0.0)),
            reflection_db=float(d.get("reflection_db", 0.0)),
            dip_deg=None if dip is None else float(dip), declination_deg=None if dec is None else float(dec),
            height_m=height, heading_deg=float(d.get("heading_deg", 0.0)),
            lat_deg=lat, lon_deg=lon, alt_m=alt,
        )

    def sheet(self, dip_rad: float, curvature_radius_m: float) -> se.FieldAlignedSheet:
        return se.FieldAlignedSheet(
            dip_rad=float(dip_rad), rotation_rad=math.radians(self.rotation_deg), width_m=self.width_m,
            length_m=self.length_m, curvature_radius_m=float(curvature_radius_m),
            tukey_fraction=self.tukey_fraction, velocity_n_mps=self.velocity_n_mps,
            velocity_u_mps=self.velocity_u_mps,
        )

    def describe(self) -> dict[str, Any]:
        curv = self.curvature if isinstance(self.curvature, str) else float(self.curvature)
        return {"id": self.facet_id, "kind": "sheet", "rotation_deg": self.rotation_deg,
                "width_m": self.width_m, "length_m": self.length_m, "curvature_radius_m": curv,
                "tukey_fraction": self.tukey_fraction, "velocity_n_mps": self.velocity_n_mps,
                "velocity_u_mps": self.velocity_u_mps, "reflection_db": self.reflection_db,
                "dip_deg": self.dip_deg, "declination_deg": self.declination_deg, "height_m": self.height_m,
                "heading_deg": self.heading_deg, "lat_deg": self.lat_deg, "lon_deg": self.lon_deg,
                "alt_m": self.alt_m}


def _track_angles(p: PlatformState, lat_deg: float, lon_deg: float) -> tuple[float, float]:
    """Along- and cross-track angles from the platform's sub-point to a point, rad.

    Along is positive ahead of the platform, along its horizontal velocity.
    """
    p1, l1 = math.radians(p.lat_deg), math.radians(p.lon_deg)
    p2, l2 = math.radians(lat_deg), math.radians(lon_deg)
    dl = l2 - l1
    cosd = math.sin(p1) * math.sin(p2) + math.cos(p1) * math.cos(p2) * math.cos(dl)
    delta = math.acos(max(-1.0, min(1.0, cosd)))
    bearing = math.atan2(math.sin(dl) * math.cos(p2),
                         math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl))
    v = np.asarray(p.velocity_ned_mps, dtype=float)
    heading = math.atan2(float(v[1]), float(v[0]))
    rel = bearing - heading
    xt = math.asin(max(-1.0, min(1.0, math.sin(delta) * math.sin(rel))))
    at = math.atan2(math.sin(delta) * math.cos(rel), math.cos(delta))
    return at, xt


@dataclass
class _Surface:
    """A facet's scatterers in its local frame at ``t = 0``, for one frequency."""

    pts: np.ndarray
    w: np.ndarray
    vel: np.ndarray
    rise_mps: float
    ref: str  # "centre" (index 0) or "shortest"


@register("l1", "facet")
class FacetTier:
    """Facet echoes per sweep from the plate and sheet models.  See the module docstring."""

    layer = "l1"
    name = "facet"

    def __init__(
        self,
        *,
        facets: Sequence[Mapping[str, Any]] | None = None,
        medium: str = "plasma",
        profile: str = "exponential",
        f_n_sat_hz: float = 1.5e6,
        scale_height_m: float = 200e3,
        geometry: str = "straight",
        v_s_mps: float = 7600.0,
        ref_alt_m: float | None = None,
        earth_radius_m: float = RE,
        doppler_dt_s: float = 1.0e-3,
        group_df_fraction: float = 1.0e-4,
        layer_echo: bool = False,
        plate_model: str = "auto",
        flat_radius_m: float = 1.0e15,
        field_model: str = "igrf",
    ) -> None:
        if medium not in ("plasma", "free"):
            raise ValueError("medium must be 'plasma' or 'free'")
        if profile not in ("exponential", "l0"):
            raise ValueError("profile must be 'exponential' or 'l0'")
        if geometry not in ("straight", "orbit"):
            raise ValueError("geometry must be 'straight' or 'orbit'")
        if plate_model not in ("auto", "2d", "3d"):
            raise ValueError("plate_model must be 'auto', '2d' or '3d'")
        if field_model not in ("igrf", "auto", "dipole"):
            raise ValueError("field_model must be 'igrf', 'auto' or 'dipole'")
        if not float(doppler_dt_s) > 0.0:
            raise ValueError("doppler_dt_s must be positive")
        if (geometry == "orbit" or (medium == "plasma" and profile == "l0")) and ref_alt_m is None:
            raise ValueError("the orbit geometry and the l0 profile need ref_alt_m")
        self.medium = medium
        self.profile_kind = profile if medium == "plasma" else "none"
        self.exponential = ExponentialTopside(f_n_sat_hz=float(f_n_sat_hz), scale_height_m=float(scale_height_m))
        self.geometry = geometry
        self.plate_model = ("3d" if geometry == "orbit" else "2d") if plate_model == "auto" else plate_model
        self.flat_radius_m = float(flat_radius_m)
        self.field_model = field_model
        self.v_s_mps = float(v_s_mps)
        self.ref_alt_m = None if ref_alt_m is None else float(ref_alt_m)
        self.earth_radius_m = float(earth_radius_m)
        self.doppler_dt_s = float(doppler_dt_s)
        self.group_df_fraction = float(group_df_fraction)
        self.layer_echo = bool(layer_echo)
        straight_profile = self.exponential if (medium == "plasma" and profile == "exponential") else None
        orbit = geometry == "orbit"
        specs: list[FacetSpec | SheetSpec] = []
        for m in facets or [DEFAULT_FACET]:
            kind = str(dict(m).get("kind", "plate"))
            if kind == "plate":
                specs.append(FacetSpec.from_mapping(m, straight_profile, orbit=orbit))
            elif kind == "sheet":
                if not orbit and medium == "plasma" and profile != "exponential":
                    raise ValueError("a straight-track sheet runs in the exponential plasma or free space")
                specs.append(SheetSpec.from_mapping(m, orbit=orbit))
            else:
                raise ValueError(f"facet kind must be 'plate' or 'sheet', not {kind!r}")
        self.facets = tuple(specs)
        if self.plate_model == "2d" and geometry == "orbit" and any(
                isinstance(f, FacetSpec) and f.halo is not None and f.halo.present for f in specs):
            raise ValueError("a halo on the orbit needs plate_model '3d'")
        #: The halo's reference field |E| of each facet at overflight, by (id, frequency).
        self._halo_ref: dict[tuple[str, float], float] = {}
        ids = [f.facet_id for f in self.facets]
        if len(set(ids)) != len(ids):
            raise ValueError("facet ids must be unique")
        #: The plate and sheet models are field free (ReportsGroupIndex hook).
        self.group_index_used: str | None = "field_free"
        #: Largest cross-track offset of the platform from a facet seen so far, m (orbit).
        self.max_cross_track_m = 0.0
        self._l0_cache: tuple[int, Any] | None = None
        self._azimuth: dict[str, float] = {}
        self._sheet_geo: dict[str, dict[str, Any]] = {}
        self._surfaces: dict[tuple[str, float], _Surface] = {}
        self._q: dict[tuple[str, float], np.ndarray] = {}

    def describe(self) -> Mapping[str, Any]:
        facets = []
        for f in self.facets:
            d = f.describe()
            if f.facet_id in self._azimuth:
                d["azimuth_deg_used"] = math.degrees(self._azimuth[f.facet_id])
            g = self._sheet_geo.get(f.facet_id)
            if g is not None:
                d.update({"dip_deg_used": g["dip_deg"], "declination_deg_used": g["declination_deg"],
                          "curvature_radius_m_used": g["curvature_radius_m"], "field_model": g["field_model"],
                          "field_provenance": g["provenance"]})
            facets.append(d)
        return {
            "name": self.name, "layer": self.layer, "medium": self.medium,
            "profile": self.profile_kind,
            "exponential": {"f_n_sat_hz": self.exponential.f_n_sat_hz,
                            "scale_height_m": self.exponential.scale_height_m},
            "geometry": self.geometry, "plate_model": self.plate_model, "flat_radius_m": self.flat_radius_m,
            "field_model": self.field_model, "v_s_mps": self.v_s_mps, "ref_alt_m": self.ref_alt_m,
            "earth_radius_m": self.earth_radius_m, "doppler_dt_s": self.doppler_dt_s,
            "group_df_fraction": self.group_df_fraction,
            "layer_echo": self.layer_echo,
            "facets": facets,
            "group_index_used": self.group_index_used,
        }

    # -- the medium ----------------------------------------------------------------

    def _base_profile(self, state: IonosphereState) -> Any:
        """The index against true depth below the reference level, or None in free space."""
        if self.medium == "free":
            return None
        if self.profile_kind == "exponential":
            return self.exponential
        key = id(state)
        if self._l0_cache is None or self._l0_cache[0] != key:
            self._l0_cache = (key, TabulatedTopside(state.alt_m, state.ne, float(self.ref_alt_m)))
        return self._l0_cache[1]

    # -- one row: Phase 2's 2-D plate ----------------------------------------------

    def _derive(self, e_field: np.ndarray, phase_path: np.ndarray, f: float,
                centre_rad: float = 0.0) -> np.ndarray:
        """Phase delay from the field and a reference phase path.

        The residual ``arg(E exp(+j 2 k0 P_c))`` is taken in the cycle centred
        on ``centre_rad``: 0 for a plate's centre ray, ``-pi/2`` for a sheet's
        shortest ray, which is the stationary-phase value of a two-dimensional
        minimum, so that the residual does not cross the cut during a pass and
        the phase delay does not jump by a cycle.
        """
        k0 = 2.0 * math.pi * f / C_LIGHT
        rel = e_field * np.exp(2j * k0 * phase_path) * np.exp(-1j * centre_rad)
        return 2.0 * (phase_path - (np.angle(rel) + centre_rad) / (2.0 * k0)) / C_LIGHT

    def _straight_row(self, facet: FacetSpec, f: float, t: float, base: Any):
        plate = facet.plate
        if base is not None:
            if float(base.mu_squared(f, 0.0)) <= 0.0:
                return None
            if float(base.mu_squared(f, plate.depth_m + 0.5 * plate.length_m)) <= 0.05:
                return None
        h = 0.5 * self.doppler_dt_s
        times = np.array([t, t - h, t + h])
        try:
            e = echo_track(times, f, plate, base, v_s_mps=self.v_s_mps)
        except ValueError:
            return None
        return float(e.group_path_m[0]), e.field, e.phase_path_m

    def _orbit_row(self, facet: FacetSpec, f: float, t: float, p: PlatformState, base: Any):
        r_ref = self.earth_radius_m + float(self.ref_alt_m)
        prof = FlattenedProfile(base if base is not None else _FreeSpace(), r_ref)
        r_f = self.earth_radius_m + float(facet.alt_m)
        r_p = self.earth_radius_m + float(p.alt_m)
        if float(prof.mu_squared(f, prof.flat_depth(r_p))) <= 0.0:
            return None
        src = facet.plate
        d_f = prof.flat_depth(r_f)
        half = 0.5 * src.length_m * r_ref / r_f
        if base is not None and float(prof.mu_squared(f, d_f + half)) <= 0.05:
            return None
        plate = Plate(x_m=0.0, depth_m=d_f, tilt_rad=src.tilt_rad, length_m=src.length_m * r_ref / r_f,
                      drift_mps=src.drift_mps * r_ref / r_p, taper=src.taper,
                      drift_x_mps=src.drift_x_mps * r_ref / r_f, tukey_fraction=src.tukey_fraction)
        if p.velocity_ned_mps is None:
            raise ValueError("the orbit geometry needs the platform velocity")
        at, xt = _track_angles(p, float(facet.lat_deg), float(facet.lon_deg))
        self.max_cross_track_m = max(self.max_cross_track_m, abs(xt) * r_f)
        v = np.asarray(p.velocity_ned_mps, dtype=float)
        v_h = float(math.hypot(v[0], v[1]))
        along0 = -r_ref * at
        extra0 = prof.flat_depth(r_p)
        d_along = r_ref * v_h / r_p
        d_extra = r_ref * float(v[2]) / r_p
        h = 0.5 * self.doppler_dt_s
        offs = np.array([0.0, -h, h])
        times = t + offs
        along = along0 + d_along * offs
        extra = extra0 + d_extra * offs
        try:
            e = echo_track(times, f, plate, prof, along_m=along, extra_depth_m=extra)
            df = self.group_df_fraction * f
            lo = echo_track(times[:1], f - df, plate, prof, along_m=along[:1], extra_depth_m=extra[:1])
            hi = echo_track(times[:1], f + df, plate, prof, along_m=along[:1], extra_depth_m=extra[:1])
        except ValueError:
            return None
        group = float(((f + df) * hi.phase_path_m[0] - (f - df) * lo.phase_path_m[0]) / (2.0 * df))
        return group, e.field, e.phase_path_m

    # -- one row: a straight-track sheet --------------------------------------------

    def _straight_sheet_row(self, spec: SheetSpec, f: float, t: float, base: Any):
        sheet = spec.sheet(math.radians(spec.dip_deg or 0.0), float(spec.curvature))
        if base is not None and float(base.mu_squared(f, 0.0)) <= 0.0:
            return None
        h = 0.5 * self.doppler_dt_s
        times = np.array([t, t - h, t + h])
        try:
            e = se.echo_track(times, f, sheet, base, height_m=float(spec.height_m),
                              heading_rad=math.radians(spec.heading_deg), v_s_mps=self.v_s_mps)
        except ValueError:
            return None
        return float(e.group_path_m[0]), e.field, e.phase_path_m

    # -- the 3-D surfaces ----------------------------------------------------------

    def sheet_geometry(self, facet_id: str, state: IonosphereState | None = None,
                       platform: PlatformState | None = None) -> dict[str, Any]:
        """The sheet's field direction and frame at its reference point, local (east, north, up).

        Dip and declination from IGRF (``field_model``) at the sheet's place and
        the epoch of ``state`` (or ``platform``), each unless given in the
        facet's mapping.  Latched on first use.
        """
        if facet_id in self._sheet_geo:
            return self._sheet_geo[facet_id]
        spec = next(f for f in self.facets if f.facet_id == facet_id)
        if not isinstance(spec, SheetSpec) or spec.lat_deg is None:
            raise ValueError(f"{facet_id!r} is not a sheet on an orbit")
        if spec.dip_deg is not None and spec.declination_deg is not None:
            dip, dec, model = spec.dip_deg, spec.declination_deg, "override"
            prov: dict[str, Any] = {"model": "override"}
        else:
            epoch = (getattr(state, "epoch", "") or "") or (getattr(platform, "epoch", "") or "")
            if not epoch and self.field_model != "dipole":
                raise ValueError("an orbit sheet needs an epoch for IGRF (or dip_deg and declination_deg)")
            from tss.l0_ionosphere.geomag import field_ned_t

            b, prov = field_ned_t(float(spec.lat_deg), float(spec.lon_deg), float(spec.alt_m), epoch,
                                  prefer=self.field_model)
            bn, be, bd = (float(v) for v in np.asarray(b)[:, 0])
            dip = math.degrees(math.atan2(bd, math.hypot(bn, be)))
            dec = math.degrees(math.atan2(be, bn))
            if spec.dip_deg is not None:
                dip = spec.dip_deg
            if spec.declination_deg is not None:
                dec = spec.declination_deg
            model = str(prov.get("model", ""))
        r_f = self.earth_radius_m + float(spec.alt_m)
        curv = (se.dipole_curvature_radius_m(math.radians(dip), r_f) if spec.curvature == "dipole"
                else float(spec.curvature))
        sheet = spec.sheet(math.radians(dip), curv)
        dr = math.radians(dec)
        # Columns: magnetic north, magnetic east, up, in local (east, north, up).
        rot = np.column_stack([[math.sin(dr), math.cos(dr), 0.0],
                               [math.cos(dr), -math.sin(dr), 0.0],
                               [0.0, 0.0, 1.0]])
        fr = se.sheet_frame(sheet)
        geo = {"dip_deg": float(dip), "declination_deg": float(dec), "curvature_radius_m": float(curv),
               "field_model": model, "provenance": dict(prov), "rotation": rot, "sheet": sheet,
               "b_local": rot @ fr["b"], "normal_local": rot @ fr["n"], "u_local": rot @ fr["u"]}
        self._sheet_geo[facet_id] = geo
        return geo

    def _surface(self, spec: FacetSpec | SheetSpec, f: float, state: IonosphereState,
                 platform: PlatformState) -> _Surface:
        key = (spec.facet_id, float(f))
        if key in self._surfaces:
            return self._surfaces[key]
        if isinstance(spec, SheetSpec):
            geo = self.sheet_geometry(spec.facet_id, state, platform)
            pts_m, w = se.sheet_points(geo["sheet"], f)
            rot = geo["rotation"]
            vel = rot @ se.sheet_velocity(geo["sheet"])
            surf = _Surface(pts=pts_m @ rot.T, w=w, vel=vel, rise_mps=float(vel[2]), ref="shortest")
        else:
            az = self._azimuth[spec.facet_id]
            width = spec.plate.length_m if spec.width_m is None else spec.width_m
            pts, w, _ = plate_surface(spec.plate, width_m=width, azimuth_rad=az,
                                      cross_tilt_rad=spec.cross_tilt_rad, f_hz=f)
            a_hat = np.array([math.sin(az), math.cos(az), 0.0])
            vel = spec.plate.drift_mps * np.array([0.0, 0.0, 1.0]) + spec.plate.drift_x_mps * a_hat
            surf = _Surface(pts=pts, w=w, vel=vel, rise_mps=float(spec.plate.drift_mps), ref="centre")
        self._surfaces[key] = surf
        return surf

    def _scatter_row(self, spec: FacetSpec | SheetSpec, f: float, t: float, p: PlatformState,
                     base: Any, state: IonosphereState):
        """One row from the 3-D ray engine (module docstring)."""
        h = 0.5 * self.doppler_dt_s
        offs = np.array([0.0, -h, h])
        if self.geometry == "orbit":
            frame = LocalFrame(float(spec.lat_deg), float(spec.lon_deg), self.earth_radius_m + float(spec.alt_m))
            if p.velocity_ned_mps is None:
                raise ValueError("the orbit geometry needs the platform velocity")
            pos = frame.from_llh(p.lat_deg, p.lon_deg, self.earth_radius_m + float(p.alt_m))
            vel_p = frame.vector_from_ned(np.asarray(p.velocity_ned_mps, dtype=float), p.lat_deg, p.lon_deg)
            r_base = float(self.ref_alt_m) - float(spec.alt_m)
            _, xt = _track_angles(p, float(spec.lat_deg), float(spec.lon_deg))
            self.max_cross_track_m = max(self.max_cross_track_m, abs(xt) * frame.r0_m)
            if isinstance(spec, FacetSpec) and spec.facet_id not in self._azimuth:
                self._azimuth[spec.facet_id] = (math.radians(spec.azimuth_deg) if spec.azimuth_deg is not None
                                                else math.atan2(float(vel_p[0]), float(vel_p[1])))
        else:
            if isinstance(spec, SheetSpec):
                raise ValueError("a straight-track sheet uses the flat sheet model")
            frame = LocalFrame(0.0, 0.0, self.flat_radius_m)
            plate = spec.plate
            self._azimuth.setdefault(spec.facet_id, 0.0)  # the track runs north
            pos = np.array([-spec.y_m, self.v_s_mps * t - plate.x_m, plate.depth_m])
            vel_p = np.array([0.0, self.v_s_mps, 0.0])
            r_base = float(plate.depth_m)  # the profile's depth is counted from the platform at t = 0
        surf = self._surface(spec, f, state, p)
        med0 = SphericalMedium(base, r_base)
        hp = float(pos[2]) + float(np.sum(pos[:2] ** 2)) / (2.0 * frame.r0_m)
        if base is not None:
            if float(med0.mu_squared(f, np.array([hp]))[0]) <= 0.0:
                return None
            z_low = float(np.min(surf.pts[:, 2])) - abs(surf.rise_mps) * (abs(t) + h)
            if float(med0.shifted(surf.rise_mps * t).mu_squared(f, np.array([z_low]))[0]) <= 0.05:
                return None
        k0 = 2.0 * math.pi * f / C_LIGHT
        fields = np.empty(3, dtype=complex)
        pcs = np.empty(3)
        group = math.nan
        qkey = (spec.facet_id, float(f))
        live = surf.w > 0
        for i, o in enumerate(offs):
            tau = t + float(o)
            pk = pos + vel_p * float(o)
            xk = surf.pts + surf.vel * tau
            med = med0.shifted(surf.rise_mps * tau)
            q0 = self._q.get(qkey)
            try:
                ph, gr, q = ray_paths(pk, xk, f, med, frame.r0_m, q0=q0,
                                      newton_iterations=8 if q0 is None else 4)
            except ValueError:
                return None
            self._q[qkey] = q
            fields[i] = np.sum(surf.w * np.exp(-2j * k0 * ph))
            j = 0 if surf.ref == "centre" else int(np.argmin(np.where(live, ph, np.inf)))
            pcs[i] = ph[j]
            if i == 0:
                group = float(gr[0]) if surf.ref == "centre" else float(np.min(gr[live]))
        return group, fields, pcs

    # -- the diffuse halo ------------------------------------------------------------

    def _halo_paths(self, spec: FacetSpec, f: float, times: np.ndarray, p: PlatformState, base: Any):
        """One-way phase and group paths to the scatterers, ``(times, points)``, and
        the facet's reference field at overflight.  ``None`` beyond reflection."""
        s, _ = spec.halo.draw(spec.plate.length_m)
        key = (spec.facet_id, float(f))
        plate = spec.plate
        if self.plate_model == "2d":  # the straight track
            try:
                phase, group = point_paths(times, f, plate, base, s, v_s_mps=self.v_s_mps)
                if key not in self._halo_ref:
                    t_over = plate.x_m / (self.v_s_mps - plate.drift_x_mps)
                    ref = echo_track(np.array([t_over]), f, plate, base, v_s_mps=self.v_s_mps)
                    self._halo_ref[key] = float(abs(ref.field[0]))
            except ValueError:
                return None
            return phase, group, self._halo_ref[key]
        # The 3-D engine: the same frame and medium as the facet's own row.
        surf = self._surfaces[(spec.facet_id, float(f))]
        az = self._azimuth[spec.facet_id]
        _, t_a = _plate_normal_and_along(plate, az, spec.cross_tilt_rad)
        pts = s[:, None] * t_a[None, :]
        if self.geometry == "orbit":
            frame = LocalFrame(float(spec.lat_deg), float(spec.lon_deg), self.earth_radius_m + float(spec.alt_m))
            vel_p = frame.vector_from_ned(np.asarray(p.velocity_ned_mps, dtype=float), p.lat_deg, p.lon_deg)
            r_base = float(self.ref_alt_m) - float(spec.alt_m)

            def platform(tau: float) -> np.ndarray:
                return frame.from_llh(p.lat_deg, p.lon_deg, self.earth_radius_m + float(p.alt_m)) + vel_p * (
                    tau - float(times[0]))

            overhead = np.array([0.0, 0.0, float(p.alt_m) - float(spec.alt_m)])
            t_over = float(times[0])
        else:
            frame = LocalFrame(0.0, 0.0, self.flat_radius_m)
            r_base = float(plate.depth_m)

            def platform(tau: float) -> np.ndarray:
                return np.array([-spec.y_m, self.v_s_mps * tau - plate.x_m, plate.depth_m])

            t_over = plate.x_m / self.v_s_mps
            overhead = platform(t_over)
        med0 = SphericalMedium(base, r_base)
        phase = np.empty((times.size, s.size))
        group = np.empty((times.size, s.size))
        try:
            for i, tau in enumerate(times):
                med = med0.shifted(surf.rise_mps * tau)
                phase[i], group[i], _ = ray_paths(platform(float(tau)), pts + surf.vel * tau, f, med, frame.r0_m)
            if key not in self._halo_ref:
                # On the orbit "overflight" is the platform overhead at its present height.
                tau_ref = 0.0 if self.geometry == "orbit" else t_over
                ph, _, _ = ray_paths(overhead, surf.pts + surf.vel * tau_ref, f,
                                     med0.shifted(surf.rise_mps * tau_ref), frame.r0_m)
                k0 = 2.0 * math.pi * f / C_LIGHT
                self._halo_ref[key] = float(abs(np.sum(surf.w * np.exp(-2j * k0 * ph))))
        except ValueError:
            return None
        return phase, group, self._halo_ref[key]

    def _halo_rows(self, spec: FacetSpec, f: float, t: float, p: PlatformState, base: Any) -> list[tuple]:
        """One row per scatterer of ``spec``'s halo (module docstring)."""
        tag = f"halo:{spec.facet_id}"
        h = 0.5 * self.doppler_dt_s
        times = np.array([t, t - h, t + h])
        got = self._halo_paths(spec, f, times, p, base)
        if got is None:
            return [(f, math.nan, math.nan, math.nan, math.nan, f"{tag};beyond_reflection")]
        phase, group, ref = got
        _, theta = spec.halo.draw(spec.plate.length_m)
        share = spec.halo.k if spec.halo.k > 0.0 else 1.0
        amplitude = ref / math.sqrt(share * spec.halo.n_points)
        lam = C_LIGHT / f
        rows = []
        for i in range(spec.halo.n_points):
            g = float(group[0, i])
            tau_p = 2.0 * float(phase[0, i]) / C_LIGHT - float(theta[i]) / (2.0 * math.pi * f)
            dopp = -f * 2.0 * float(phase[2, i] - phase[1, i]) / C_LIGHT / self.doppler_dt_s
            amp = (20.0 * math.log10(lam / (4.0 * math.pi * 2.0 * g)) + spec.reflection_db
                   + (20.0 * math.log10(amplitude) if amplitude > 0.0 else -math.inf))
            rows.append((f, 2.0 * g / C_LIGHT, tau_p, dopp, amp, tag))
        return rows

    # -- the layer echo ------------------------------------------------------------

    @staticmethod
    def layer_paths(state: IonosphereState, platform_alt_m: float, f_hz: float,
                    n_nodes: int = 96) -> tuple[float, float, float] | None:
        """(group path, phase path, mu at the platform) of the vertical O echo, m.

        Field free, from the platform down to where ``f_N = f``; ``None`` when
        the wave is evanescent at the platform or there is no reflection.
        """
        prof = TabulatedTopside(state.alt_m, state.ne, float(platform_alt_m))
        f = float(f_hz)
        if float(prof.mu_squared(f, 0.0)) <= 0.0:
            return None
        depth_max = float(platform_alt_m) - float(np.min(state.alt_m))
        grid = np.linspace(0.0, depth_max, 4001)
        m2 = prof.mu_squared(f, grid)
        below = np.nonzero(m2 <= 0.0)[0]
        if below.size == 0:
            return None
        j = int(below[0])
        lo, hi = float(grid[j - 1]), float(grid[j])
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if float(prof.mu_squared(f, mid)) > 0.0:
                lo = mid
            else:
                hi = mid
        d_r = 0.5 * (lo + hi)
        x, w = np.polynomial.legendre.leggauss(int(n_nodes))
        s_max = math.sqrt(d_r)
        s = 0.5 * s_max * (x + 1.0)
        ws = 0.5 * s_max * w
        depth = d_r - s * s
        mu = np.sqrt(np.clip(prof.mu_squared(f, depth), 1e-30, None))
        group = float(np.sum(ws * 2.0 * s / mu))
        phase = float(np.sum(ws * 2.0 * s * mu))
        return group, phase, float(math.sqrt(float(prof.mu_squared(f, 0.0))))

    def _layer_row(self, state: IonosphereState, p: PlatformState, f: float):
        got = self.layer_paths(state, p.alt_m, f)
        if got is None:
            return (f, math.nan, math.nan, math.nan, math.nan, "layer:O;no_echo")
        group, phase, mu_sat = got
        v_down = 0.0 if p.velocity_ned_mps is None else float(np.asarray(p.velocity_ned_mps)[2])
        dopp = 2.0 * f * mu_sat * v_down / C_LIGHT
        lam = C_LIGHT / f
        amp = 20.0 * math.log10(lam / (4.0 * math.pi * 2.0 * group))
        return (f, 2.0 * group / C_LIGHT, 2.0 * phase / C_LIGHT, dopp, amp, "layer:O")

    # -- the contract --------------------------------------------------------------

    def _row(self, spec: FacetSpec | SheetSpec, f: float, t: float, platform: PlatformState, base: Any,
             state: IonosphereState):
        if isinstance(spec, SheetSpec):
            if self.geometry == "straight":
                return self._straight_sheet_row(spec, f, t, base)
            return self._scatter_row(spec, f, t, platform, base, state)
        if self.plate_model == "2d":
            if self.geometry == "straight":
                return self._straight_row(spec, f, t, base)
            return self._orbit_row(spec, f, t, platform, base)
        return self._scatter_row(spec, f, t, platform, base, state)

    def echoes(
        self,
        state: IonosphereState,
        platform: PlatformState,
        freq_hz: np.ndarray,
        rng: np.random.Generator,
    ) -> EchoSet:
        """One row per facet and frequency at the platform's ``t_s`` (0 if unset)."""
        t = 0.0 if getattr(platform, "t_s", None) is None else float(platform.t_s)
        freqs = np.atleast_1d(np.asarray(freq_hz, dtype=float))
        base = self._base_profile(state)
        rows: list[tuple[float, float, float, float, float, str]] = []
        for facet in self.facets:
            tag = f"facet:{facet.facet_id}" + (";sheet" if isinstance(facet, SheetSpec) else "")
            for f in freqs:
                f = float(f)
                got = self._row(facet, f, t, platform, base, state)
                halo = facet.halo if isinstance(facet, FacetSpec) else None
                if got is None:
                    rows.append((f, math.nan, math.nan, math.nan, math.nan, f"{tag};beyond_reflection"))
                    continue
                diffuse = (self._halo_rows(facet, f, t, platform, base)
                           if halo is not None and halo.present else [])
                if halo is not None and halo.k == 0.0:  # the facet is absent; only its halo is received
                    rows.extend(diffuse)
                    continue
                group, fields, pcs = got
                centre = -0.5 * math.pi if isinstance(facet, SheetSpec) else 0.0
                tau_p = float(self._derive(np.asarray(fields[:1]), np.asarray(pcs[:1]), f, centre)[0])
                d_tau = pair_phase_delay_difference(fields[1], fields[2], pcs[1], pcs[2], f)
                dopp = -f * d_tau / self.doppler_dt_s
                lam = C_LIGHT / f
                mag = abs(fields[0])
                amp = (20.0 * math.log10(lam / (4.0 * math.pi * 2.0 * group)) + facet.reflection_db
                       + (20.0 * math.log10(mag) if mag > 0.0 else -math.inf))
                rows.append((f, 2.0 * group / C_LIGHT, tau_p, dopp, amp, tag))
                rows.extend(diffuse)
        if self.layer_echo:
            for f in freqs:
                rows.append(self._layer_row(state, platform, float(f)))
        if rows:
            f_a, tg, tp, fd, am, fl = (np.array(c) for c in zip(*rows, strict=True))
        else:
            f_a = tg = tp = fd = am = np.array([], dtype=float)
            fl = np.array([], dtype=object)
        n = f_a.size
        return EchoSet(
            freq_hz=f_a.astype(float), mode=np.array(["O"] * n, dtype=object),
            group_delay_s=tg.astype(float), amplitude_db=am.astype(float),
            absorption_db=np.zeros(n), doppler_hz=fd.astype(float),
            flags=np.array(fl, dtype=object), phase_delay_s=tp.astype(float),
            source=(f"l1/facet ({self.medium}, {self.profile_kind}, {self.geometry}, {self.plate_model}, "
                    f"t={t:.6g} s)"),
        )
