"""Layer: L1
WP: WP5
Claims: C3
Status: partial
Requirements: SCI-15, INS-09, INS-15

The facet tier (``l1/facet``): the verified physical-optics plate model as an
L1 tier, one echo set per sweep (WP5 Phase 2, decision D4).

WHAT IT DOES
------------
For each sweep of a drift-mode sequence the pipeline calls :meth:`FacetTier.echoes`
with the platform state of that sweep, whose ``t_s`` is the time from the start
of the sequence.  For every facet and every sounding frequency the tier calls
:func:`tss.l1_propagation.plate_echo.echo_track` and returns one EchoSet row:

* ``group_delay_s``: twice the group path to the plate's centre over ``c``;
* ``phase_delay_s``: the phase delay of the summed field, ``2 P_eff / c`` with
  ``P_eff = P_c - arg(E exp(+j 2 k0 P_c)) / (2 k0)``, ``P_c`` the centre ray's
  phase path; ``2 pi f tau_p`` is then ``-arg(E)`` to rounding, in L2's sign;
* ``doppler_hz``: ``-f d tau_p / dt`` by a central difference over
  ``doppler_dt_s`` (closing positive, as :mod:`tss.l1_propagation.doppler`);
* ``amplitude_db``: the free-space spreading of the analytic tier over twice
  the group path, plus ``reflection_db``, plus ``20 log10 |E|``.  ``|E|`` is 1
  at the specular peak of a uniform plate, so ``reflection_db`` sets the echo
  strength there;
* ``flags``: ``facet:<id>``, the name L2 keys fading and sequence state on.

A facet the frequency cannot reach (below the platform's plasma frequency, or
at or past the reflection level) is a masked row, ``facet:<id>;beyond_reflection``,
with NaN delays (contract rule 3).  The plate model is called unchanged (D4).

STRAIGHT GEOMETRY (verified)
----------------------------
``geometry="straight"`` is the plate model's own straight level track at
``v_s_mps``; the platform state supplies only the time, and a facet is placed
by ``depth_m`` (or ``mu_facet`` at ``mu_at_hz``) and ``x_m`` relative to the
platform at ``t = 0``.  ``medium="plasma"`` with ``profile="exponential"`` is the
reference topside of Technical Volume Chapter 8 (``f_n_sat_hz`` 1.5 MHz,
``scale_height_m`` 200 km), the c3_verify profile.  This is the configuration
the Phase 2 exit criterion was met on.

ORBIT GEOMETRY AND THE L0 PROFILE (new, unverified; Andrew, 2026-11-25)
-----------------------------------------------------------------------
``geometry="orbit"`` takes the platform's position and velocity from each
platform state (an orbit, :func:`tss.platform.orbits.platform_state_at`) and
places each facet by its sub-point (``lat_deg``, ``lon_deg``) and altitude
(``alt_m``).  The plate model is flat, so the sphere enters through the
Earth-flattening map about the reference radius ``R = R_E + ref_alt_m``:

    x = R psi,    d = R ln(R / r),    mu'^2(d) = mu^2(r) (r / R)^2,

with ``psi`` the along-track angle and ``r`` the radius.  The map is conformal,
so tilts are kept, and Snell's law in the flat frame is Bouguer's law
``r mu sin(theta) = const`` on the sphere, so every ray and its phase path
``int mu ds`` are kept exactly (Kerr 1951).  The plate model's group-path
integral assumes ``mu_g = 1/mu``, which the flattened index does not obey, so
in this geometry the group path is taken as ``d(f P)/df`` by a central
difference in frequency.  What the sphere does to route one is real: near
overflight the range's second derivative is ``v^2 (r_f / r_p) / D``, not the
flat ``v^2 / D``, about 7 % less at 300 km below a 600 km orbit
(``tests/unit/test_l1_facet_orbit.py``), so an estimator that assumes the flat
slope reads drift about 7 % wrong.  The plate is two-dimensional: the facet is
taken in the platform's track plane and the platform's cross-track offset from
it is ignored; the largest seen is in :attr:`FacetTier.max_cross_track_m`.
Drift is the medium rising past the platform (the plate model's frame), so it
is mapped with the scale at the platform's height; along-track drift is mapped
with the scale at the facet.

``profile="l0"`` takes the electron density from the L0 state the pipeline
passes (``f_N^2 = 80.616 N_e``, interpolated in ``log N_e``), with
``ref_alt_m`` as the level depths are counted from.  The 64-node quadrature of
the plate model assumes a smooth integrand; on a tabulated profile it is held
to the exponential on the same grid by the second difference of the phase
across sweeps (same test file).

THE LAYER ECHO
--------------
With ``layer_echo=True`` each frequency also gets the vertical echo of the
layer below, from the L0 state the pipeline passes, flagged ``layer:O``: the
field-free group path ``int dz / mu`` and phase path ``int mu dz`` from the
platform to the reflection level (by the substitution ``z = z_r - s^2``, which
takes the square-root singularity out), its spreading over twice the group
path, and the Doppler of the platform's vertical velocity.  No absorption, no
magnetoionic split, no layer drift.  A frequency below the platform's plasma
frequency or above the layer's critical frequency gives a masked row.

WHAT IS LEFT OUT
----------------
Everything the plate model leaves out (the magnetic field, absorption,
reflectivity variation across the plate); the sheet model (``sheet_echo``) is
not wrapped yet; there is no diffuse halo (Phase 3).
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from tss.constants import C_LIGHT
from tss.interfaces import EchoSet, IonosphereState, PlatformState
from tss.l1_propagation.plate_echo import ExponentialTopside, Plate, echo_track
from tss.platform.orbits import RE
from tss.registry import register

__all__ = [
    "FacetSpec",
    "FacetTier",
    "DEFAULT_FACET",
    "TabulatedTopside",
    "FlattenedProfile",
    "FN2_PER_NE",
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
    """Earth-flattened index ``mu'^2(d) = mu^2(r) (r/R)^2`` (module docstring).

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


@dataclass(frozen=True)
class FacetSpec:
    """One facet: its id, its plate, its reflection (dB) and, on an orbit, its place."""

    facet_id: str
    plate: Plate
    reflection_db: float = 0.0
    lat_deg: float | None = None
    lon_deg: float | None = None
    alt_m: float | None = None

    @classmethod
    def from_mapping(cls, m: Mapping[str, Any], profile: Any, *, orbit: bool) -> FacetSpec:
        """From a configuration mapping.

        Keys: ``id``; on the straight track ``depth_m``, or ``mu_facet`` with
        ``mu_at_hz`` (inside the plasma), and ``x_m``; on an orbit ``lat_deg``,
        ``lon_deg`` and ``alt_m``; and ``tilt_deg``, ``length_m``,
        ``drift_mps``, ``drift_x_mps``, ``taper``, ``reflection_db``.
        """
        d = dict(m)
        lat = lon = alt = None
        if orbit:
            if not all(k in d for k in ("lat_deg", "lon_deg", "alt_m")):
                raise ValueError("on an orbit a facet needs lat_deg, lon_deg and alt_m")
            lat, lon, alt = float(d["lat_deg"]), float(d["lon_deg"]), float(d["alt_m"])
            depth = 0.0  # set by the tier in the flattened frame
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
        return cls(facet_id=str(d.get("id", "f1")), plate=plate,
                   reflection_db=float(d.get("reflection_db", 0.0)), lat_deg=lat, lon_deg=lon, alt_m=alt)

    def describe(self) -> dict[str, Any]:
        p = self.plate
        return {"id": self.facet_id, "x_m": p.x_m, "depth_m": p.depth_m,
                "tilt_deg": math.degrees(p.tilt_rad), "length_m": p.length_m,
                "drift_mps": p.drift_mps, "drift_x_mps": p.drift_x_mps, "taper": p.taper,
                "reflection_db": self.reflection_db, "lat_deg": self.lat_deg,
                "lon_deg": self.lon_deg, "alt_m": self.alt_m}


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


@register("l1", "facet")
class FacetTier:
    """Facet echoes per sweep from the plate model.  See the module docstring."""

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
    ) -> None:
        if medium not in ("plasma", "free"):
            raise ValueError("medium must be 'plasma' or 'free'")
        if profile not in ("exponential", "l0"):
            raise ValueError("profile must be 'exponential' or 'l0'")
        if geometry not in ("straight", "orbit"):
            raise ValueError("geometry must be 'straight' or 'orbit'")
        if not float(doppler_dt_s) > 0.0:
            raise ValueError("doppler_dt_s must be positive")
        if (geometry == "orbit" or (medium == "plasma" and profile == "l0")) and ref_alt_m is None:
            raise ValueError("the orbit geometry and the l0 profile need ref_alt_m")
        self.medium = medium
        self.profile_kind = profile if medium == "plasma" else "none"
        self.exponential = ExponentialTopside(f_n_sat_hz=float(f_n_sat_hz), scale_height_m=float(scale_height_m))
        self.geometry = geometry
        self.v_s_mps = float(v_s_mps)
        self.ref_alt_m = None if ref_alt_m is None else float(ref_alt_m)
        self.earth_radius_m = float(earth_radius_m)
        self.doppler_dt_s = float(doppler_dt_s)
        self.group_df_fraction = float(group_df_fraction)
        self.layer_echo = bool(layer_echo)
        straight_profile = self.exponential if (medium == "plasma" and profile == "exponential") else None
        self.facets = tuple(
            FacetSpec.from_mapping(m, straight_profile, orbit=(geometry == "orbit"))
            for m in (facets or [DEFAULT_FACET])
        )
        ids = [f.facet_id for f in self.facets]
        if len(set(ids)) != len(ids):
            raise ValueError("facet ids must be unique")
        #: The plate model is field free (ReportsGroupIndex hook).
        self.group_index_used: str | None = "field_free"
        #: Largest cross-track offset of the platform from a facet seen so far, m
        #: (orbit geometry; the plate model is two-dimensional).
        self.max_cross_track_m = 0.0
        self._l0_cache: tuple[int, Any] | None = None

    def describe(self) -> Mapping[str, Any]:
        return {
            "name": self.name, "layer": self.layer, "medium": self.medium,
            "profile": self.profile_kind,
            "exponential": {"f_n_sat_hz": self.exponential.f_n_sat_hz,
                            "scale_height_m": self.exponential.scale_height_m},
            "geometry": self.geometry, "v_s_mps": self.v_s_mps, "ref_alt_m": self.ref_alt_m,
            "earth_radius_m": self.earth_radius_m, "doppler_dt_s": self.doppler_dt_s,
            "group_df_fraction": self.group_df_fraction,
            "layer_echo": self.layer_echo,
            "facets": [f.describe() for f in self.facets],
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

    # -- one row -----------------------------------------------------------------

    def _derive(self, e_field: np.ndarray, phase_path: np.ndarray, f: float) -> np.ndarray:
        k0 = 2.0 * math.pi * f / C_LIGHT
        rel = e_field * np.exp(2j * k0 * phase_path)
        return 2.0 * (phase_path - np.angle(rel) / (2.0 * k0)) / C_LIGHT

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
        tau_p = self._derive(e.field, e.phase_path_m, f)
        return float(e.group_path_m[0]), tau_p, complex(e.field[0])

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
        tau_p = self._derive(e.field, e.phase_path_m, f)
        return group, tau_p, complex(e.field[0])

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
            for f in freqs:
                f = float(f)
                got = (self._straight_row(facet, f, t, base) if self.geometry == "straight"
                       else self._orbit_row(facet, f, t, platform, base))
                if got is None:
                    rows.append((f, math.nan, math.nan, math.nan, math.nan,
                                 f"facet:{facet.facet_id};beyond_reflection"))
                    continue
                group, tau_p, field0 = got
                dopp = -f * float(tau_p[2] - tau_p[1]) / self.doppler_dt_s
                lam = C_LIGHT / f
                mag = abs(field0)
                amp = (20.0 * math.log10(lam / (4.0 * math.pi * 2.0 * group)) + facet.reflection_db
                       + (20.0 * math.log10(mag) if mag > 0.0 else -math.inf))
                rows.append((f, 2.0 * group / C_LIGHT, float(tau_p[0]), dopp, amp, f"facet:{facet.facet_id}"))
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
            source=f"l1/facet ({self.medium}, {self.profile_kind}, {self.geometry}, t={t:.6g} s)",
        )

