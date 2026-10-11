"""Layer: study
WP: WP5, WP7
Claims: C3
Status: implemented
Requirements: SCI-15, SCI-23, SCI-27

Route two, phase 2: drift and tilt over an area from the oblique links of a
networked constellation, at one instant, for nine sensors and for 128.

What it asks
------------
The two-node study (:mod:`tss.studies.route_two`) measured drift and tilt by
looking at the same place twice.  A constellation whose sounders hear each
other does it differently: every pair of satellites in range links obliquely,
at every frequency that reflects, and each link's reflection point lands
somewhere in the ionosphere between them.  A patch of ionosphere hit by links
of different obliquity and different motion is measured completely, all at
once, with no wait for a second pass (:mod:`tss.l3_signal.novel.network_tilt_estimator`).
This study asks how much of the ionosphere that covers, how well, and how the
answer changes from the nine-satellite cluster of Chartier et al. to a
128-satellite constellation.

Method
------
Each scenario is a constellation, a velocity assignment, an ionosphere and a
measurement-error model.  For each snapshot every pair of satellites in range
and every satellite on its own (the vertical sounding) is evaluated at 2-20
MHz in 1 MHz steps with the tabulated stratified-ray model of
:mod:`tss.l1_propagation.link_table`.  A link counts when its reflected ray
exists, it lies at least 1% of the branch inside both its skip distance and its
maximum range, and its signal-to-noise ratio, scaled from 10 dB at 1000 km of
group path as the inverse square of the group path, is at least 0 dB.  Its
group path and line-of-sight rate errors scale with the same amplitude:

    sigma = sigma_ref * max(P' / 1000 km, 1).

Each link's reflection point falls in one patch of an equal-area grid (2
degrees of latitude, and the same width in kilometres in longitude, about 220
km square).  The Fisher information of each patch in ``(w, gE, gN, s, scale,
top_scale)`` is the sum over its links.  The patch is *measured* when the
information is full rank, and meets the goal when the one-sigma drift is
within 25 m/s (the SCI-23 goal scale) and each tilt component within 0.29
degrees.  The threshold variant uses 50 m/s.

Three versions of every patch are scored:

* all links, with the profile scale and topside scale height solved for;
* all links, with the profile taken as known;
* the vertical soundings alone, which is route one.

The ray model, the shifted-sphere forward model and the joint profile solve
were verified against the Haselgrove ray tracer of the Chartier reproduction
before this study was run (``tests/unit/test_route_two_network.py`` keeps the
table against the ray-by-ray model; the ray-tracer check is recorded in the
project notes because it needs ``nisc-verify``).

Not modelled: the Earth's magnetic field, absorption, focusing and defocusing
near the caustics (excluded by the margin rather than modelled), horizontal
gradients other than the patch's own tilt, fading, multipath separation of the
reflected echo from the direct line-of-sight signal, clock offsets between
satellites (a reciprocal pair of links cancels them; only one direction of
each pair is counted) and any information from revisiting a patch later.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

R_EARTH_M = 6371.0e3
MU = 3.986004418e14
OMEGA_E = 7.2921150e-5
FREQS_HZ = tuple(float(f) * 1e6 for f in range(2, 21))
P_REF_M = 1000e3
SNR_REF_DB = 10.0
SNR_MIN_DB = 0.0
SIGMA_RATE_REF = 5.0          # m/s at 1000 km, the route-one budget without fading
SIGMA_GROUP_REF = 100.0       # m at 1000 km, the two-node study's range error
GOAL_DRIFT = 25.0
THRESH_DRIFT = 50.0
GOAL_TILT_DEG = 0.29
PATCH_DEG = 2.0
#: Priors for the drift-only score: tilt of a degree per component, and the
#: profile errors the verification showed the joint solve absorbing (3% in
#: density, 5% in topside scale height), as from the nodes' own inversions.
PRIOR_TILT_RAD = math.radians(1.0)
PRIOR_SCALE = 0.03
PRIOR_TOP = 0.05
SNAPSHOT_STEP_S = 60.0
N_SNAPSHOTS = 101             # 100 min, about one orbit


# ---------------------------------------------------------------------------
# Constellations
# ---------------------------------------------------------------------------

def _enu(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    e = np.array([-np.sin(lo), np.cos(lo), 0.0])
    n = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    u = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    return e, n, u


def grid_cluster(lats, lons, alt_m, headings) -> tuple[np.ndarray, np.ndarray]:
    """Fixed grid of satellites on circular orbits; ``headings`` maps (lat, lon) -> (inclination, ascending).

    Velocities are Earth-fixed: the inertial circular velocity along the
    orbit's heading at that latitude, less the Earth's rotation.
    """
    r = R_EARTH_M + alt_m
    v = math.sqrt(MU / r)
    pos, vel = [], []
    for la in lats:
        for lo in lons:
            inc, asc = headings(la, lo)
            e, n, u = _enu(la, lo)
            s = math.cos(math.radians(inc)) / math.cos(math.radians(la))
            s = max(-1.0, min(1.0, s))
            c = math.sqrt(1.0 - s * s) * (1.0 if asc else -1.0)
            vi = v * (s * e + c * n)
            p = r * u
            pos.append(p)
            vel.append(vi - np.cross([0.0, 0.0, OMEGA_E], p))
    return np.array(pos), np.array(vel)


def walker(n_planes: int, per_plane: int, inc_deg: float, alt_m: float, raan_span_deg: float, phase_f: int,
           t_s: float) -> tuple[np.ndarray, np.ndarray]:
    """Walker constellation, circular, spherical Earth, Earth-fixed states at ``t_s``."""
    r = R_EARTH_M + alt_m
    nmo = math.sqrt(MU / r**3)
    inc = math.radians(inc_deg)
    total = n_planes * per_plane
    theta = OMEGA_E * t_s
    rot = np.array([[math.cos(theta), math.sin(theta), 0.0], [-math.sin(theta), math.cos(theta), 0.0],
                    [0.0, 0.0, 1.0]])
    pos, vel = [], []
    for p in range(n_planes):
        node = math.radians(p * raan_span_deg / n_planes)
        for k in range(per_plane):
            u = math.radians(k * 360.0 / per_plane + p * phase_f * 360.0 / total) + nmo * t_s
            cu, su, co, so = math.cos(u), math.sin(u), math.cos(node), math.sin(node)
            ri = r * np.array([co * cu - so * su * math.cos(inc), so * cu + co * su * math.cos(inc), su * math.sin(inc)])
            vi = r * nmo * np.array([-co * su - so * cu * math.cos(inc), -so * su + co * cu * math.cos(inc),
                                     cu * math.sin(inc)])
            re = rot @ ri
            pos.append(re)
            vel.append(rot @ vi - np.cross([0.0, 0.0, OMEGA_E], re))
    return np.array(pos), np.array(vel)


CHARTIER_LATS = (10.0, 20.0, 30.0)
CHARTIER_LONS = (-150.0, -140.0, -130.0)


def _cluster_headings(variant: str):
    if variant == "comoving_polar":
        return lambda la, lo: (86.4, True)
    if variant == "mixed_polar":
        return lambda la, lo: (86.4, lo != -140.0)
    if variant == "mixed_53":
        return lambda la, lo: (53.0, lo != -140.0)
    raise ValueError(variant)


@dataclass(frozen=True)
class Scenario:
    name: str
    kind: str                 # "cluster9", "walker"
    alt_m: float
    profile: str = "day"
    variant: str = ""
    walker: tuple = ()        # (planes, per_plane, inc, raan_span, F)
    noise_factor: float = 1.0
    snapshots: int = 1

    def states(self, t_s: float):
        if self.kind == "cluster9":
            return grid_cluster(CHARTIER_LATS, CHARTIER_LONS, self.alt_m, _cluster_headings(self.variant))
        p, k, inc, span, f = self.walker
        return walker(p, k, inc, self.alt_m, span, f, t_s)


SCENARIOS = (
    Scenario("cluster9_comoving_polar", "cluster9", 890e3, variant="comoving_polar"),
    Scenario("cluster9_mixed_polar", "cluster9", 890e3, variant="mixed_polar"),
    Scenario("cluster9_mixed_53", "cluster9", 890e3, variant="mixed_53"),
    Scenario("walker9_star", "walker", 1000e3, walker=(3, 3, 86.4, 180.0, 1), snapshots=N_SNAPSHOTS),
    Scenario("star128", "walker", 1000e3, walker=(8, 16, 86.4, 180.0, 4), snapshots=N_SNAPSHOTS),
    Scenario("delta128_53", "walker", 1000e3, walker=(8, 16, 53.0, 360.0, 1), snapshots=N_SNAPSHOTS),
    Scenario("cluster9_mixed_polar_night", "cluster9", 890e3, profile="night", variant="mixed_polar"),
    Scenario("star128_night", "walker", 1000e3, profile="night", walker=(8, 16, 86.4, 180.0, 4),
             snapshots=N_SNAPSHOTS),
    Scenario("star128_diurnal", "walker", 1000e3, profile="diurnal", walker=(8, 16, 86.4, 180.0, 4),
             snapshots=N_SNAPSHOTS),
    Scenario("delta128_53_diurnal", "walker", 1000e3, profile="diurnal", walker=(8, 16, 53.0, 360.0, 1),
             snapshots=N_SNAPSHOTS),
    Scenario("star64", "walker", 1000e3, walker=(8, 8, 86.4, 180.0, 4), snapshots=N_SNAPSHOTS),
    Scenario("star256", "walker", 1000e3, walker=(16, 16, 86.4, 180.0, 8), snapshots=N_SNAPSHOTS),
    Scenario("star512", "walker", 1000e3, walker=(16, 32, 86.4, 180.0, 8), snapshots=N_SNAPSHOTS),
    Scenario("cluster9_mixed_polar_noise2", "cluster9", 890e3, variant="mixed_polar", noise_factor=2.0),
    Scenario("star128_noise2", "walker", 1000e3, walker=(8, 16, 86.4, 180.0, 4), noise_factor=2.0,
             snapshots=N_SNAPSHOTS),
)


def profiles():
    from tss.l1_propagation.stratified_link import QUIET_DAYTIME, ChapmanProfile
    night = ChapmanProfile(layers=((2.5e11, 350e3, 50e3, 110e3, 0.5), (4.0e9, 110e3, 10e3, 10e3, 0.5)))
    return {"day": QUIET_DAYTIME, "night": night}


# ---------------------------------------------------------------------------
# Patches
# ---------------------------------------------------------------------------

class PatchGrid:
    """Equal-area grid: bands of ``dlat`` degrees, each split into cells about as wide as they are tall."""

    def __init__(self, dlat: float = PATCH_DEG):
        self.dlat = dlat
        self.nb = int(round(180.0 / dlat))
        lat_c = -90.0 + dlat * (np.arange(self.nb) + 0.5)
        self.nlon = np.maximum(1, np.round(360.0 * np.cos(np.radians(lat_c)) / dlat)).astype(int)
        self.offset = np.concatenate([[0], np.cumsum(self.nlon)[:-1]])
        self.n = int(self.nlon.sum())
        self.lat_c = lat_c

    def index(self, lat, lon):
        b = np.clip(np.floor((np.asarray(lat) + 90.0) / self.dlat).astype(int), 0, self.nb - 1)
        w = 360.0 / self.nlon[b]
        j = np.floor(((np.asarray(lon) + 180.0) % 360.0) / w).astype(int)
        return self.offset[b] + np.minimum(j, self.nlon[b] - 1)

    def centre(self, pid):
        pid = np.asarray(pid)
        b = np.searchsorted(self.offset, pid, side="right") - 1
        j = pid - self.offset[b]
        w = 360.0 / self.nlon[b]
        return self.lat_c[b], -180.0 + (j + 0.5) * w

    def area_km2(self, pid):
        pid = np.asarray(pid)
        b = np.searchsorted(self.offset, pid, side="right") - 1
        la1 = np.radians(-90.0 + self.dlat * b)
        la2 = la1 + np.radians(self.dlat)
        return (R_EARTH_M / 1e3) ** 2 * (2.0 * np.pi / self.nlon[b]) * (np.sin(la2) - np.sin(la1))


# ---------------------------------------------------------------------------
# One snapshot
# ---------------------------------------------------------------------------

def _llh(v):
    r = np.linalg.norm(v, axis=1)
    return np.degrees(np.arcsin(v[:, 2] / r)), np.degrees(np.arctan2(v[:, 1], v[:, 0]))


def _merge(ev_day: dict, ev_night: dict, day: np.ndarray) -> dict:
    out = {}
    for k, v in ev_day.items():
        if isinstance(v, np.ndarray) and v.shape[:1] == day.shape:
            sel = day.reshape((-1,) + (1,) * (v.ndim - 1))
            out[k] = np.where(sel, v, ev_night[k])
        else:
            out[k] = v
    return out


def local_time_h(lon_deg, t_s: float):
    """Local solar time at an equinox, with the Sun over longitude 0 at t = 0."""
    return (12.0 + t_s / 3600.0 + np.asarray(lon_deg) / 15.0) % 24.0


def snapshot(table, pos, vel, grid: PatchGrid, noise_factor: float = 1.0, *, night_table=None,
             t_s: float = 0.0) -> dict[str, Any]:
    """Per-patch Fisher information from all the links of one snapshot.

    With ``night_table`` the ionosphere is diurnal: a link uses ``table`` when
    its reflection point is between 06 and 18 local time and ``night_table``
    otherwise.
    """
    from tss.l1_propagation.link_table import rows_for_links

    n = len(pos)
    ii, jj = np.triu_indices(n)                     # includes i == j, the vertical soundings
    uh = pos / np.linalg.norm(pos, axis=1, keepdims=True)
    ang = np.arccos(np.clip(np.sum(uh[ii] * uh[jj], axis=1), -1.0, 1.0))
    theta_max = max(table.theta_max, night_table.theta_max if night_table is not None else 0.0)
    keep = ang <= theta_max + 1e-9
    ii, jj = ii[keep], jj[keep]
    vertical = ii == jj
    a, b, va, vb = pos[ii], pos[jj], vel[ii], vel[jj]
    fim = {}
    stats = {"pairs_in_range": int(np.sum(~vertical)), "links_oblique": 0, "links_vertical": 0,
             "per_freq_oblique": {}, "per_freq_vertical": {}}
    freqs = sorted(set(table.tables) | (set(night_table.tables) if night_table is not None else set()))
    for f in freqs:
        ev = table.evaluate(f, a, b, va, vb)
        if night_table is not None:
            evn = night_table.evaluate(f, a, b, va, vb)
            if "mid" not in ev and "mid" not in evn:
                continue
            ref = ev if "mid" in ev else evn
            _, lon_m = _llh(ref["mid"])
            day = (local_time_h(lon_m, t_s) >= 6.0) & (local_time_h(lon_m, t_s) < 18.0)
            if "mid" not in ev:
                ev = {**evn, "ok": evn["ok"] & ~day}
            elif "mid" not in evn:
                ev = {**ev, "ok": ev["ok"] & day}
            else:
                ev = _merge(ev, evn, day)
        ok = ev["ok"]
        if not np.any(ok):
            continue
        snr = SNR_REF_DB - 20.0 * np.log10(np.where(ok, ev["group"], np.inf) / P_REF_M)
        ok = ok & (snr >= SNR_MIN_DB)
        idx = np.flatnonzero(ok)
        if idx.size == 0:
            continue
        stats["links_oblique"] += int(np.sum(~vertical[idx]))
        stats["links_vertical"] += int(np.sum(vertical[idx]))
        stats["per_freq_oblique"][f / 1e6] = int(np.sum(~vertical[idx]))
        stats["per_freq_vertical"][f / 1e6] = int(np.sum(vertical[idx]))
        lat, lon = _llh(ev["mid"][idx])
        pid = grid.index(lat, lon)
        clat, clon = grid.centre(pid)
        grp, rate = rows_for_links(ev, idx, clat, clon)
        amp = np.maximum(ev["group"][idx] / P_REF_M, 1.0) * noise_factor
        grp = grp / (SIGMA_GROUP_REF * amp)[:, None]
        rate = rate / (SIGMA_RATE_REF * amp)[:, None]
        outer = grp[:, :, None] * grp[:, None, :] + rate[:, :, None] * rate[:, None, :]
        for k, p in enumerate(pid):
            ent = fim.setdefault(int(p), {"all": np.zeros((6, 6)), "vert": np.zeros((6, 6)), "n_obl": 0, "n_vert": 0})
            ent["all"] += outer[k]
            if vertical[idx[k]]:
                ent["vert"] += outer[k]
                ent["n_vert"] += 1
            else:
                ent["n_obl"] += 1
    return {"fim": fim, "stats": stats}


def _sigmas(f: np.ndarray, rel: float = 1e-9) -> np.ndarray:
    """(sigma_w, sigma_gE deg, sigma_gN deg), inf where the information is rank deficient.

    The parameters have very different units, so the test for rank is made on
    the information matrix scaled to unit diagonal.
    """
    d = np.sqrt(np.diag(f))
    if np.any(d <= 0.0):
        return np.array([np.inf, np.inf, np.inf])
    fs = f / np.outer(d, d)
    w, v = np.linalg.eigh(fs)
    if w[0] <= rel * w[-1]:
        return np.array([np.inf, np.inf, np.inf])
    cov = ((v / w) @ v.T) / np.outer(d, d)
    s = np.sqrt(np.diag(cov)[:3])
    return np.array([s[0], np.degrees(s[1]), np.degrees(s[2])])


PRIOR_INFO = np.diag([0.0, PRIOR_TILT_RAD**-2, PRIOR_TILT_RAD**-2, 0.0, PRIOR_SCALE**-2, PRIOR_TOP**-2])


def score(fim: dict, grid: PatchGrid) -> list[dict[str, Any]]:
    rows = []
    for pid, ent in fim.items():
        lat, lon = grid.centre(pid)
        full = _sigmas(ent["all"])
        known = _sigmas(ent["all"][:4, :4])
        vert = _sigmas(ent["vert"][:4, :4]) if ent["n_vert"] else np.array([np.inf] * 3)
        prior = _sigmas(ent["all"] + PRIOR_INFO)
        prior_v = _sigmas(ent["vert"] + PRIOR_INFO) if ent["n_vert"] else np.array([np.inf] * 3)
        rec = {"patch": pid, "lat": float(lat), "lon": float(lon), "area_km2": float(grid.area_km2(pid)),
               "n_oblique": ent["n_obl"], "n_vertical": ent["n_vert"]}
        for tag, s in (("", full), ("_known", known), ("_vert", vert)):
            rec[f"sigma_w{tag}"] = float(s[0])
            rec[f"sigma_ge_deg{tag}"] = float(s[1])
            rec[f"sigma_gn_deg{tag}"] = float(s[2])
            tilt_ok = max(s[1], s[2]) <= GOAL_TILT_DEG
            rec[f"goal{tag}"] = bool(s[0] <= GOAL_DRIFT and tilt_ok)
            rec[f"thresh{tag}"] = bool(s[0] <= THRESH_DRIFT and tilt_ok)
            rec[f"measured{tag}"] = bool(np.isfinite(s[0]))
        rec["sigma_w_prior"] = float(prior[0])
        rec["drift_goal_prior"] = bool(prior[0] <= GOAL_DRIFT)
        rec["sigma_w_prior_vert"] = float(prior_v[0])
        rec["drift_goal_prior_vert"] = bool(prior_v[0] <= GOAL_DRIFT)
        rows.append(rec)
    return rows


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

EARTH_AREA_KM2 = 4.0 * np.pi * (R_EARTH_M / 1e3) ** 2


def run_scenario(sc: Scenario, tables: dict, grid: PatchGrid) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if sc.profile == "diurnal":
        table, night = tables[("day", sc.alt_m)], tables[("night", sc.alt_m)]
    else:
        table, night = tables[(sc.profile, sc.alt_m)], None
    snap_rows, patch_rows = [], []
    union_goal, union_drift, union_vert = set(), set(), set()
    for s in range(sc.snapshots):
        t = s * SNAPSHOT_STEP_S
        pos, vel = sc.states(t)
        out = snapshot(table, pos, vel, grid, sc.noise_factor, night_table=night, t_s=t)
        sc_rows = score(out["fim"], grid)
        st = out["stats"]
        rec = {"scenario": sc.name, "snapshot": s, "t_s": t, "n_sats": len(pos),
               "pairs_in_range": st["pairs_in_range"], "links_oblique": st["links_oblique"],
               "links_vertical": st["links_vertical"], "patches_hit": len(sc_rows)}
        for tag in ("", "_known", "_vert"):
            for key in ("measured", "goal", "thresh"):
                area = sum(r["area_km2"] for r in sc_rows if r[f"{key}{tag}"])
                rec[f"{key}{tag}_km2"] = area
                rec[f"{key}{tag}_n"] = sum(1 for r in sc_rows if r[f"{key}{tag}"])
        for key in ("drift_goal_prior", "drift_goal_prior_vert"):
            rec[f"{key}_km2"] = sum(r["area_km2"] for r in sc_rows if r[key])
        union_goal |= {r["patch"] for r in sc_rows if r["goal"]}
        union_drift |= {r["patch"] for r in sc_rows if r["drift_goal_prior"]}
        union_vert |= {r["patch"] for r in sc_rows if r["drift_goal_prior_vert"]}
        rec["drift_prior_vert_union_km2"] = (float(np.sum(grid.area_km2(np.array(sorted(union_vert), int))))
                                             if union_vert else 0.0)
        rec["goal_union_km2"] = float(np.sum(grid.area_km2(np.array(sorted(union_goal), int)))) if union_goal else 0.0
        rec["drift_prior_union_km2"] = (float(np.sum(grid.area_km2(np.array(sorted(union_drift), int))))
                                        if union_drift else 0.0)
        snap_rows.append(rec)
        for r in sc_rows:
            r.update({"scenario": sc.name, "snapshot": s})
        patch_rows.extend(sc_rows)
    return snap_rows, patch_rows


def summarise(sc: Scenario, snaps, patches) -> dict[str, Any]:
    def mean(key):
        return float(np.mean([r[key] for r in snaps]))
    fin = [p for p in patches if np.isfinite(p["sigma_w"])]
    def med(k):
        return float(np.median([p[k] for p in fin])) if fin else float("nan")
    tilt = [max(p["sigma_ge_deg"], p["sigma_gn_deg"]) for p in fin]
    return {
        "scenario": sc.name, "n_sats": snaps[0]["n_sats"], "altitude_km": sc.alt_m / 1e3,
        "profile": sc.profile, "noise_factor": sc.noise_factor, "snapshots": len(snaps),
        "pairs_in_range": mean("pairs_in_range"), "links_oblique": mean("links_oblique"),
        "links_vertical": mean("links_vertical"), "patches_hit": mean("patches_hit"),
        "measured_km2": mean("measured_km2"), "goal_km2": mean("goal_km2"), "thresh_km2": mean("thresh_km2"),
        "goal_known_km2": mean("goal_known_km2"), "measured_vert_km2": mean("measured_vert_km2"),
        "goal_vert_km2": mean("goal_vert_km2"),
        "drift_prior_km2": mean("drift_goal_prior_km2"), "drift_prior_vert_km2": mean("drift_goal_prior_vert_km2"),
        "goal_earth_pct": 100.0 * mean("goal_km2") / EARTH_AREA_KM2,
        "goal_union_earth_pct": 100.0 * snaps[-1]["goal_union_km2"] / EARTH_AREA_KM2,
        "drift_prior_union_earth_pct": 100.0 * snaps[-1]["drift_prior_union_km2"] / EARTH_AREA_KM2,
        "drift_prior_vert_union_earth_pct": 100.0 * snaps[-1]["drift_prior_vert_union_km2"] / EARTH_AREA_KM2,
        "median_sigma_w": med("sigma_w"),
        "median_sigma_tilt_deg": float(np.median(tilt)) if tilt else float("nan"),
        "median_sigma_w_known": float(np.median([p["sigma_w_known"] for p in fin])) if fin else float("nan"),
    }


def build_tables(scenarios, freqs=FREQS_HZ):
    from tss.l1_propagation.link_table import LinkTable
    prof = profiles()
    out = {}
    for sc in scenarios:
        for p in (("day", "night") if sc.profile == "diurnal" else (sc.profile,)):
            key = (p, sc.alt_m)
            if key not in out:
                out[key] = LinkTable.build(prof[p], R_EARTH_M + sc.alt_m, freqs)
    return out


def _write_csv(path: Path, rows: list[dict[str, Any]], header: str = "") -> Path:
    fields = list(rows[0].keys())
    with open(path, "w", encoding="utf-8") as fh:
        if header:
            fh.write(f"# {header}\n")
        fh.write(",".join(fields) + "\n")
        for r in rows:
            vals = []
            for k in fields:
                v = r.get(k, "")
                if isinstance(v, (bool, np.bool_)):
                    vals.append(str(bool(v)))
                elif isinstance(v, str):
                    vals.append(f'"{v}"')
                elif isinstance(v, (float, np.floating)):
                    vals.append(f"{float(v):.6g}")
                else:
                    vals.append(str(v))
            fh.write(",".join(vals) + "\n")
    return path


def _table_rows(tables) -> list[dict[str, Any]]:
    rows = []
    for (prof, alt), t in tables.items():
        for f, ft in sorted(t.tables.items()):
            rows.append({"profile": prof, "altitude_km": alt / 1e3, "f_mhz": f / 1e6,
                         "theta_lo_km": ft.theta_lo * R_EARTH_M / 1e3, "theta_hi_km": ft.theta_hi * R_EARTH_M / 1e3,
                         "drift_coeff_vertical": float(-ft.d(0.0)) if ft.theta_lo == 0.0 else float("nan")})
    return rows


def _figure(summary, patches, scenarios, out_dir: Path) -> list[Path]:
    try:
        import matplotlib.patches as mpatches
    except ImportError:  # pragma: no cover
        return []
    from tss.figures import style

    style.use()
    fig, axes = style.figure(2, 2, figsize=(10.0, 7.6))
    col_goal, col_drift, col_none = style.SERIES[0], style.SERIES[3], "#d9d8d4"

    def ylat(lat):
        return np.sin(np.radians(lat))

    def patch_map(ax, name, title, extent=None, sats=None):
        rows = [p for p in patches if p["scenario"] == name and p["snapshot"] == 0]
        for p in rows:
            g = PATCH_DEG
            w = 360.0 / max(1, round(360.0 * np.cos(np.radians(p["lat"])) / g))
            col = col_goal if p["goal"] else (col_drift if p["drift_goal_prior"] else col_none)
            y0, y1 = ylat(p["lat"] - g / 2), ylat(p["lat"] + g / 2)
            ax.add_patch(mpatches.Rectangle((p["lon"] - w / 2, y0), w, y1 - y0, color=col, lw=0))
        if sats is not None:
            lat, lon = _llh(sats)
            ax.plot(lon, ylat(lat), "k^", ms=3)
        lo, hi, la, lb = extent if extent else (-180, 180, -90, 90)
        ax.set_xlim(lo, hi)
        ax.set_ylim(ylat(la), ylat(lb))
        ticks = [t for t in (-90, -60, -30, 0, 10, 20, 30, 40, 60, 90) if la <= t <= lb]
        if not extent:
            ticks = [-90, -60, -30, 0, 30, 60, 90]
        ax.set_yticks([ylat(t) for t in ticks])
        ax.set_yticklabels([f"{t}" for t in ticks])
        ax.set_xlabel("longitude (deg)")
        ax.set_ylabel("latitude (deg, equal-area scale)")
        ax.set_title(title, fontsize=9)

    sc = {s.name: s for s in scenarios}
    pos, _ = sc["cluster9_mixed_polar"].states(0.0)
    patch_map(axes[0, 0], "cluster9_mixed_polar", "(a) Nine satellites: Chartier grid, mixed headings",
              (-160, -120, 0, 40), pos)
    pos, _ = sc["star128"].states(0.0)
    patch_map(axes[0, 1], "star128", "(b) 128 satellites, polar star, one instant", None, pos)
    handles = [mpatches.Patch(color=col_goal, label="drift <= 25 m/s and both tilts <= 0.29 deg"),
               mpatches.Patch(color=col_drift, label="drift <= 25 m/s given a 1 deg tilt prior"),
               mpatches.Patch(color=col_none, label="links reflect here, neither")]
    axes[0, 0].legend(handles=handles, fontsize=6, frameon=True, loc="lower left")

    ax = axes[1, 0]
    order = [s for s in summary if "noise2" not in s["scenario"]]
    names = [s["scenario"] for s in order]
    y = np.arange(len(names))
    ax.barh(y + 0.27, [s["goal_union_earth_pct"] for s in order], 0.27, color=col_goal,
            label="drift and tilt, all links")
    ax.barh(y, [s["drift_prior_union_earth_pct"] for s in order], 0.27, color=col_drift,
            label="drift with tilt prior, all links")
    ax.barh(y - 0.27, [s["drift_prior_vert_union_earth_pct"]
                       for s in order], 0.27, color=style.SERIES[2], label="drift with tilt prior, vertical only")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=6)
    ax.set_xscale("symlog", linthresh=0.1)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Earth measured to the goal at least once in 100 min (%)")
    ax.set_title("(c) Coverage in one orbit", fontsize=9)
    ax.legend(fontsize=6, frameon=False, loc="lower right")

    ax = axes[1, 1]
    for k, name in enumerate(("cluster9_mixed_polar", "cluster9_comoving_polar", "star128", "delta128_53",
                              "star256")):
        vals = np.array([p["sigma_w_prior"] for p in patches if p["scenario"] == name])
        if vals.size:
            fin = np.sort(vals[np.isfinite(vals)])
            ax.semilogx(fin, np.arange(1, fin.size + 1) / vals.size, "-", color=style.SERIES[k], label=name)
    ax.axvline(GOAL_DRIFT, color=style.INK_MUTED, ls="--", lw=0.8)
    ax.axvline(THRESH_DRIFT, color=style.INK_MUTED, ls=":", lw=0.8)
    ax.set_xlabel("one-sigma drift with a 1 deg tilt prior (m/s)")
    ax.set_ylabel("fraction of patches with a link")
    ax.set_title("(d) Drift precision where links reflect", fontsize=9)
    ax.legend(fontsize=6, frameon=False, loc="upper left")
    fig.tight_layout()
    return style.savefig(fig, out_dir / "route_two_network", formats=("png", "pdf"))


@register_study(
    "route_two_network",
    description="Route two, phase 2: drift and tilt over an area from a constellation's oblique links, 9 vs 128",
    work_packages=["WP5", "WP7"],
    claims=["C3"],
    status="implemented",
    requirements=["SCI-15", "SCI-23", "SCI-27"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    opts = dict(config.get("route_two_network", {}))
    only = opts.get("scenarios")
    scenarios = [s for s in SCENARIOS if not only or s.name in only]
    if "snapshots" in opts:
        scenarios = [Scenario(**{**s.__dict__, "snapshots": min(s.snapshots, int(opts["snapshots"]))})
                     for s in scenarios]
    grid = PatchGrid(PATCH_DEG)
    tables = build_tables(scenarios)
    summary, all_snaps, all_patches = [], [], []
    for sc in scenarios:
        snaps, patches = run_scenario(sc, tables, grid)
        all_snaps.extend(snaps)
        all_patches.extend(patches)
        summary.append(summarise(sc, snaps, patches))
    hdr = (f"links 2-20 MHz; SNR {SNR_REF_DB} dB at {P_REF_M/1e3:.0f} km, >= {SNR_MIN_DB} dB; "
           f"sigma {SIGMA_RATE_REF} m/s, {SIGMA_GROUP_REF} m at 1000 km; patches {PATCH_DEG} deg equal area; "
           f"goal drift {GOAL_DRIFT} m/s and tilt {GOAL_TILT_DEG} deg per component")
    outputs = [
        _write_csv(out_dir / "route_two_network_summary.csv", summary, hdr),
        _write_csv(out_dir / "route_two_network_snapshots.csv", all_snaps, hdr),
        _write_csv(out_dir / "route_two_network_patches.csv", [p for p in all_patches if p["snapshot"] == 0],
                   hdr + "; first snapshot only"),
        _write_csv(out_dir / "route_two_network_branches.csv", _table_rows(tables),
                   "usable central-angle range of the reflected branch, after the 1% margins"),
    ]
    if any(s.name == "cluster9_mixed_polar" for s in scenarios) and any(s.name == "star128" for s in scenarios):
        outputs += _figure(summary, all_patches, scenarios, out_dir)
    by = {s["scenario"]: s for s in summary}
    metrics = {s["scenario"]: {k: s[k] for k in ("goal_km2", "thresh_km2", "goal_vert_km2", "median_sigma_w",
                                                  "median_sigma_tilt_deg", "links_oblique")} for s in summary}
    lines = ["| scenario | links (obl/vert) | drift+tilt, instant (10^6 km^2) | drift+tilt, 100 min (% Earth) | "
             "drift w/ 1 deg prior, instant (10^6 km^2) | same, 100 min (% Earth) | vertical only, 100 min (%) |",
             "|---|---|---|---|---|---|---|"]
    for s in summary:
        lines.append(f"| {s['scenario']} | {s['links_oblique']:.0f}/{s['links_vertical']:.0f} | "
                     f"{s['goal_km2']/1e6:.3g} | {s['goal_union_earth_pct']:.3g} | {s['drift_prior_km2']/1e6:.3g} | "
                     f"{s['drift_prior_union_earth_pct']:.3g} | {s['drift_prior_vert_union_earth_pct']:.3g} |")
    text = ("Route two, phase 2: drift and tilt over an area from oblique links at one instant.\n\n"
            + "\n".join(lines) + "\n\n" + hdr + ".\n")
    if "star128" in by and "cluster9_mixed_polar" in by:
        text += (f"\n128 polar star: {by['star128']['goal_earth_pct']:.2f}% of the Earth meets the goal at any "
                 f"instant; the nine-satellite cluster {by['cluster9_mixed_polar']['goal_km2']/1e6:.3g} million km^2.\n")
    return StudyResult(outputs=outputs, summary=text, metrics=metrics)
