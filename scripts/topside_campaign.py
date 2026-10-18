#!/usr/bin/env python3
"""Map and trace the topside difference between PyIRI and Fortran IRI-2020.

Reproduces ``docs/topside-model-campaign.md`` and the Technical Volume section
on the two implementations of IRI-2020's NeQuick topside.

    python scripts/topside_campaign.py grid    OUT_DIR
    python scripts/topside_campaign.py isis    OUT_DIR --crc CRC.asc --t1 T1.csv
    python scripts/topside_campaign.py analyse OUT_DIR

``grid`` runs both models over a global grid: latitude -80 to 80 deg in 10 deg
steps, longitude -180 to 150 deg in 30 deg steps, every 3 h of universal time,
on the 15th of March, June, September and December, at F10.7 = 70, 120 and
180 sfu.  That is 19,584 profiles per model, 100 to 2000 km in 10 km steps.
PyIRI runs here; Fortran IRI-2020 runs through PHaRLAP's iri2020 mex in a
MATLAB batch session (``scripts/topside_campaign_iri.m``), with the full
iono_extra vector kept so that the topside can be rebuilt from the model's own
parameters.

``isis`` evaluates PyIRI at the 5,118 CRC hand-scaled ISIS-1 profiles that
``Dissertation/ISISValidation/t1_compare.py`` paired with IRI-2020, at the same
time, place, altitude and solar level, and fits each measured profile with the
NeQuick topside shape.

``analyse`` writes ``summary.json`` and the figures.

Both models are run with IRI-2020's default options: URSI foF2, SHU-2015 hmF2,
ABT-2009 bottomside, NeQuick topside.  PHaRLAP takes a user R12, set here from
the F10.7 by inverting the same relation the mex uses to go the other way,
F10.7 = 63.75 + 0.728 R12 + 0.00089 R12^2, which is also the relation PyIRI
uses internally (``PyIRI.main_library.F107_2_R12``).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

YEAR, DAY = 2024, 15
MONTHS = (3, 6, 9, 12)
F107S = (70.0, 120.0, 180.0)
UTS = np.arange(0.0, 24.0, 3.0)
LATS = np.arange(-80.0, 81.0, 10.0)
LONS = np.arange(-180.0, 180.0, 30.0)
ALT_KM = np.arange(100.0, 2000.0 + 0.1, 10.0)
REPORT_KM = (500.0, 600.0, 700.0, 800.0, 900.0)
PHARLAP = r"C:\Users\Andrew\PhD\Chartier Verification\pharlap_4.7.4"
PYIRI_OPTS = dict(foF2_coeff="URSI", hmF2_model="SHU2015", coord="GEO", old_output=False)

# iono_extra columns (1-based in PHaRLAP's help, 0-based here)
EX = dict(nmf2=0, hmf2=1, b0=9, sza=22, dip=24, diplat=25, modip=26, r12=32,
          b1=34, m3000=35, ig12=38, f107=40)


# ---------------------------------------------------------------- solar index
def f107_from_r12(r12):
    """F10.7 from R12 as the PHaRLAP mex and PyIRI both compute it."""
    r12 = np.asarray(r12, dtype=float)
    return 63.75 + r12 * (0.728 + 0.00089 * r12)


def r12_from_f107(f107):
    """Inverse of :func:`f107_from_r12`."""
    a, b, c = 0.00089, 0.728, 63.75 - np.asarray(f107, dtype=float)
    return (-b + np.sqrt(b * b - 4.0 * a * c)) / (2.0 * a)


# ------------------------------------------------------- the NeQuick topside
G, RFAC = 0.125, 100.0


def b2bot(fof2_mhz, m3000):
    """NeQuick bottomside thickness, km (irisub.for itopn=2; PyIRI thickness_F2)."""
    dndhmx = 0.01 * np.exp(-3.467 + 1.714 * np.log(fof2_mhz) + 2.02 * np.log(m3000))
    return 0.04774 * fof2_mhz ** 2 / dndhmx


def shape_factor(fof2_mhz, hmf2_km, b2bot_km, r12):
    """Coisson et al. (2006) shape factor k, as both codes write it."""
    return (3.22 - 0.0538 * fof2_mhz - 0.00664 * hmf2_km
            + 0.113 * hmf2_km / b2bot_km + 0.00257 * r12)


def h0_iri(fof2_mhz, m3000, hmf2_km, r12):
    """IRI-2020 (irisub.for, itopn=2): H0 = k' B2bot, k smoothly floored at 1."""
    bb = b2bot(fof2_mhz, m3000)
    k = shape_factor(fof2_mhz, hmf2_km, bb, r12)
    ee = np.exp(2.0 * (k - 1.0))
    return (k * ee + 1.0) / (ee + 1.0) * bb


def h0_pyiri(fof2_mhz, m3000, hmf2_km, r12):
    """PyIRI 0.1.7 (sh_library.thickness_F2): H0 = k B2bot / v(x).

    x = (k B2bot - 150)/100 and v = 0.041163 x^2 - 0.183981 x + 1.424472 is the
    auxiliary step of the NeQuick topside thickness as NeQuick G specifies it
    (European Commission 2016; JRC reference code,
    NeQuickG_JRC_iono_F2_layer.c, F2_layer_exosphere_adjust).
    """
    bb = b2bot(fof2_mhz, m3000)
    k = shape_factor(fof2_mhz, hmf2_km, bb, r12)
    x = (k * bb - 150.0) / 100.0
    return (100.0 * x + 150.0) / (0.041163 * x * x - 0.183981 * x + 1.424472)


def v_factor(kb2bot_km):
    """The divisor v of the NeQuick G step, as a function of k B2bot."""
    x = (np.asarray(kb2bot_km, dtype=float) - 150.0) / 100.0
    return 0.041163 * x * x - 0.183981 * x + 1.424472


def topside_shape(h_km, hmf2_km, h0_km):
    """Ne/NmF2 of the NeQuick semi-Epstein topside (TOPQ in irifun.for)."""
    dh = np.asarray(h_km, dtype=float) - hmf2_km
    g1 = G * dh
    z = dh / (h0_km * (1.0 + RFAC * g1 / (RFAC * h0_km + g1)))
    e = np.exp(-np.abs(z))
    return 4.0 * e / (1.0 + e) ** 2


def wpercentile(x, w, q):
    """Weighted percentile (q in percent) of the finite values of x."""
    x = np.asarray(x, dtype=float)
    w = np.broadcast_to(np.asarray(w, dtype=float), x.shape).ravel()
    x = x.ravel()
    m = np.isfinite(x) & np.isfinite(w)
    x, w = x[m], w[m]
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return float(np.interp(q / 100.0, c, x))


def _pyiri():
    import PyIRI
    import PyIRI.sh_library as sh
    return PyIRI, sh


# ---------------------------------------------------------------------- grid
def grid_sites():
    lat, lon = np.meshgrid(LATS, LONS, indexing="ij")
    return lat.ravel(), lon.ravel()


def cmd_grid(out: Path, matlab: str = "matlab", skip_iri: bool = False) -> None:
    out.mkdir(parents=True, exist_ok=True)
    PyIRI, sh = _pyiri()
    lat, lon = grid_sites()
    nM, nS, nT, nG, nV = len(MONTHS), len(F107S), UTS.size, lat.size, ALT_KM.size
    ne = np.empty((nM, nS, nT, nG, nV), dtype=np.float32)
    keys = ("Nm", "fo", "M3000", "hm", "B_top", "B_bot", "B0", "B1")
    par = {k: np.empty((nM, nS, nT, nG)) for k in keys}
    modip = None
    t0 = time.time()
    for i, month in enumerate(MONTHS):
        for j, f107 in enumerate(F107S):
            F2, _F1, _E, _Es, _sun, mag, edp = sh.IRI_density_1day(
                YEAR, month, DAY, UTS, lon, lat, ALT_KM, f107, **PYIRI_OPTS)
            ne[i, j] = np.transpose(np.asarray(edp), (0, 2, 1))
            for k in keys:
                par[k][i, j] = np.asarray(F2[k])
            modip = np.asarray(mag["modip"], dtype=float).ravel()
    print(f"PyIRI grid: {nM * nS * nT * nG} profiles in {time.time() - t0:.1f} s")
    np.savez_compressed(out / "pyiri_grid.npz", ne=ne, modip=modip, lat=lat, lon=lon,
                        months=np.array(MONTHS), f107=np.array(F107S), uts=UTS,
                        alt_km=ALT_KM, pyiri_version=PyIRI.__version__, **par)

    rows = []
    for month in MONTHS:
        for f107 in F107S:
            r12 = float(r12_from_f107(f107))
            for ut in UTS:
                for la, lo in zip(lat, lon, strict=True):
                    rows.append((YEAR, month, DAY, ut, 0.0, la, lo, r12))
    np.savetxt(out / "iri_points.csv", np.array(rows), delimiter=",", fmt="%.6f")
    np.savetxt(out / "alt_km.csv", ALT_KM, fmt="%.1f")
    if not skip_iri:
        run_matlab(out / "iri_points.csv", out / "alt_km.csv", out / "iri_grid.mat", matlab)


def run_matlab(points: Path, alt: Path, out_mat: Path, matlab: str = "matlab") -> None:
    sd = str(ROOT / "scripts").replace("\\", "/")
    p, a, o = (str(x).replace("\\", "/") for x in (points, alt, out_mat))
    root = PHARLAP.replace("\\", "/")
    call = f"addpath('{sd}'); topside_campaign_iri('{p}','{a}','{o}','{root}')"
    t0 = time.time()
    subprocess.run([matlab, "-batch", call], check=True)
    print(f"IRI-2020 (PHaRLAP) grid written in {time.time() - t0:.1f} s")


# ---------------------------------------------------------------------- isis
def load_isis_selection(crc: Path, t1: Path):
    """The T1 selection and order, with profiles, checked against t1_results."""
    sys.path.insert(0, str(crc.parent))
    import crc_read  # Dissertation/ISISValidation/crc_read.py

    _, recs = crc_read.read_crc(str(crc), want_profiles=True)
    sel = [x for x in recs if 500 <= x["h_sat_km"] <= 900 and 0 <= x["iquality"] <= 7
           and np.isfinite(x["fn_sat_mhz"]) and x["rz12"] is not None]
    for x in sel:
        x["rz_bucket"] = int(round(x["rz12"] / 5.0) * 5)
    order = [x for b in sorted({y["rz_bucket"] for y in sel}) for x in sel if x["rz_bucket"] == b]
    T = np.genfromtxt(t1, delimiter=",", names=True)
    if len(T) != len(order):
        raise RuntimeError(f"t1_results has {len(T)} rows, selection has {len(order)}")
    if not np.allclose(T["h_sat_km"], [x["h_sat_km"] for x in order], atol=0.06):
        raise RuntimeError("t1_results is not in the order the selection gives")
    return order, T


def fit_h0(h_km, ne, hmf2_km, nmf2):
    """Least-squares H0 of the NeQuick topside through one measured profile."""
    from scipy.optimize import minimize_scalar

    def cost(h0):
        r = np.log(ne) - np.log(nmf2 * topside_shape(h_km, hmf2_km, h0))
        return float(np.sum(r * r))

    res = minimize_scalar(cost, bounds=(5.0, 500.0), method="bounded",
                          options={"xatol": 0.05})
    rms = np.sqrt(res.fun / max(1, h_km.size))
    return float(res.x), float(rms)


def cmd_isis(out: Path, crc: Path, t1: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    PyIRI, sh = _pyiri()
    order, T = load_isis_selection(crc, t1)
    n = len(order)
    res = {k: np.full(n, np.nan) for k in
           ("ne_sat_pyiri", "Nm", "fo", "M3000", "hm", "B_top", "modip", "f107")}
    groups: dict = {}
    for i, x in enumerate(order):
        d = dt.date(x["year"], 1, 1) + dt.timedelta(days=int(x["doy"]) - 1)
        groups.setdefault((d, x["rz_bucket"]), []).append(i)
    t0 = time.time()
    for gi, ((d, rzb), idx) in enumerate(groups.items()):
        f107 = float(f107_from_r12(rzb))
        ut = np.array([order[i]["ut_hour"] for i in idx])
        la = np.array([order[i]["lat"] for i in idx], dtype=float)
        lo = (np.array([order[i]["lon"] for i in idx], dtype=float) + 180.0) % 360.0 - 180.0
        hs = np.array([order[i]["h_sat_km"] for i in idx], dtype=float)
        ut_u, ti = np.unique(ut, return_inverse=True)
        h_u, hi = np.unique(hs, return_inverse=True)
        F2, _F1, _E, _Es, _sun, mag, edp = sh.IRI_density_1day(
            d.year, d.month, d.day, ut_u, lo, la, h_u, f107, **PYIRI_OPTS)
        g = np.arange(len(idx))
        ix = np.asarray(idx)
        res["ne_sat_pyiri"][ix] = np.asarray(edp)[ti, hi, g]
        for k in ("Nm", "fo", "M3000", "hm", "B_top"):
            res[k][ix] = np.asarray(F2[k])[ti, g]
        res["modip"][ix] = np.asarray(mag["modip"], dtype=float).ravel()[g]
        res["f107"][ix] = f107
        if gi % 100 == 0:
            print(f"  isis group {gi + 1}/{len(groups)}  {time.time() - t0:.0f} s", flush=True)
    print(f"PyIRI at ISIS-1 points: {n} in {len(groups)} calls, {time.time() - t0:.0f} s")

    # the measured profiles: observed peak, and the NeQuick H0 that fits each
    obs = {k: np.full(n, np.nan) for k in ("nmf2_obs", "hmf2_obs", "h0_fit", "fit_rms", "n_pts", "depth_km")}
    for i, x in enumerate(order):
        h, ne = np.asarray(x["h_km"], float), np.asarray(x["ne_cm3"], float) * 1.0e6
        j = int(np.argmax(ne))
        obs["nmf2_obs"][i], obs["hmf2_obs"][i] = ne[j], h[j]
        obs["depth_km"][i] = h[0] - h[j]
        use = (h >= h[j] + 30.0) & np.isfinite(ne) & (ne > 0)
        obs["n_pts"][i] = use.sum()
        if use.sum() >= 4 and obs["depth_km"][i] >= 150.0:
            obs["h0_fit"][i], obs["fit_rms"][i] = fit_h0(h[use], ne[use], h[j], ne[j])

    cols = dict(
        year=np.array([x["year"] for x in order]), doy=np.array([x["doy"] for x in order]),
        ut=np.array([x["ut_hour"] for x in order]), lt=np.array([x["lt_hour"] for x in order]),
        lat=np.array([x["lat"] for x in order], float), lon=np.array([x["lon"] for x in order], float),
        h_sat=np.array([x["h_sat_km"] for x in order], float),
        rz12=np.array([x["rz12"] for x in order], float),
        rz_bucket=np.array([x["rz_bucket"] for x in order], float),
        quality=np.array([x["iquality"] for x in order], float),
        ne_sat_obs=np.array([x["ne_sat_cm3"] for x in order], float) * 1.0e6,
        fn_obs=T["fn_obs_MHz"], fn_iri2020=T["fn_NeQuick_MHz"],
        hmf2_hdr=np.array([x["hmf2_iri"] if x["hmf2_iri"] else np.nan for x in order], float),
    )
    np.savez_compressed(out / "isis_pyiri.npz", pyiri_version=PyIRI.__version__,
                        **cols, **res, **obs)


# ------------------------------------------------------------------ analyse
def _wstats(x, w):
    x = np.asarray(x, dtype=float)
    w = np.broadcast_to(np.asarray(w, dtype=float), x.shape)
    return {"median": wpercentile(x, w, 50), "p10": wpercentile(x, w, 10),
            "p90": wpercentile(x, w, 90), "n": int(np.isfinite(x).sum())}


def _wmean(x, w):
    x = np.asarray(x, dtype=float)
    w = np.broadcast_to(np.asarray(w, dtype=float), x.shape)
    m = np.isfinite(x)
    return float(np.sum(x[m] * w[m]) / np.sum(w[m]))


def fp_mhz(ne_m3):
    return 8.98e-6 * np.sqrt(np.clip(np.asarray(ne_m3, dtype=float), 0.0, None))


def load_grid(out: Path):
    import scipy.io as sio

    P = dict(np.load(out / "pyiri_grid.npz", allow_pickle=False))
    S = sio.loadmat(out / "iri_grid.mat")
    shp = P["Nm"].shape
    Q = {"ne": S["ne"].astype(float).reshape(shp + (-1,)),
         "n_failed": int(np.asarray(S["n_failed"]).ravel()[0]),
         "elapsed_s": float(np.asarray(S["elapsed_s"]).ravel()[0])}
    ex = S["extra"].reshape(shp + (-1,))
    for k, c in EX.items():
        Q[k] = ex[..., c]
    Q["fo"] = np.sqrt(Q["nmf2"] / 1.24e10)
    return P, Q


def grid_analysis(out: Path) -> dict:
    P, Q = load_grid(out)
    alt = P["alt_km"]
    lat, lon = P["lat"], P["lon"]
    nM, nS, nT, nG = P["Nm"].shape
    w = np.broadcast_to(np.cos(np.deg2rad(lat)), (nM, nS, nT, nG))
    r12 = r12_from_f107(P["f107"])[None, :, None, None] * np.ones((nM, nS, nT, nG))
    lt = (P["uts"][:, None] + lon[None, :] / 15.0) % 24.0
    lt = np.broadcast_to(lt, (nM, nS, nT, nG))
    amod = np.broadcast_to(np.abs(P["modip"]), (nM, nS, nT, nG))
    s = {"grid": {"profiles_per_model": int(nM * nS * nT * nG), "months": P["months"].tolist(),
                  "f107": P["f107"].tolist(), "r12": r12_from_f107(P["f107"]).round(2).tolist(),
                  "iri_failed": Q["n_failed"], "iri_elapsed_s": Q["elapsed_s"],
                  "pyiri_version": str(P["pyiri_version"])}}

    # 1. each code's topside rebuilt from its own parameters
    h0i = h0_iri(Q["fo"], Q["m3000"], Q["hmf2"], Q["r12"])
    h0p_rebuilt = h0_pyiri(P["fo"], P["M3000"], P["hm"], r12)
    above_i = alt[None, None, None, None, :] > Q["hmf2"][..., None] + 50.0
    above_p = alt[None, None, None, None, :] > P["hm"][..., None] + 50.0
    rb_i = Q["nmf2"][..., None] * topside_shape(alt, Q["hmf2"][..., None], h0i[..., None])
    rb_p = P["Nm"][..., None] * topside_shape(alt, P["hm"][..., None], P["B_top"][..., None])
    qi = (rb_i / Q["ne"])[above_i]
    qp = (rb_p / P["ne"])[above_p]
    s["rebuild"] = {
        "iri_rebuilt_over_actual": [float(np.nanpercentile(qi, q)) for q in (0.5, 50, 99.5)],
        "pyiri_rebuilt_over_actual": [float(np.nanpercentile(qp, q)) for q in (0.5, 50, 99.5)],
        "pyiri_Btop_over_formula": [float(np.nanpercentile(P["B_top"] / h0p_rebuilt, q)) for q in (0.5, 50, 99.5)],
        "iri_r12_used": float(np.nanmedian(Q["r12"])),
    }

    # 2. the parameters behind the topside
    bb_i, bb_p = b2bot(Q["fo"], Q["m3000"]), b2bot(P["fo"], P["M3000"])
    k_i = shape_factor(Q["fo"], Q["hmf2"], bb_i, Q["r12"])
    k_p = shape_factor(P["fo"], P["hm"], bb_p, r12)
    h0p_iriform = h0_iri(P["fo"], P["M3000"], P["hm"], r12)
    par = {
        "nmf2_ratio": P["Nm"] / Q["nmf2"], "hmf2_diff_km": P["hm"] - Q["hmf2"],
        "m3000_ratio": P["M3000"] / Q["m3000"], "b2bot_ratio": bb_p / bb_i,
        "k_ratio": k_p / k_i, "h0_ratio": P["B_top"] / h0i,
        "h0_formula_ratio": P["B_top"] / h0p_iriform, "h0_inputs_ratio": h0p_iriform / h0i,
        "v_factor": v_factor(k_p * bb_p), "kb2bot_pyiri_km": k_p * bb_p,
        "h0_iri_km": h0i, "h0_pyiri_km": P["B_top"], "b2bot_iri_km": bb_i, "k_iri": k_i,
    }
    s["parameters"] = {name: {str(int(f)): _wstats(v[:, j], w[:, j]) for j, f in enumerate(P["f107"])}
                       for name, v in par.items()}
    s["parameters_all"] = {name: _wstats(v, w) for name, v in par.items()}

    # 3. the density at the orbit altitudes, and 4. its decomposition
    s["altitudes"] = {}
    ratio_prof = {}
    for hk in REPORT_KM + tuple(np.arange(400.0, 2000.1, 100.0)):
        iz = int(np.argmin(np.abs(alt - hk)))
        nei, nep = Q["ne"][..., iz], P["ne"][..., iz]
        ratio = nep / nei
        S_p = topside_shape(hk, P["hm"], P["B_top"])
        S_pi = topside_shape(hk, P["hm"], h0p_iriform)
        S_i = topside_shape(hk, Q["hmf2"], h0i)
        S_mix = topside_shape(hk, P["hm"], h0i)
        terms = {
            "peak_density": np.log(P["Nm"] / Q["nmf2"]),
            "thickness_formula": np.log(S_p / S_pi),
            "hmf2_shift": np.log(S_mix / S_i),
            "thickness_inputs": np.log(S_pi / S_mix),
        }
        corrected = P["Nm"] * S_pi / nei
        ratio_prof[hk] = {str(int(f)): _wstats(ratio[:, j], w[:, j]) for j, f in enumerate(P["f107"])}
        ratio_prof[hk]["corrected"] = {str(int(f)): _wstats(corrected[:, j], w[:, j]) for j, f in enumerate(P["f107"])}
        if hk not in REPORT_KM:
            continue
        d = {"ratio": {str(int(f)): _wstats(ratio[:, j], w[:, j]) for j, f in enumerate(P["f107"])},
             "ratio_all": _wstats(ratio, w),
             "corrected_ratio_all": _wstats(corrected, w),
             "corrected_ratio": {str(int(f)): _wstats(corrected[:, j], w[:, j]) for j, f in enumerate(P["f107"])},
             "log_ratio_gm": float(np.exp(_wmean(np.log(ratio), w))),
             "terms_gm": {k: float(np.exp(_wmean(v, w))) for k, v in terms.items()},
             "terms_gm_by_f107": {str(int(f)): {k: float(np.exp(_wmean(v[:, j], w[:, j]))) for k, v in terms.items()}
                                  for j, f in enumerate(P["f107"])},
             "fp_p05_mhz": {}, "fp_p50_mhz": {}}
        fpi, fpp, fpc = fp_mhz(nei), fp_mhz(nep), fp_mhz(P["Nm"] * S_pi)
        for j, f in enumerate(P["f107"]):
            key = str(int(f))
            d["fp_p05_mhz"][key] = {"iri2020": wpercentile(fpi[:, j], w[:, j], 5),
                                    "pyiri": wpercentile(fpp[:, j], w[:, j], 5),
                                    "pyiri_iri_h0": wpercentile(fpc[:, j], w[:, j], 5)}
            d["fp_p50_mhz"][key] = {"iri2020": wpercentile(fpi[:, j], w[:, j], 50),
                                    "pyiri": wpercentile(fpp[:, j], w[:, j], 50)}
        d["fp_p05_mhz"]["all"] = {"iri2020": wpercentile(fpi, w, 5), "pyiri": wpercentile(fpp, w, 5),
                                  "pyiri_iri_h0": wpercentile(fpc, w, 5)}
        # where it is largest
        by = {}
        for name, m in {"day_08_16": (lt >= 8) & (lt < 16), "night_20_04": (lt >= 20) | (lt < 4),
                        "modip_0_20": amod < 20, "modip_20_50": (amod >= 20) & (amod < 50),
                        "modip_50_90": amod >= 50}.items():
            by[name] = {str(int(f)): wpercentile(ratio[:, j][m[:, j]], w[:, j][m[:, j]], 50)
                        for j, f in enumerate(P["f107"])}
        for i, mo in enumerate(P["months"]):
            by[f"month_{int(mo)}"] = {str(int(f)): wpercentile(ratio[i, j], w[i, j], 50)
                                      for j, f in enumerate(P["f107"])}
        d["median_by"] = by
        s["altitudes"][str(int(hk))] = d
    s["_ratio_profile"] = {str(int(k)): v for k, v in ratio_prof.items()}
    return s, P, Q, dict(h0i=h0i, h0p_iriform=h0p_iriform, k_p=k_p, bb_p=bb_p, w=w, lt=lt, amod=amod)


BANDS = ((500.0, 600.0), (600.0, 700.0), (700.0, 800.0), (800.0, 900.0))


def isis_analysis(out: Path):
    J = dict(np.load(out / "isis_pyiri.npz", allow_pickle=False))
    hs, rz, rzb = J["h_sat"], J["rz12"], J["rz_bucket"]
    ne_obs = J["ne_sat_obs"]
    fn_obs, fn_iri = J["fn_obs"], J["fn_iri2020"]
    fn_py = fp_mhz(J["ne_sat_pyiri"])
    h0_py_rebuilt = h0_pyiri(J["fo"], J["M3000"], J["hm"], rzb)
    h0_py_iri = h0_iri(J["fo"], J["M3000"], J["hm"], rzb)
    fn_pyc = fp_mhz(J["Nm"] * topside_shape(hs, J["hm"], h0_py_iri))
    fo_o = np.sqrt(J["nmf2_obs"] / 1.24e10)
    bb_o = b2bot(fo_o, J["M3000"])
    k_o = shape_factor(fo_o, J["hmf2_obs"], bb_o, rzb)
    h0o_i = h0_iri(fo_o, J["M3000"], J["hmf2_obs"], rzb)
    h0o_p = h0_pyiri(fo_o, J["M3000"], J["hmf2_obs"], rzb)
    ne_oi = J["nmf2_obs"] * topside_shape(hs, J["hmf2_obs"], h0o_i)
    ne_op = J["nmf2_obs"] * topside_shape(hs, J["hmf2_obs"], h0o_p)
    fit_ok = np.isfinite(J["h0_fit"]) & (J["depth_km"] >= 150.0)
    lt = J["lt"]
    night = (lt >= 20) | (lt < 4)
    day = (lt >= 8) & (lt < 16)
    regimes = {"solar_max": rz >= 100, "solar_min": rz <= 30, "all": np.ones_like(rz, bool)}

    def db(model_fn):
        return 20.0 * np.log10(model_fn / fn_obs)

    s = {"n": int(hs.size), "n_profile_fits": int(fit_ok.sum()),
         "checks": {
             "pyiri_Btop_over_formula_median": float(np.nanmedian(J["B_top"] / h0_py_rebuilt)),
             "pyiri_with_iri_h0_over_iri2020_fn_median": float(np.nanmedian(fn_pyc / fn_iri)),
             "pyiri_over_iri2020_ne_median": float(np.nanmedian((fn_py / fn_iri) ** 2)),
             "hmf2_obs_minus_iri_header_km_median": float(np.nanmedian(J["hmf2_obs"] - J["hmf2_hdr"])),
             "profile_depth_km_median": float(np.nanmedian(J["depth_km"])),
             "fit_rms_ln_median": float(np.nanmedian(J["fit_rms"][fit_ok])),
             "year_counts": {str(int(y)): int((J["year"] == y).sum()) for y in np.unique(J["year"])},
         },
         "bands": {}}
    for name, reg in regimes.items():
        s["bands"][name] = {}
        for lo, hi in BANDS:
            m = reg & (hs >= lo) & ((hs < hi) if hi < 900 else (hs <= hi))
            if m.sum() < 20:
                continue
            mo = m & fit_ok
            s["bands"][name][f"{int(lo)}-{int(hi)}"] = {
                "n": int(m.sum()), "n_fit": int(mo.sum()),
                "median_db": {"iri2020": float(np.median(db(fn_iri)[m])),
                              "pyiri": float(np.median(db(fn_py)[m])),
                              "pyiri_iri_h0": float(np.median(db(fn_pyc)[m])),
                              "obs_peak_iri_h0": float(np.median(10 * np.log10(ne_oi / ne_obs)[mo])),
                              "obs_peak_nequickg_h0": float(np.median(10 * np.log10(ne_op / ne_obs)[mo]))},
                "p05_mhz": {"observed": float(np.percentile(fn_obs[m], 5)),
                            "iri2020": float(np.percentile(fn_iri[m], 5)),
                            "pyiri": float(np.percentile(fn_py[m], 5)),
                            "pyiri_iri_h0": float(np.percentile(fn_pyc[m], 5))},
            }
    r_i, r_p = J["h0_fit"] / h0o_i, J["h0_fit"] / h0o_p
    groups = {"all": fit_ok, "solar_max": fit_ok & (rz >= 100), "solar_min": fit_ok & (rz <= 30),
              "day_08_16": fit_ok & day, "night_20_04": fit_ok & night,
              "abs_lat_lt_30": fit_ok & (np.abs(J["lat"]) < 30),
              "abs_lat_30_60": fit_ok & (np.abs(J["lat"]) >= 30) & (np.abs(J["lat"]) < 60),
              "abs_lat_ge_60": fit_ok & (np.abs(J["lat"]) >= 60)}
    s["h0_fit"] = {}
    for name, m in groups.items():
        if m.sum() < 10:
            continue
        s["h0_fit"][name] = {
            "n": int(m.sum()),
            "fit_over_iri_law": [float(np.percentile(r_i[m], q)) for q in (10, 50, 90)],
            "fit_over_nequickg_law": [float(np.percentile(r_p[m], q)) for q in (10, 50, 90)],
            "fraction_closer_to_iri_law": float(np.mean(np.abs(np.log(r_i[m])) < np.abs(np.log(r_p[m])))),
            "h0_fit_km_median": float(np.median(J["h0_fit"][m])),
            "kb2bot_km_median": float(np.median((k_o * bb_o)[m])),
        }
    aux = dict(J=J, fn_py=fn_py, fn_pyc=fn_pyc, k_o=k_o, bb_o=bb_o, h0o_i=h0o_i, h0o_p=h0o_p,
               fit_ok=fit_ok, ne_oi=ne_oi, ne_op=ne_op)
    return s, aux


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else round(float(o), 5)
    if isinstance(o, (np.integer,)):
        return int(o)
    return o


def _save(fig, fd: Path, stem: str) -> None:
    fig.savefig(fd / f"{stem}.pdf")
    fig.savefig(fd / f"{stem}.png", dpi=150)


def make_figures(out: Path, s: dict, P: dict, Q: dict, aux: dict, iaux: dict | None) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fd = out / "figures"
    fd.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                         "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7})
    alt = P["alt_km"]
    nla, nlo = LATS.size, LONS.size
    f107s = [int(f) for f in P["f107"]]
    colours = {70: "#1b6ca8", 120: "#e07b00", 180: "#b2182b"}

    # 1. maps: ratio and plasma frequency difference at 900 km
    iz = int(np.argmin(np.abs(alt - 900.0)))
    ratio = P["ne"][..., iz] / Q["ne"][..., iz]
    dfp = fp_mhz(P["ne"][..., iz]) - fp_mhz(Q["ne"][..., iz])
    lat_e = np.r_[LATS - 5.0, LATS[-1] + 5.0]
    lon_e = np.r_[LONS - 15.0, LONS[-1] + 15.0]
    modip = P["modip"].reshape(nla, nlo)
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 3.9), sharex=True, sharey=True,
                            constrained_layout=True)
    for j, f in enumerate(f107s):
        r = np.median(ratio[:, j].reshape(-1, nla, nlo), axis=0)
        d = np.median(dfp[:, j].reshape(-1, nla, nlo), axis=0)
        m1 = axs[0, j].pcolormesh(lon_e, lat_e, r, cmap="viridis", vmin=0.35, vmax=0.55, shading="flat")
        m2 = axs[1, j].pcolormesh(lon_e, lat_e, d, cmap="magma", vmin=-1.0, vmax=0.0, shading="flat")
        for ax in (axs[0, j], axs[1, j]):
            ax.contour(LONS, LATS, modip, levels=[-45, 0, 45], colors="w",
                            linewidths=[0.6, 1.0, 0.6], linestyles=["--", "-", "--"])
            ax.set_xlim(-195, 165)
            ax.set_ylim(-85, 85)
        axs[0, j].set_title(f"F10.7 = {f} sfu")
        axs[1, j].set_xlabel("longitude (deg)")
    axs[0, 0].set_ylabel("latitude (deg)")
    axs[1, 0].set_ylabel("latitude (deg)")
    fig.colorbar(m1, ax=axs[0, :], shrink=0.9, label="Ne ratio, PyIRI / IRI-2020")
    fig.colorbar(m2, ax=axs[1, :], shrink=0.9, label="fp difference (MHz)")
    _save(fig, fd, "topside_maps")
    plt.close(fig)

    # 2. ratio against altitude, and the decomposition
    prof = s["_ratio_profile"]
    hk = np.array(sorted(float(k) for k in prof))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.2), constrained_layout=True,
                                 gridspec_kw={"width_ratios": [1.0, 1.2]})
    for f in f107s:
        med = [prof[str(int(h))][str(f)]["median"] for h in hk]
        lo = [prof[str(int(h))][str(f)]["p10"] for h in hk]
        hi = [prof[str(int(h))][str(f)]["p90"] for h in hk]
        cor = [prof[str(int(h))]["corrected"][str(f)]["median"] for h in hk]
        a1.fill_betweenx(hk, lo, hi, color=colours[f], alpha=0.18, lw=0)
        a1.plot(med, hk, color=colours[f], lw=1.4, label=f"PyIRI / IRI-2020, F10.7 {f}")
        a1.plot(cor, hk, color=colours[f], lw=1.0, ls="--")
    a1.plot([], [], color="0.3", ls="--", lw=1.0, label="PyIRI with the IRI-2020 H0 law")
    a1.axvline(1.0, color="0.5", lw=0.6)
    a1.set_xlim(0.3, 1.15)
    a1.set_ylim(400, 2000)
    a1.set_xlabel("Ne ratio")
    a1.set_ylabel("altitude (km)")
    a1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False, fontsize=6)
    a1.set_title("(a) PyIRI / IRI-2020, median and 10th-90th percentile")
    names = [("peak_density", "NmF2"), ("hmf2_shift", "hmF2"),
             ("thickness_inputs", "H0 inputs"), ("thickness_formula", "H0 law")]
    hs_ = [500, 700, 900]
    width = 0.17
    for i, (key, lab) in enumerate(names):
        vals = [s["altitudes"][str(h)]["terms_gm_by_f107"]["120"][key] for h in hs_]
        a2.bar(np.arange(3) + (i - 2) * width, vals, width, label=lab)
    tot = [s["altitudes"][str(h)]["terms_gm_by_f107"]["120"] for h in hs_]
    tot = [np.prod(list(t.values())) for t in tot]
    a2.bar(np.arange(3) + 2 * width, tot, width, color="k", label="product")
    a2.axhline(1.0, color="0.5", lw=0.6)
    a2.set_xticks(np.arange(3))
    a2.set_xticklabels([f"{h} km" for h in hs_])
    a2.set_ylim(0.3, 1.1)
    a2.set_ylabel("factor in PyIRI / IRI-2020 (geometric mean)")
    a2.legend(ncol=3, frameon=False, loc="upper center")
    a2.set_title("(b) what the ratio is made of, F10.7 = 120")
    _save(fig, fd, "topside_ratio_altitude")
    plt.close(fig)

    # 3. the H0 law, with the ISIS-1 profile fits
    fig, ax = plt.subplots(figsize=(4.6, 3.4), constrained_layout=True)
    x = np.linspace(30.0, 220.0, 200)
    ax.plot(x, x, color="#1b6ca8", lw=1.4, label="IRI-2020 and NeQuick 2: H0 = k B2bot")
    ax.plot(x, x / v_factor(x), color="#b2182b", lw=1.4,
            label="PyIRI 0.1.7: H0 = k B2bot / v (NeQuick G step)")
    kb = aux["k_p"] * aux["bb_p"]
    ax.axvspan(wpercentile(kb, aux["w"], 5), wpercentile(kb, aux["w"], 95), color="0.85",
               zorder=0, label="k B2bot in the campaign, 5th-95th percentile")
    if iaux is not None:
        ok = iaux["fit_ok"]
        xo = (iaux["k_o"] * iaux["bb_o"])[ok]
        yo = iaux["J"]["h0_fit"][ok]
        ax.scatter(xo, yo, s=2, color="k", alpha=0.15, lw=0, label="ISIS-1, H0 fitted to each profile")
        bins = np.arange(40.0, 110.0, 5.0)
        mid, med = [], []
        for b0, b1 in zip(bins[:-1], bins[1:], strict=True):
            m = (xo >= b0) & (xo < b1)
            if m.sum() >= 30:
                mid.append(0.5 * (b0 + b1))
                med.append(np.median(yo[m]))
        ax.plot(mid, med, "o-", color="k", ms=3, lw=1.0, label="ISIS-1, binned median")
    ax.set_xlim(30, 160)
    ax.set_ylim(0, 160)
    ax.set_xlabel("k B2bot (km)")
    ax.set_ylabel("topside thickness H0 (km)")
    ax.legend(frameon=False, loc="upper left")
    _save(fig, fd, "topside_h0_law")
    plt.close(fig)

    # 4. both against ISIS-1, solar maximum
    if iaux is not None:
        b = s["isis"]["bands"]["solar_max"]
        keys = list(b)
        xc = np.arange(len(keys))
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
        series = [("iri2020", "IRI-2020 (PHaRLAP)", "#1b6ca8", "o", "-"),
                  ("pyiri", "PyIRI 0.1.7", "#b2182b", "s", "-"),
                  ("obs_peak_iri_h0", "measured peak, IRI-2020 H0 law", "#1b6ca8", "o", ":"),
                  ("obs_peak_nequickg_h0", "measured peak, PyIRI H0 law", "#b2182b", "s", ":")]
        for key, lab, c, mk, ls in series:
            a1.plot(xc, [b[k]["median_db"][key] for k in keys], marker=mk, color=c, ls=ls,
                    ms=4, lw=1.1, label=lab)
        a1.axhline(0.0, color="0.5", lw=0.6)
        a1.set_xticks(xc)
        a1.set_xticklabels(keys)
        a1.set_xlabel("altitude band (km)")
        a1.set_ylabel("median model / measured, dB")
        a1.set_title("(a) bias at the spacecraft")
        a1.legend(frameon=False, loc="center left", fontsize=6)
        for key, lab, c, mk in [("observed", "ISIS-1 measured", "k", "D"),
                                ("iri2020", "IRI-2020 (PHaRLAP)", "#1b6ca8", "o"),
                                ("pyiri", "PyIRI 0.1.7", "#b2182b", "s")]:
            a2.plot(xc, [b[k]["p05_mhz"][key] for k in keys], marker=mk, color=c, ms=4, lw=1.1, label=lab)
        a2.set_xticks(xc)
        a2.set_xticklabels(keys)
        a2.set_xlabel("altitude band (km)")
        a2.set_ylabel("5th percentile of fp at the spacecraft (MHz)")
        a2.set_title("(b) the lower tail")
        a2.legend(frameon=False)
        _save(fig, fd, "topside_isis")
        plt.close(fig)


def cmd_analyse(out: Path) -> None:
    s, P, Q, aux = grid_analysis(out)
    if (out / "isis_pyiri.npz").is_file():
        s["isis"], iaux = isis_analysis(out)
    else:
        iaux = None
    (out / "summary.json").write_text(json.dumps(_jsonable(s), indent=1))
    make_figures(out, s, P, Q, aux, iaux)


# ------------------------------------------------------------------- main
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("grid")
    g.add_argument("out_dir", type=Path)
    g.add_argument("--matlab", default="matlab")
    g.add_argument("--skip-iri", action="store_true")
    i = sub.add_parser("isis")
    i.add_argument("out_dir", type=Path)
    i.add_argument("--crc", type=Path, required=True)
    i.add_argument("--t1", type=Path, required=True)
    a = sub.add_parser("analyse")
    a.add_argument("out_dir", type=Path)
    args = ap.parse_args(argv)
    if args.cmd == "grid":
        cmd_grid(args.out_dir, args.matlab, args.skip_iri)
    elif args.cmd == "isis":
        cmd_isis(args.out_dir, args.crc, args.t1)
    else:
        cmd_analyse(args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
