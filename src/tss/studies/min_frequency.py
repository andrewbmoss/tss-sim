"""The lower edge of the sounding band, and the profile completeness behind it.

Layer:        study
WP:           WP0, WP1
Claims:       C1, C2
Status:       implemented
Requirements: INS-01a, SCI-02, MIS-01

A sounder cannot radiate below the local plasma frequency at its own position,
so the lower edge of the band is not a design choice.  It is a statistical
property of the ionosphere at the orbit altitude, and it is the parameter that
drives antenna size and therefore the whole small-platform argument.  Alouette I
began at 0.45 MHz, the Radio Plasma Imager at 3 kHz and the TOPside Automated
Sounder at 0.5 MHz, all round numbers with no supporting distribution behind
them.  This study derives the number instead.

PROVENANCE
----------
This is the band-edge harness of ``Dissertation/MinFreqStudy`` (SRM Open
Item 1), ported on 2026-10-19 so that WP0's numbers are produced, tested and
recorded inside this repository.  ``profile_parameters``, ``start_altitude``,
``completeness``, ``run_case`` and ``summarize`` are that harness's functions,
carried across without change to their arithmetic; the orbit sampling is
:class:`tss.platform.orbits.CircularOrbit`, which draws the same samples for the
same seed.  The run that set INS-01a in Rev G (1.2 million IRI-2020 samples
over 60 cases, 2026-09-07) is kept in ``data/reference/band_edge`` and a test
holds this module to it.  The placeholder that carried this study's name before
the port sampled a coarse grid at 400 to 800 km and, because it built the IRI
tier without a solar index, silently ran on the analytic stand-ins only.

WHAT A RUN IS
-------------
A *plan*, a YAML file under ``config/studies``, lists the cases: altitudes,
orbit classes (O1 80 deg, O2 30 deg), solar cases, topside models, samples per
case, seed and the altitude grid.  ``band_edge_reference`` reproduces the Rev G
run.  A solar case is either a user R12 (the harness's route, storm model off)
or ``mode: epoch``, which drives IRI from its own index files over a stated
window with the storm switches the plan gives.  Select a plan and override any
key from the command line::

    tss study min_frequency --set min_frequency.plan=band_edge_reference
    tss study min_frequency --set min_frequency.samples=2000

Each case writes one row of ``minfreq_summary.csv`` (the harness's columns,
then the additions listed in :data:`EXTRA_SUMMARY_COLUMNS`), its completeness
curve, and an ``.npz`` of per-sample scalars (time, position, local time, fN at
the spacecraft, NmF2, hmF2 and the IRI extras with the indices it ran on) for
the SCI-02 and solar-minimum analyses.  Profiles are not kept.

THE COMBINATION DIRECTION IS THE WHOLE RESULT
---------------------------------------------
The per-model fifth percentiles must be combined by taking the MINIMUM, not the
maximum.  Taking the maximum picks the most optimistic model at each point and
produces 1.9 MHz, which delivered 53.4 per cent profile coverage against the 95
per cent claimed.  The minimum is the conservative envelope across a
disagreeing ensemble and produces 0.33 MHz.  An earlier revision of this work
took the maximum, passed every internal check, and was caught only by the ISIS-1
archive comparison.  That is the evidence for the second governing principle of
Section 3.1, and :func:`combine_ensemble` carries a test that pins the
direction so the error cannot come back.

THE TWO METRICS
---------------
Coverage.  If f_min <= fN_sat the whole topside profile is reachable, because
the sweep effectively starts at fN_sat.  If f_min > fN_sat the region between
the spacecraft and the altitude where fN = f_min is lost.  So coverage is
P(fN_sat >= f_min), and 95 per cent coverage is the 5th PERCENTILE of fN_sat.

Completeness.  Coverage scores a profile truncated by 1 km the same as one
truncated by 300 km.  Completeness measures what survives::

    altitude completeness = (h(f_start) - hmF2) / (h_sat - hmF2)
    density  completeness = ln(foF2 / f_start) / ln(foF2 / fN_sat)

with f_start = max(f_min, fN_sat).  The altitude form is the one SCI-02's goal
is scored by.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import shutil
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyInputMissing, StudyResult, register_study

#: Coverage fractions INS-01a is stated at.  Threshold then goal.
COVERAGE_LEVELS = (0.95, 0.99)

#: The harness's summary columns, in its order, so a summary can be compared
#: with ``data/reference/band_edge/minfreq_summary.csv`` column for column.
SUMMARY_COLUMNS = (
    "orbit_class", "inclination_deg", "altitude_km", "solar", "r12", "f107", "topside",
    "backend", "provenance", "f_min_for_90pct_coverage_MHz", "f_min_for_95pct_coverage_MHz",
    "f_min_for_99pct_coverage_MHz", "samples_used", "samples_dropped", "fn_sat_median_MHz",
    "fn_sat_p05_MHz", "fn_sat_p50_MHz", "fn_sat_p95_MHz", "fof2_median_MHz", "hmf2_median_km",
    "band_ratio_median",
)

#: Columns this port adds after the harness's, for the closure analyses.
EXTRA_SUMMARY_COLUMNS = (
    "solar_mode", "epoch_start", "duration_days", "fn_sat_p01_MHz", "hmf2_p001_km",
    "hmf2_p999_km", "frac_hsat_minus_hmf2_lt_100km", "samples_file",
)

#: The default plan, which reproduces the Rev G band-edge run.
DEFAULT_PLAN = "band_edge_reference"


# ---------------------------------------------------------------------------
# The metrics (MinFreqStudy/minfreq_study.py, unchanged in arithmetic)
# ---------------------------------------------------------------------------
def plasma_frequency(ne_m3: np.ndarray | float) -> np.ndarray:
    """Plasma frequency, Hz, from electron density, m^-3."""
    from tss.constants import plasma_frequency as _fp

    return np.asarray(_fp(np.asarray(ne_m3, dtype=float)), dtype=float)


def combine_ensemble(per_model_percentiles: np.ndarray) -> float:
    """Combine per-model percentiles into one band edge.

    Takes the minimum across models, which is the conservative envelope.  Do not
    change this to a maximum.  See the module docstring: the maximum gave
    1.9 MHz and 53.4 per cent real coverage against 95 per cent claimed.
    """
    values = np.asarray(per_model_percentiles, dtype=float)
    if values.size == 0:
        raise ValueError("cannot combine an empty ensemble")
    return float(np.nanmin(values))


def coverage_fraction(fn_sat_hz: np.ndarray, f_min_hz: float) -> float:
    """Fraction of samples for which the FULL topside profile is reachable."""
    fn = np.asarray(fn_sat_hz, dtype=float)
    fn = fn[np.isfinite(fn)]
    if fn.size == 0:
        return float("nan")
    return float(np.mean(fn >= f_min_hz))


def f_min_for_coverage(fn_sat_hz: np.ndarray, coverage: float) -> float:
    """Lowest start frequency that reaches ``coverage`` of samples.

    Note the inversion: 95 per cent coverage is the 5th percentile of fN_sat.
    """
    fn = np.asarray(fn_sat_hz, dtype=float)
    fn = fn[np.isfinite(fn)]
    if fn.size == 0:
        return float("nan")
    return float(np.percentile(fn, 100.0 * (1.0 - coverage)))


def profile_parameters(ne, alt_grid_m, h_sat, extra=None):
    """Reduce each profile to fN at the spacecraft, foF2, hmF2 and NmF2.

    ``ne`` is [n_samples, n_alt] on the shared ``alt_grid_m``.  The peak is
    taken over the grid below the spacecraft, which is where a topside sounder
    can see, and IRI's own NmF2 and hmF2 replace it where the backend supplies
    them (``extra`` columns 0 and 1), because they are not quantised to the grid
    step.
    """
    alt = np.asarray(alt_grid_m, dtype=float)
    below = alt <= h_sat
    if not below.any():
        raise ValueError("altitude grid does not reach the spacecraft altitude")

    i_sat = int(np.argmin(np.abs(alt - h_sat)))
    ne_sat = ne[:, i_sat]

    # IRI returns -1 where it cannot produce a value; the backend maps that to
    # NaN.  nanargmax would raise on an all-NaN row, so fill first.
    sub = np.where(np.isfinite(ne[:, below]), ne[:, below], -np.inf)
    k = np.argmax(sub, axis=1)
    nmf2 = sub[np.arange(sub.shape[0]), k]
    hmf2 = alt[below][k]
    bad = ~np.isfinite(nmf2)
    nmf2 = np.where(bad, np.nan, nmf2)
    hmf2 = np.where(bad, np.nan, hmf2)

    if extra is not None:
        extra = np.asarray(extra)
        if extra.shape[0] != ne.shape[0]:
            raise ValueError(
                f"extra has {extra.shape[0]} rows but ne has {ne.shape[0]}; "
                f"peak parameters and profiles are indexed together")
        nmf2_iri = extra[:, 0]
        hmf2_iri = extra[:, 1] * 1e3
        use = np.isfinite(nmf2_iri) & np.isfinite(hmf2_iri)
        nmf2 = np.where(use, nmf2_iri, nmf2)
        hmf2 = np.where(use, hmf2_iri, hmf2)

    return dict(fn_sat=plasma_frequency(ne_sat), ne_sat=ne_sat,
                fof2=plasma_frequency(nmf2), nmf2=nmf2, hmf2=hmf2, i_sat=i_sat)


def start_altitude(ne, alt_grid_m, i_sat, f_start):
    """Altitude at which the plasma frequency first equals ``f_start``.

    Searches downward from the spacecraft; this is where the profile actually
    begins when f_min truncates it.  Linear interpolation in plasma frequency
    between grid points.
    """
    alt = np.asarray(alt_grid_m, dtype=float)
    fn = plasma_frequency(ne)
    n = ne.shape[0]
    out = np.full(n, np.nan)
    for r in range(n):
        col = fn[r, : i_sat + 1][::-1]
        h = alt[: i_sat + 1][::-1]
        idx = np.nonzero(col >= f_start[r])[0]
        if idx.size == 0:
            continue
        j = idx[0]
        if j == 0:
            out[r] = h[0]
            continue
        f0, f1 = col[j - 1], col[j]
        h0, h1 = h[j - 1], h[j]
        w = 0.0 if f1 == f0 else (f_start[r] - f0) / (f1 - f0)
        out[r] = h0 + w * (h1 - h0)
    return out


def completeness(f_min, params, ne, alt_grid_m, h_sat):
    """Altitude and density completeness of the retrievable profile at ``f_min``."""
    fn_sat = params["fn_sat"]
    fof2 = params["fof2"]
    hmf2 = params["hmf2"]

    f_start = np.maximum(f_min, fn_sat)
    usable = f_start < fof2

    h_start = start_altitude(ne, alt_grid_m, params["i_sat"], f_start)

    span = np.maximum(h_sat - hmf2, 1.0)
    alt_c = np.clip((h_start - hmf2) / span, 0.0, 1.0)
    alt_c = np.where(usable & np.isfinite(h_start), alt_c, 0.0)

    denom = np.log(np.maximum(fof2 / np.maximum(fn_sat, 1.0), 1.0 + 1e-12))
    den_c = np.log(np.maximum(fof2 / np.maximum(f_start, 1.0), 1.0)) / denom
    den_c = np.clip(np.where(usable, den_c, 0.0), 0.0, 1.0)

    return dict(f_start=f_start, h_start=h_start, usable=usable,
                altitude_completeness=alt_c, density_completeness=den_c)


def run_case(backend, orbit, epoch_start, n_samples, alt_grid_m,
             duration_days=365.25, seed=0, chunk=20000):
    """Sample one orbit and evaluate the model along it.

    ``chunk`` bounds the points handed to the backend in one call.  Keep it
    large: the MATLAB backend starts a fresh ``matlab -batch`` per call, and
    start-up is about 10 s against about 2.3 ms of model time per sample.

    A backend that reports per-point extras overwrites ``backend.extra`` on
    every call, so it holds the LAST chunk only.  It is collected chunk by
    chunk here; reading it once after the loop would pair one chunk of peak
    parameters with every profile, which the harness once did.
    """
    s = orbit.sample(n_samples, duration_days=duration_days, seed=seed)
    when = [epoch_start + dt.timedelta(seconds=float(t)) for t in s["t_seconds"]]

    ne = np.empty((int(n_samples), len(alt_grid_m)))
    extra_parts, n_chunks = [], 0
    for a in range(0, int(n_samples), chunk):
        b = min(a + chunk, int(n_samples))
        n_chunks += 1
        ne[a:b] = backend.profiles(s["lat_deg"][a:b], s["lon_deg"][a:b],
                                   when[a:b], alt_grid_m)
        ex = getattr(backend, "extra", None)
        if ex is not None:
            ex = np.asarray(ex)
            if ex.shape[0] != b - a:
                raise RuntimeError(
                    f"backend returned {ex.shape[0]} extra rows for a chunk of "
                    f"{b - a} profiles; the two are indexed together and must agree")
            extra_parts.append(ex)

    if not extra_parts:
        extra = None
    elif len(extra_parts) != n_chunks:
        raise RuntimeError(
            f"backend supplied extra for {len(extra_parts)} of {n_chunks} chunks; "
            "it must be all or none, or the rows do not align")
    else:
        extra = np.vstack(extra_parts)

    params = profile_parameters(ne, alt_grid_m, orbit.altitude, extra=extra)
    params["extra"] = extra
    return s, ne, params


def summarize(params, ne, alt_grid_m, h_sat, coverages=(0.90, 0.95, 0.99), f_min_grid=None):
    """Percentile table and completeness curve for one case."""
    fn = np.asarray(params["fn_sat"], dtype=float)
    good = np.isfinite(fn) & np.isfinite(params["fof2"])
    n_bad = int((~good).sum())
    fn = fn[good]
    if fn.size == 0:
        raise RuntimeError("every sample failed -- check the IRI binding")
    rows: dict[str, Any] = {f"f_min_for_{int(c * 100)}pct_coverage_MHz":
                            f_min_for_coverage(fn, c) / 1e6 for c in coverages}
    rows["samples_used"] = int(fn.size)
    rows["samples_dropped"] = n_bad
    rows.update({
        "fn_sat_median_MHz": float(np.median(fn)) / 1e6,
        "fn_sat_p05_MHz": float(np.percentile(fn, 5)) / 1e6,
        "fn_sat_p50_MHz": float(np.percentile(fn, 50)) / 1e6,
        "fn_sat_p95_MHz": float(np.percentile(fn, 95)) / 1e6,
        "fof2_median_MHz": float(np.nanmedian(params["fof2"])) / 1e6,
        "hmf2_median_km": float(np.nanmedian(params["hmf2"])) / 1e3,
        "band_ratio_median": float(np.nanmedian(params["fof2"][good] / np.maximum(fn, 1.0))),
    })

    if f_min_grid is None:
        f_min_grid = np.logspace(np.log10(0.1e6), np.log10(4e6), 40)
    curve = []
    for fm in f_min_grid:
        c = completeness(fm, params, ne, alt_grid_m, h_sat)
        curve.append(dict(
            f_min_MHz=fm / 1e6,
            coverage=coverage_fraction(fn, fm),
            mean_altitude_completeness=float(np.mean(c["altitude_completeness"])),
            mean_density_completeness=float(np.mean(c["density_completeness"])),
            p05_altitude_completeness=float(np.percentile(c["altitude_completeness"], 5)),
        ))
    return rows, curve


# ---------------------------------------------------------------------------
# Plans and backends
# ---------------------------------------------------------------------------
def load_plan(config: dict[str, Any]) -> dict[str, Any]:
    """The plan named by ``min_frequency.plan``, with any other key overriding it."""
    from tss.config import load_yaml, resolve_path

    opts = dict(config.get("min_frequency", {}) or {})
    name = str(opts.pop("plan", DEFAULT_PLAN))
    plan = dict(load_yaml(resolve_path(name, subdir="studies")))
    plan.update(opts)
    plan.setdefault("name", name)
    for key in ("altitudes_km", "orbit_classes", "topsides", "solar", "samples", "alt_grid_km"):
        if key not in plan:
            raise ValueError(f"band-edge plan {name!r} has no {key!r}")
    lo, hi, step = (float(v) for v in plan["alt_grid_km"])
    if hi < max(float(a) for a in plan["altitudes_km"]):
        raise ValueError(f"plan {name!r}: the altitude grid stops at {hi:g} km, below the "
                         "highest spacecraft altitude")
    return plan


def _when(text: str) -> dt.datetime:
    return dt.datetime.fromisoformat(str(text).replace("Z", "")).replace(tzinfo=None)


def _slug(text: str) -> str:
    return str(text).replace(" ", "-").replace(".", "").lower()


def build_backend(plan: dict[str, Any], topside: str, solar: dict[str, Any],
                  workdir: Path, probe: bool):
    """The L0 tier for one topside and one solar case."""
    from tss.l0_ionosphere.base import r12_to_f107
    from tss.registry import build

    if str(plan.get("backend", "pharlap")) == "analytic":
        if solar.get("mode", "r12") != "r12":
            raise ValueError("the analytic stand-in takes only a user R12")
        return build("l0", "analytic", f107=float(r12_to_f107(float(solar["r12"]))))
    kw: dict[str, Any] = dict(
        model="pharlap", topside=topside, workdir=str(workdir),
        timeout=float(plan.get("timeout_s", 14400)), study="band_edge", progress=0,
        matlab=str(plan.get("matlab", "matlab")), probe=probe,
    )
    if plan.get("pharlap_root"):
        kw["pharlap_root"] = str(plan["pharlap_root"])
    storm = dict(solar.get("iri_options", {}) or {})
    if storm:
        kw["iri_options"] = storm
    if solar.get("mode", "r12") == "epoch":
        kw["epoch_indices"] = True
    else:
        kw["r12"] = float(solar["r12"])
    return build("l0", "iri", **kw)


def _cases(plan: dict[str, Any]):
    """(topside, solar) pairs in the harness's loop order."""
    for topside in plan["topsides"]:
        for solar in plan["solar"]:
            yield topside, dict(solar)


# ---------------------------------------------------------------------------
# The study
# ---------------------------------------------------------------------------
@register_study(
    "min_frequency",
    description=("Lower band edge and profile completeness from IRI-2020 sampled along the "
                 "orbit (the MinFreqStudy harness)"),
    work_packages=["WP0", "WP1"],
    claims=["C1", "C2"],
    status="implemented",
    requirements=["INS-01a", "SCI-02", "MIS-01"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    from tss.l0_ionosphere.base import r12_to_f107
    from tss.platform.orbits import orbit_classes

    plan = load_plan(config)
    lo, hi, step = (float(v) for v in plan["alt_grid_km"])
    alt_grid = np.arange(lo * 1e3, hi * 1e3 + 1.0, step * 1e3)
    n = int(plan["samples"])
    seed = int(plan.get("seed", 17))
    chunk = int(plan.get("chunk", 20000))
    keep_samples = bool(plan.get("keep_samples", True))
    work = out_dir / "iri_work"
    sample_dir = out_dir / "samples"
    sample_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    curves: list[list[Any]] = []
    probed = False
    for topside, solar in _cases(plan):
        mode = str(solar.get("mode", "r12"))
        epoch0 = _when(solar.get("epoch_start", plan.get("epoch_start", "2015-01-01T00:00:00")))
        days = float(solar.get("duration_days", plan.get("duration_days", 365.25)))
        try:
            backend = build_backend(plan, topside, solar, work / f"{_slug(topside)}_{solar['name']}",
                                    probe=not probed)
        except (ImportError, RuntimeError, FileNotFoundError) as exc:
            raise StudyInputMissing(f"the IRI backend could not be built: {exc}") from exc
        probed = True
        r12 = float(solar["r12"]) if mode == "r12" else float("nan")
        for alt_km in plan["altitudes_km"]:
            h_sat = float(alt_km) * 1e3
            for cls, orb in orbit_classes(h_sat).items():
                if cls not in plan["orbit_classes"]:
                    continue
                s, ne, p = run_case(backend, orb, epoch0, n, alt_grid, duration_days=days,
                                    seed=seed, chunk=chunk)
                summ, curve = summarize(p, ne, alt_grid, h_sat)
                hm = p["hmf2"]
                fn = p["fn_sat"]
                fname = f"{cls}_{int(alt_km)}km_{solar['name']}_{_slug(topside)}.npz"
                if keep_samples:
                    ex = p.get("extra")
                    np.savez_compressed(
                        sample_dir / fname,
                        t_seconds=s["t_seconds"], lat_deg=s["lat_deg"].astype(np.float32),
                        lon_deg=s["lon_deg"].astype(np.float32),
                        local_time_hours=s["local_time_hours"].astype(np.float32),
                        fn_sat_hz=fn.astype(np.float32), ne_sat=p["ne_sat"].astype(np.float32),
                        nmf2=p["nmf2"].astype(np.float32), hmf2_m=hm.astype(np.float32),
                        extra=(np.asarray(ex, dtype=np.float32) if ex is not None
                               else np.zeros((n, 0), np.float32)),
                        extra_columns=np.array(_extra_columns()),
                        epoch_start=str(epoch0), h_sat_m=h_sat,
                    )
                rec = dict(
                    orbit_class=cls, inclination_deg=orb.inclination_deg, altitude_km=h_sat / 1e3,
                    solar=solar["name"], r12=r12,
                    f107=round(float(r12_to_f107(r12)), 1) if mode == "r12" else float("nan"),
                    topside=topside, backend=backend.label, provenance=backend.provenance,
                    **summ,
                    solar_mode=mode, epoch_start=epoch0.isoformat(), duration_days=days,
                    fn_sat_p01_MHz=float(np.nanpercentile(fn, 1)) / 1e6,
                    hmf2_p001_km=float(np.nanpercentile(hm, 0.1)) / 1e3,
                    hmf2_p999_km=float(np.nanpercentile(hm, 99.9)) / 1e3,
                    frac_hsat_minus_hmf2_lt_100km=float(np.nanmean((h_sat - hm) < 100e3)),
                    samples_file=f"samples/{fname}" if keep_samples else "",
                )
                rows.append(rec)
                for c in curve:
                    curves.append([cls, h_sat / 1e3, solar["name"], topside, f"{c['f_min_MHz']:.4f}",
                                   f"{c['coverage']:.4f}", f"{c['mean_altitude_completeness']:.4f}",
                                   f"{c['mean_density_completeness']:.4f}",
                                   f"{c['p05_altitude_completeness']:.4f}"])
                print(f"{cls} {h_sat / 1e3:5.0f} km  solar {solar['name']:10s} {topside:20s} "
                      f"f_min(95%) {summ['f_min_for_95pct_coverage_MHz']:6.3f} MHz", flush=True)
    if bool(plan.get("clean_workdir", True)):
        shutil.rmtree(work, ignore_errors=True)

    summary_csv = out_dir / "minfreq_summary.csv"
    with open(summary_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(SUMMARY_COLUMNS) + list(EXTRA_SUMMARY_COLUMNS))
        w.writeheader()
        w.writerows(rows)
    curves_csv = out_dir / "completeness_curves.csv"
    with open(curves_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["orbit_class", "altitude_km", "solar", "topside", "f_min_MHz", "coverage",
                    "mean_alt_completeness", "mean_density_completeness",
                    "p05_alt_completeness"])
        w.writerows(curves)

    ens = ensemble_table(rows)
    ens_csv = out_dir / "ensemble.csv"
    with open(ens_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ens[0].keys()))
        w.writeheader()
        w.writerows(ens)
    gov = min(ens, key=lambda r: r["f_min_95_ensemble_MHz"])
    plan_json = out_dir / "plan.json"
    plan_json.write_text(json.dumps(plan, indent=1, default=str), encoding="utf-8")

    metrics = {
        "plan": plan["name"],
        "n_cases": len(rows),
        "samples_per_case": n,
        "combination": "minimum of per-model 5th percentiles",
        "governing_case": {k: gov[k] for k in ("orbit_class", "altitude_km", "solar")},
        "f_min_mhz_at_95pc": gov["f_min_95_ensemble_MHz"],
        "f_min_mhz_at_99pc": gov["f_min_99_ensemble_MHz"],
        "unsafe_maximum_combination_mhz": gov["f_min_95_unsafe_max_MHz"],
        "requirement": "INS-01a",
    }
    summary = (
        f"Plan {plan['name']}: {len(rows)} cases of {n} samples, sampled along the orbit.\n\n"
        f"The lowest ensemble band edge is {gov['f_min_95_ensemble_MHz']:.3f} MHz at 95 per cent "
        f"coverage ({gov['f_min_99_ensemble_MHz']:.3f} MHz at 99 per cent), in "
        f"{gov['orbit_class']} at {gov['altitude_km']:.0f} km, solar {gov['solar']}, taking the "
        "minimum of the per-model 5th percentiles.  Taking the maximum would give "
        f"{gov['f_min_95_unsafe_max_MHz']:.3f} MHz there, which is the error an earlier revision "
        "made.\n\nPer-case rows are in minfreq_summary.csv, the ensemble in ensemble.csv, "
        "completeness curves in completeness_curves.csv and per-sample scalars under samples/."
    )
    outputs = [summary_csv, curves_csv, ens_csv, plan_json]
    return StudyResult(outputs=outputs, summary=summary, metrics=metrics)


def _extra_columns() -> tuple[str, ...]:
    from tss.l0_ionosphere.backends import EXTRA_COLUMNS

    return EXTRA_COLUMNS


def ensemble_table(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Per (orbit class, altitude, solar case): the minimum across topsides."""
    keys: dict[tuple, list[dict[str, Any]]] = {}
    for r in rows:
        keys.setdefault((r["orbit_class"], r["altitude_km"], r["solar"]), []).append(r)
    out = []
    for (cls, alt, solar), rs in keys.items():
        f95 = [r["f_min_for_95pct_coverage_MHz"] for r in rs]
        f99 = [r["f_min_for_99pct_coverage_MHz"] for r in rs]
        out.append(dict(
            orbit_class=cls, altitude_km=alt, solar=solar, n_models=len(rs),
            f_min_95_ensemble_MHz=combine_ensemble(np.array(f95)),
            f_min_99_ensemble_MHz=combine_ensemble(np.array(f99)),
            f_min_95_unsafe_max_MHz=float(np.nanmax(f95)),
            governing_model_95=rs[int(np.nanargmin(f95))]["topside"],
        ))
    return out
