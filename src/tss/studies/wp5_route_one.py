"""Layer: studies
WP: WP5
Claims: C3
Status: implemented
Requirements: SCI-15, INS-15, INS-26

WP5 Phase 3: route one at signal level, the study ``wp5_route_one``.

The Phase 0 scenes (:mod:`tss.studies.wp5_campaign`) through L1, L2 and the
drift tier, with noise, fading, the oscillator, leakage and the canceller, and
the diffuse halo.  The layout is the one Andrew approved on 2026-11-27.

THE BLOCKS
----------
``scenes``
    Every Phase 0 scene once, without noise, through a receiver without
    impairments, quantisation or leakage.  Facets, tilted facets, extended and
    long plates, field-aligned sheets and the two budget points.  The
    parameter-level result of ``wp5_campaign`` for the same scene is recorded
    beside it (``param_err_mps``).
``monte_carlo``
    Nine scenes (the untilted facet at nadir, 25 m/s, in each of the seven
    ``c3_verify`` environments, and the two budget points) at 0, 5, 10, 15 and
    20 dB per sweep and at the link budget's own SNR, without fading and with
    Rayleigh specular fading at 0.1, 0.5 and 1.0 s, the full instrument with the
    oscillator at the INS-26 threshold.
``oscillator``
    The same nine scenes without fading at 10 dB and the link budget, with the
    oscillator at the INS-26 goal.
``halo``
    The plasma 5 MHz facet without a halo and with halos of two facet lengths at
    K = 10 and 3, and the patch alone (K = 0), at 10 and 20 dB; a new halo each
    seed.
``orbit``
    The three plasma facets on a 600 km orbit (3-D plates), at 10 and 20 dB,
    processed with a prior whose scale height is 100, 200 (the medium's) and
    300 km.

NOISE AND SNR
-------------
The SNR per sweep is the route-one budget's (decision (l)): the per-sweep value
at the envelope's peak, after the gate, over its noise.  For a stated level the
full instrument runs with its thermal noise off and white noise is added to the
records, set against the same scene's unfaded record with the oscillator at the
threshold (the residual of the impairments is more than 20 dB below the highest
level).  The ``link`` level is the full instrument with its thermal noise on and
nothing added.  ``snr_db`` in every noisy row is measured from the detector's
noise reference of that record against the unfaded peak, so it is the level's
own check and, at ``link``, the drift-mode SNR per sweep (TV assumption 1).

WHAT IS RECORDED
----------------
One row per estimator (``pinned``, ``joint``, ``joint_fading``) and echo; the
facet's echo is the one within two range cells of its beat, and a row is kept
for it when it is not detected.  Every echo carries route one's
``peak_contrast`` and the drift tier's ``far_power_over_noise`` and
``edge_power_over_noise`` (``docs/drift-tier.md`` section 16: recorded, no
rejection rule).  The summary
reports availability and error against level, fading, oscillator, halo
strength and scale height, with the matching Cramer-Rao bound at the nominal
SNR.

PRIORS AND TRUTH
----------------
Straight scenes use the truth profile as the prior; the orbit scenes the true
``f_N`` at the spacecraft with the assumed scale height (decision (p)).  The
onboard ``f_N`` from an ionogram is not part of this study.  Truth is the
facet's vertical drift (a tilted facet's predicted bias is recorded beside it),
for a sheet the line-of-sight velocity of the parameter-level scene, and for
the patch alone the drift it moves with.

NOT DONE
--------
Phase-wander and Rician fading, the range-correlation credit with and without
(open item 12), clutter and sky noise beyond the thermal level, a halo with
fading, a halo that rises and falls with the pass, and Poole's comparator
(run separately, ``scripts/wp5_poole.py``).

RUNNING
-------
``tss study wp5_route_one --out results/wp5_route_one --strict``.  Hours: about
five on 14 processes.  Each job writes ``parts/<job>.json`` as it ends, so a
run that stops is resumed by the same command with the same ``--out`` from the
same commit.  ``progress.json`` holds the count and the time left.  The
``wp5_route_one`` configuration block (``--set wp5_route_one.seeds=30``; lists
as JSON) sets :class:`Settings`' fields.
"""

from __future__ import annotations

import csv
import dataclasses
import gzip
import json
import math
import os
import pickle
import time
import zlib
from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tss.constants import C_LIGHT
from tss.studies import StudyResult, register_study

__all__ = [
    "BLOCKS",
    "ROW_FIELDS",
    "Settings",
    "Job",
    "Setup",
    "settings",
    "jobs",
    "scene_setup",
    "needed_caches",
    "prepare",
    "run_job",
    "summarise",
]

BLOCKS = ("scenes", "monte_carlo", "oscillator", "halo", "orbit")
ESTIMATORS = ("pinned", "joint", "joint_fading")
V_S = 7600.0
REPETITION_HZ = 45.0
SWEEP_S = 1.0 / REPETITION_HZ
PLATFORM_ALT_M = 600e3
EPOCH = "2026-01-01T12:00:00Z"
#: The instrument configuration the campaigns use, on one 600 km platform.
CONFIGURATION = {"name": "leo_600", "orbit": {"altitude_m": 6.0e5, "inclination_deg": 97.4},
                 "antenna_count": 1, "antenna_length_m": 5.0}
#: Two range cells of drift-mode beat, Hz: an echo this close to the facet's beat is its echo.
NEAR_HZ = 132.0
MC_SCENES = (
    "A_free_3MHz_500km_x0_v25", "A_free_5MHz_500km_x0_v25", "A_free_10MHz_500km_x0_v25",
    "A_free_5MHz_1000km_x0_v25", "A_plasma_3MHz_mu0.5_x0_v25", "A_plasma_5MHz_mu0.5_x0_v25",
    "A_plasma_10MHz_mu0.8_x0_v25", "E_plasma_5MHz_P500km_v25", "E_plasma_10MHz_P500km_v25",
)
HALO_SCENE = "A_plasma_5MHz_mu0.5_x0_v25"
ORBIT_ENVIRONMENTS = ("plasma_3MHz_mu0.5", "plasma_5MHz_mu0.5", "plasma_10MHz_mu0.8")
ORBIT_INCLINATION_DEG, ORBIT_U0_DEG, ORBIT_DRIFT_MPS = 60.0, 20.0, 25.0

ROW_FIELDS = (
    "block", "scene_id", "group", "echo_class_truth", "medium", "f_mhz", "level_db", "tau_c_s",
    "oscillator", "halo", "h_assumed_km", "seed", "estimator", "echo_id", "facet_echo", "detected",
    "formed", "echo_class", "truth_mps", "bias_mps", "drift_mps", "err_mps", "param_err_mps",
    "drift_sigma_mps", "drift_pinned_mps", "slope_ratio", "peak_contrast", "far_power_over_noise",
    "edge_power_over_noise", "tau_c_fit_s", "snr_db", "range_m", "blanked_fraction", "n_sweeps", "w_s", "p_group_m", "mu_sat",
    "flags",
)


@dataclass(frozen=True)
class Settings:
    """The study's numbers; the ``wp5_route_one`` configuration block overrides them."""

    seed: int = 20261127
    seeds: int = 100
    processes: int = 14
    blocks: tuple[str, ...] = BLOCKS
    scene_ids: tuple[str, ...] = ()  # the "scenes" block: empty is every Phase 0 scene
    mc_scenes: tuple[str, ...] = MC_SCENES
    levels: tuple[Any, ...] = (0.0, 5.0, 10.0, 15.0, 20.0, "link")
    tau_c_s: tuple[Any, ...] = (None, 0.1, 0.5, 1.0)
    estimators: tuple[str, ...] = ESTIMATORS
    oscillator_main: float = 5.8e-8  # INS-26 threshold
    oscillator_second: float = 3.2e-9  # INS-26 goal
    oscillator_levels: tuple[Any, ...] = (10.0, "link")
    halo_cases: Mapping[str, float | None] = dataclasses.field(
        default_factory=lambda: {"none": None, "K10": 10.0, "K3": 3.0, "K0": 0.0})
    halo_ratio: float = 2.0
    halo_points: int = 200
    halo_levels: tuple[Any, ...] = (10.0, 20.0)
    orbit_environments: tuple[str, ...] = ORBIT_ENVIRONMENTS
    orbit_levels: tuple[Any, ...] = (10.0, 20.0)
    orbit_scale_heights_km: tuple[float, ...] = (100.0, 200.0, 300.0)


def settings(config: Mapping[str, Any]) -> Settings:
    """:class:`Settings` from the study configuration."""
    given = dict(config.get("wp5_route_one", {}) or {})
    known = {f.name for f in dataclasses.fields(Settings)}
    unknown = set(given) - known
    if unknown:
        raise ValueError(f"unknown wp5_route_one keys {sorted(unknown)}")
    for key, value in given.items():
        if isinstance(getattr(Settings(), key), tuple) and isinstance(value, str):
            raise ValueError(f'wp5_route_one.{key} needs a JSON list, such as ["a"]')
        if isinstance(value, list):
            given[key] = tuple(value)
    for key in ("levels", "oscillator_levels", "halo_levels", "orbit_levels"):
        if key in given:
            given[key] = tuple(v if v == "link" else float(v) for v in given[key])
    if "tau_c_s" in given:
        given["tau_c_s"] = tuple(None if v in (None, "none") else float(v) for v in given["tau_c_s"])
    out = Settings(**given)
    bad = set(out.blocks) - set(BLOCKS)
    if bad:
        raise ValueError(f"unknown blocks {sorted(bad)}")
    return out


@dataclass(frozen=True)
class Job:
    """One unit of work: a scene (or orbit environment), a seed, and what varies."""

    block: str
    scene_id: str
    seed: int = 0
    tau_c_s: float | None = None
    halo: str = "none"

    @property
    def job_id(self) -> str:
        tau = "none" if self.tau_c_s is None else f"{self.tau_c_s:g}"
        return f"{self.block}__{self.scene_id}__tau{tau}__halo{self.halo}__s{self.seed}"


def jobs(s: Settings) -> list[Job]:
    """Every job of the study, in block order."""
    from tss.studies.wp5_campaign import SCENES

    out: list[Job] = []
    if "scenes" in s.blocks:
        out += [Job("scenes", i) for i in (s.scene_ids or tuple(sc.scene_id for sc in SCENES))]
    if "monte_carlo" in s.blocks:
        out += [Job("monte_carlo", i, n, tau_c_s=tau)
                for i in s.mc_scenes for tau in s.tau_c_s for n in range(s.seeds)]
    if "oscillator" in s.blocks:
        out += [Job("oscillator", i, n) for i in s.mc_scenes for n in range(s.seeds)]
    if "halo" in s.blocks:
        out += [Job("halo", HALO_SCENE, n, halo=case) for case in s.halo_cases for n in range(s.seeds)]
    if "orbit" in s.blocks:
        out += [Job("orbit", f"O_{env}", n) for env in s.orbit_environments for n in range(s.seeds)]
    return out


# ---------------------------------------------------------------------------
# Scenes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Setup:
    """A scene as the signal-level chain sees it."""

    scene_id: str
    group: str
    echo_class_truth: str
    medium: str
    f_hz: float
    facet: dict[str, Any]
    l1_kwargs: dict[str, Any]
    t_s: np.ndarray  # sweep start times
    truth_mps: float
    bias_mps: float
    w_s: float  # the facet's envelope half-power width, for the bound


def _starts(t_centre: float, half_s: float) -> np.ndarray:
    """Sweep start times whose centres cover ``t_centre +/- half_s``."""
    n = int(math.ceil(half_s * REPETITION_HZ))
    return t_centre + np.arange(-n, n + 1) * SWEEP_S - 0.5 * SWEEP_S


def scene_setup(scene: Any) -> Setup:
    """One ``wp5_campaign`` scene as a facet, a medium and a time window.

    Facets take ``c3_verify``'s geometry and window (centred where the platform
    is on the facet's normal, 2.5 envelope widths each side, or the scene's own);
    sheets take ``field_aligned_verify``'s sheet and window.
    """
    from tss.l1_propagation.plate_echo import Plate
    from tss.studies import c3_verify as cv
    from tss.studies.route_one import envelope_width_s

    if scene.group == "S_sheet":
        from tss.studies import field_aligned_verify as fa
        from tss.studies.wp5_campaign import run_scene

        plasma = scene.medium == "plasma"
        half = 1.0 if plasma else 2.0
        sheet = fa._sheet(scene.dip_deg, scene.rotation_deg, 25.0, 30.0, 1e3)
        th = fa.se.heading_vector(np.radians(scene.heading_deg))
        along = abs(float(np.dot(th, fa.se.sheet_frame(sheet)["b"]))) * fa.V_S * half
        length = max(30e3, 2.0 * (along + 12e3) / 0.7)
        facet = {"id": "s", "kind": "sheet", "dip_deg": scene.dip_deg, "rotation_deg": scene.rotation_deg,
                 "heading_deg": scene.heading_deg, "length_m": length, "width_m": sheet.width_m,
                 "curvature_radius_m": float(sheet.curvature_radius_m), "tukey_fraction": sheet.tukey_fraction,
                 "velocity_n_mps": 25.0, "velocity_u_mps": 30.0,
                 "height_m": fa.PROFILE.depth_for_mu(fa.F_HZ, 0.5) if plasma else fa.HEIGHT_M}
        truth = float(run_scene(scene)["truth_mps"])
        return Setup(scene.scene_id, scene.group, scene.echo_class_truth, scene.medium, fa.F_HZ, facet,
                     {} if plasma else {"medium": "free"}, _starts(0.0, half), truth, math.nan, math.nan)

    known = {e.name: e for e in cv.ENVIRONMENTS}
    env = known.get(scene.environment) or cv.Environment(
        scene.environment, scene.f_mhz * 1e6, scene.medium == "plasma", mu_facet=float(scene.mu_facet))
    tilt = math.radians(scene.tilt_deg)
    depth = env.plate_depth()
    probe = cv._echo(env, Plate(x_m=0.0, depth_m=depth, length_m=1.0), np.array([0.0]))
    p_vertical, mu_s = float(probe.group_path_m[0]), probe.mu_sat
    extended = scene.group == "C_extended"
    w_own = float(envelope_width_s(env.f_hz * mu_s, p_vertical, min(scene.plate_m, 2.0e3) if extended
                                   else scene.plate_m))
    w_env = max(w_own, float(envelope_width_s(env.f_hz * mu_s, p_vertical, 2.0e3)))
    half = 2.5 * w_env if scene.half_window_s is None else float(scene.half_window_s)
    t_centre = (scene.x_m + cv.normal_crossing_offset_m(env, tilt)) / V_S
    facet = {"id": "f", "depth_m": depth, "x_m": scene.x_m, "tilt_deg": scene.tilt_deg,
             "length_m": scene.plate_m, "drift_mps": scene.drift_mps,
             "taper": "tukey" if extended else "uniform"}
    bias = float(cv.predicted_tilt_bias_mps(env, tilt)) if scene.tilt_deg else 0.0
    return Setup(scene.scene_id, scene.group, scene.echo_class_truth, scene.medium, env.f_hz, facet,
                 {} if env.plasma else {"medium": "free"}, _starts(t_centre, half), float(scene.drift_mps),
                 bias, float(envelope_width_s(env.f_hz * mu_s, p_vertical, scene.plate_m)))


def _prior(medium: str, scale_height_m: float = 200e3, source: str = "truth") -> Any:
    from tss.l3_signal.drift import OnboardPrior

    f_n = 1.5e6 if medium == "plasma" else 1.0e3  # free space: no plasma to speak of
    return OnboardPrior.exponential(f_n_sat_hz=f_n, scale_height_m=scale_height_m,
                                    ref_alt_m=PLATFORM_ALT_M, source=source)


class _ExponentialL0:
    """The orbit scenes' medium: ``c3_verify``'s exponential topside as a state."""

    name = "exponential"
    layer = "l0"

    def state_at(self, lat: float, lon: float, epoch: str, alt_m: np.ndarray) -> Any:
        from tss.interfaces import IonosphereState

        alt = np.asarray(alt_m, dtype=float)
        ne = 1.5e6**2 * np.exp((PLATFORM_ALT_M - alt) / 200e3) / 80.6163
        return IonosphereState(alt_m=alt, ne=ne)

    def describe(self) -> dict[str, Any]:
        return {"name": self.name, "layer": self.layer, "f_n_sat_hz": 1.5e6, "scale_height_m": 200e3}


def _instrument(oscillator: float | None, *, thermal: bool = False, tau_c_s: float | None = None) -> Any:
    """The bare receiver (``oscillator`` None) or the campaigns' instrument with it."""
    from tss.config import load_campaign
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.l2_instrument.fading import FadingModel
    from tss.studies import instrument_config

    if oscillator is None:
        return BaselineInstrument(quantise=False, thermal_noise=False)
    cfg = instrument_config(load_campaign("baseline").config, CONFIGURATION)
    cfg = {**cfg, "oscillator": {**dict(cfg.get("oscillator") or {}), "allan_deviation_1s": oscillator,
                                 "allan_deviation_2s": oscillator}}
    fading = None
    if tau_c_s is not None:
        fading = FadingModel(k_factor=math.inf, specular="rayleigh", tau_i_s=float(tau_c_s),
                             tau_log_sigma=0.0)
    return BaselineInstrument.from_config(cfg, thermal_noise=thermal, fading=fading)


def _gate_values(records: list[Any], fb: float) -> tuple[np.ndarray, int]:
    """The per-sweep value at beat ``fb`` after the first delay, and the samples summed."""
    x = np.stack([r.samples[0] for r in records])
    rec = records[0]
    tt = np.arange(x.shape[1]) / rec.fs_hz
    i0 = int(math.ceil(fb / rec.sweep.k_f_hz_per_s * rec.fs_hz)) + 3
    z = x[:, i0:] @ np.exp(-2j * math.pi * fb * (tt[i0:] - tt[i0:].mean()))
    return z, x.shape[1] - i0


def _simulate(setup: Setup, tiers: Mapping[str, Any], platforms: list[Any], truth_at=None) -> Any:
    from tss import sweep_law
    from tss.pipeline import Pipeline

    sweep = sweep_law.drift_mode(setup.f_hz)
    pipe = Pipeline(dict(tiers), seed=1)
    alt = np.linspace(50e3, 700e3, 1301)
    kw = {} if truth_at is None else {"truth_at": truth_at}
    return pipe.run_sequence(setup.scene_id, sweep=sweep, t_s=setup.t_s, platforms=platforms, alt_m=alt,
                             epoch=EPOCH, **kw), sweep


def _straight(setup: Setup) -> dict[str, Any]:
    """L1 and the bare receiver for a straight-track scene."""
    from tss.interfaces import PlatformState
    from tss.registry import build

    platforms = [PlatformState(alt_m=PLATFORM_ALT_M, lat_deg=0.0, lon_deg=0.0) for _ in setup.t_s]
    tiers = {"l0": build("l0", "chapman"), "l1": build("l1", "facet", facets=[setup.facet], **setup.l1_kwargs),
             "l2": _instrument(None)}
    run, sweep = _simulate(setup, tiers, platforms)
    return {"setup": setup, "echoes": list(run.echoes), "base": run.sequence, "sweep": sweep}


def _orbit(env_name: str) -> dict[str, Any]:
    """L1 and the bare receiver for one orbit scene (a 3-D plate under the track)."""
    from tss.l1_propagation.plate_echo import Plate
    from tss.platform.orbits import CircularOrbit, platform_state_at
    from tss.registry import build
    from tss.studies import c3_verify as cv
    from tss.studies.route_one import envelope_width_s

    env = {e.name: e for e in cv.ENVIRONMENTS}[env_name]
    depth = env.plate_depth()
    probe = cv._echo(env, Plate(x_m=0.0, depth_m=depth, length_m=1.0), np.array([0.0]))
    width = float(envelope_width_s(env.f_hz * probe.mu_sat, float(probe.group_path_m[0]), 2.0e3))
    orbit = CircularOrbit(PLATFORM_ALT_M, ORBIT_INCLINATION_DEG, u0_deg=ORBIT_U0_DEG)
    under = platform_state_at(orbit, 0.0)
    facet = {"id": "o", "lat_deg": under.lat_deg, "lon_deg": under.lon_deg,
             "alt_m": PLATFORM_ALT_M - depth, "drift_mps": ORBIT_DRIFT_MPS, "length_m": 2.0e3}
    setup = Setup(f"O_{env_name}", "O_orbit", "compact_facet", "plasma", env.f_hz, facet,
                  {"geometry": "orbit", "ref_alt_m": PLATFORM_ALT_M}, _starts(0.0, 2.5 * width),
                  ORBIT_DRIFT_MPS, 0.0, width)
    platforms = [platform_state_at(orbit, float(t) + 0.5 * SWEEP_S) for t in setup.t_s]
    tiers = {"l0": _ExponentialL0(), "l1": build("l1", "facet", facets=[facet], **setup.l1_kwargs),
             "l2": _instrument(None)}
    run, sweep = _simulate(setup, tiers, platforms, truth_at=(under.lat_deg, under.lon_deg))
    return {"setup": setup, "echoes": list(run.echoes), "base": run.sequence, "sweep": sweep}


def _facet_echo(prod: Any, fb: float, k: float) -> int | None:
    """The index of the echo within two range cells of beat ``fb``, or None."""
    if not len(prod.echo_id):
        return None
    beats = 2.0 * k * np.asarray(prod.products["range_beat_m"], dtype=float) / C_LIGHT
    j = int(np.argmin(np.abs(beats - fb)))
    return j if abs(beats[j] - fb) < NEAR_HZ else None


def _reference(got: dict[str, Any], s: Settings) -> dict[str, Any]:
    """Add the facet's beat and the unfaded peak against which every level is set.

    The beat comes from the tier on the bare record; the peak from the full
    instrument without thermal noise, with the oscillator at the threshold.
    """
    from tss.l3_signal.drift import DriftProcessing

    setup, base = got["setup"], got["base"]
    prod = DriftProcessing(prior=_prior(setup.medium)).process_sequence(base)
    if not len(prod.echo_id):
        raise RuntimeError(f"{setup.scene_id}: no echo in the noise-free record")
    j = int(np.argmax(prod.flags["detection"]["excess_db"]))
    k = got["sweep"].k_f_hz_per_s
    fb = 2.0 * k * float(prod.products["range_beat_m"][j]) / C_LIGHT
    records = _instrument(s.oscillator_main).receive_sequence(got["echoes"], got["sweep"], setup.t_s,
                                                              np.random.default_rng(s.seed))
    z, n_gate = _gate_values(records, fb)
    p = {n: v[j] for n, v in prod.products.items()}
    return {**got, "fb": fb, "ref_peak": float(np.max(np.abs(z) ** 2)), "n_gate": n_gate,
            "p_group_m": float(p["range_m"]), "mu_sat": float(p["mu_sat"])}


def needed_caches(work: list[Job], s: Settings) -> list[str]:
    """The prepared scenes the jobs read: the L1 run and the SNR reference of each."""
    keys: list[str] = []
    for j in work:
        if j.block in ("monte_carlo", "oscillator", "halo"):
            keys.append(f"straight:{j.scene_id}")
        elif j.block == "orbit":
            keys.append(f"orbit:{j.scene_id[2:]}")
    return list(dict.fromkeys(keys))


def prepare(key: str, s: Settings) -> dict[str, Any]:
    """One prepared scene (:func:`needed_caches`)."""
    from tss.studies.wp5_campaign import SCENES

    kind, name = key.split(":", 1)
    if kind == "straight":
        setup = scene_setup({sc.scene_id: sc for sc in SCENES}[name])
        return _reference(_straight(setup), s)
    return _reference(_orbit(name), s)


# ---------------------------------------------------------------------------
# One job
# ---------------------------------------------------------------------------


def _seed(s: Settings, *parts: Any) -> int:
    return zlib.crc32("|".join([str(s.seed), *(str(p) for p in parts)]).encode())


def _rows(job: Job, s: Settings, got: dict[str, Any], sequence: Any, level: Any, oscillator: str,
          priors: list[tuple[float, Any]], param_err: float) -> list[dict[str, Any]]:
    """Every estimator (and prior) on one received sequence."""
    from tss.l3_signal.drift import DriftProcessing
    from tss.l3_signal.novel.route_one_joint import sample_noise_variance

    setup, fb = got["setup"], got.get("fb", math.nan)
    k = got["sweep"].k_f_hz_per_s
    out: list[dict[str, Any]] = []
    for h_km, prior in priors:
        for est in s.estimators:
            prod = DriftProcessing(prior=prior, estimator=est).process_sequence(sequence)
            beat = fb
            if not math.isfinite(beat) and len(prod.echo_id):
                # A scene without a reference: its strongest echo is the target.
                j = int(np.argmax(prod.flags["detection"]["excess_db"]))
                beat = 2.0 * k * float(prod.products["range_beat_m"][j]) / C_LIGHT
            mine = _facet_echo(prod, beat, k)
            snr = math.nan
            det = prod.flags["detection"]
            if level != "none" and det.get("noise_power_per_bin") and "ref_peak" in got:
                var = sample_noise_variance(det["noise_power_per_bin"], det["segment_samples"])
                snr = 10.0 * math.log10(got["ref_peak"] / (got["n_gate"] * var))
            base = {
                "block": job.block, "scene_id": setup.scene_id, "group": setup.group,
                "echo_class_truth": setup.echo_class_truth, "medium": setup.medium,
                "f_mhz": setup.f_hz / 1e6, "level_db": level,
                "tau_c_s": math.nan if job.tau_c_s is None else job.tau_c_s, "oscillator": oscillator,
                "halo": job.halo, "h_assumed_km": h_km, "seed": job.seed, "estimator": est,
                "truth_mps": setup.truth_mps, "bias_mps": setup.bias_mps, "param_err_mps": param_err,
                "snr_db": snr, "n_sweeps": int(setup.t_s.size), "w_s": setup.w_s,
                "p_group_m": got.get("p_group_m", math.nan), "mu_sat": got.get("mu_sat", math.nan),
            }
            blank = {n: math.nan for n in ("drift_mps", "err_mps", "drift_sigma_mps", "drift_pinned_mps",
                                           "slope_ratio", "peak_contrast", "far_power_over_noise",
                                           "edge_power_over_noise", "tau_c_fit_s", "range_m",
                                           "blanked_fraction")}
            for j, echo_id in enumerate(prod.echo_id):
                v = {n: prod.products[n][j] for n in prod.products}

                def num(name: str, v: dict = v) -> float:
                    return float(v[name]) if name in v else math.nan

                drift = num("drift_mps")
                row = {**base, **blank, "echo_id": echo_id, "facet_echo": j == mine, "detected": True,
                       "formed": math.isfinite(drift), "echo_class": str(v.get("echo_class", "")),
                       "drift_mps": drift,
                       "err_mps": drift - setup.truth_mps if math.isfinite(drift) else math.nan,
                       "flags": ";".join(prod.flags["echo_flags"][j])}
                for name in ("drift_sigma_mps", "drift_pinned_mps", "slope_ratio", "peak_contrast",
                             "far_power_over_noise", "edge_power_over_noise", "tau_c_fit_s", "range_m",
                             "blanked_fraction"):
                    row[name] = num(name)
                out.append(row)
            if mine is None:
                out.append({**base, **blank, "echo_id": "", "facet_echo": True, "detected": False,
                            "formed": False, "echo_class": "", "flags": ""})
    return [{name: r[name] for name in ROW_FIELDS} for r in out]


def _noisy(sequence: Any, got: dict[str, Any], records: list[Any], level: float, seed: int) -> Any:
    """``records`` with white noise setting the facet's SNR per sweep to ``level`` dB."""
    sig2 = got["ref_peak"] / (got["n_gate"] * 10.0 ** (level / 10.0))
    rng = np.random.default_rng(seed)
    out = [dataclasses.replace(r, samples=r.samples + math.sqrt(sig2 / 2.0) * (
        rng.standard_normal(r.samples.shape) + 1j * rng.standard_normal(r.samples.shape))) for r in records]
    return dataclasses.replace(sequence, records=out)


def _levels(job: Job, s: Settings) -> tuple[Any, ...]:
    return {"monte_carlo": s.levels, "oscillator": s.oscillator_levels, "halo": s.halo_levels,
            "orbit": s.orbit_levels}[job.block]


def run_job(job: Job, s: Settings, cache: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The rows of one job.  ``cache`` holds the prepared scenes (:func:`prepare`)."""
    from tss.studies.wp5_campaign import SCENES, run_scene

    if job.block == "scenes":
        scene = {sc.scene_id: sc for sc in SCENES}[job.scene_id]
        setup = scene_setup(scene)
        got = _straight(setup)
        param = run_scene(scene)
        # For a sheet drift_err_mps is the error with the rotation about B assumed zero.
        param_err = float(param["drift_err_mps"])
        prior = _prior(setup.medium)
        return _rows(job, s, got, got["base"], "none", "none", [(math.nan, prior)], param_err)

    cached = cache[needed_caches([job], s)[0]]
    got = cached
    if job.block == "halo":
        k_factor = s.halo_cases[job.halo]
        if k_factor is not None:
            # A new halo each seed; the level is still set by the facet alone.
            facet = {**cached["setup"].facet, "halo": {
                "k": k_factor, "ratio": s.halo_ratio, "n_points": s.halo_points,
                "seed": _seed(s, "halo", job.halo, job.seed)}}
            run = _straight(dataclasses.replace(cached["setup"], facet=facet))
            got = {**cached, "echoes": run["echoes"], "base": run["base"]}
    setup = got["setup"]
    oscillator = "goal" if job.block == "oscillator" else "threshold"
    sigma = s.oscillator_second if job.block == "oscillator" else s.oscillator_main
    if job.block == "orbit":
        priors = [(h, _prior("plasma", h * 1e3, "assumed")) for h in s.orbit_scale_heights_km]
    else:
        priors = [(math.nan, _prior(setup.medium))]
    levels = _levels(job, s)
    rows: list[dict[str, Any]] = []
    stated = [lv for lv in levels if lv != "link"]
    if stated:
        records = _instrument(sigma, tau_c_s=job.tau_c_s).receive_sequence(
            got["echoes"], got["sweep"], setup.t_s, np.random.default_rng(_seed(s, job.job_id, "stated")))
        for level in stated:
            seq = _noisy(got["base"], got, records, float(level), _seed(s, job.job_id, level))
            rows += _rows(job, s, got, seq, float(level), oscillator, priors, math.nan)
    if "link" in levels:
        records = _instrument(sigma, thermal=True, tau_c_s=job.tau_c_s).receive_sequence(
            got["echoes"], got["sweep"], setup.t_s, np.random.default_rng(_seed(s, job.job_id, "stated")))
        seq = dataclasses.replace(got["base"], records=records)
        rows += _rows(job, s, got, seq, "link", oscillator, priors, math.nan)
    return rows


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

_GROUP = ("block", "scene_id", "level_db", "tau_c_s", "oscillator", "halo", "h_assumed_km", "estimator")


def _bound(first: Mapping[str, Any], level: Any, snr_db: float) -> float:
    """The Cramer-Rao bound of the group's facet at its nominal (or measured) SNR."""
    from tss.l3_signal.novel import drift_crb

    snr = float(level) if isinstance(level, (int, float)) else snr_db
    need = [float(first[n]) for n in ("w_s", "p_group_m", "mu_sat", "f_mhz")]
    if not (math.isfinite(snr) and all(math.isfinite(v) for v in need)):
        return math.nan
    n = int(first["n_sweeps"])
    t = (np.arange(n) - n // 2) / REPETITION_HZ
    tau = float(first["tau_c_s"])
    got = drift_crb.route_one_crb(drift_crb.fisher_route_one(
        t, snr_db=snr, width_s=need[0], f_hz=need[3] * 1e6, mu_sat=need[2], p_prime_m=need[1],
        drift_mps=25.0, tau_c_s=tau if math.isfinite(tau) else None))
    return float(got.sigma["v"])


def _median(values: list[float]) -> float:
    v = [x for x in values if math.isfinite(x)]
    return float(np.median(v)) if v else math.nan


def summarise(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Availability and error of the facet's echo for every combination."""
    groups: dict[tuple, list[dict[str, Any]]] = {}
    for r in rows:
        key = tuple("none" if isinstance(r[n], float) and math.isnan(r[n]) else r[n] for n in _GROUP)
        groups.setdefault(key, []).append(r)
    bounds: dict[tuple, float] = {}
    out = []
    for key, sel in groups.items():
        mine = [r for r in sel if r["facet_echo"]]
        if not mine:
            continue
        first = mine[0]
        formed = [r for r in mine if r["formed"]]
        err = np.array([r["err_mps"] for r in formed], dtype=float)
        snr = _median([float(r["snr_db"]) for r in mine])
        bound = math.nan
        if first["block"] != "scenes" and first["halo"] in ("none",) and first["echo_class_truth"] == "compact_facet":
            level = first["level_db"]
            where = (first["scene_id"], level, str(first["tau_c_s"]), round(snr, 1) if level == "link" else 0.0)
            if where not in bounds:
                bounds[where] = _bound(first, level, snr)
            bound = bounds[where]
        out.append({
            **dict(zip(_GROUP, key, strict=True)), "group": first["group"],
            "n": len(mine), "detected": sum(bool(r["detected"]) for r in mine), "formed": len(formed),
            "rms_mps": float(np.sqrt(np.mean(err**2))) if err.size else math.nan,
            "robust_mps": float(1.4826 * np.median(np.abs(err - np.median(err)))) if err.size else math.nan,
            "bias_mps": float(err.mean()) if err.size else math.nan,
            "bound_mps": bound,
            "outliers": float(np.mean(np.abs(err) > 3.0 * bound)) if err.size and math.isfinite(bound)
            else math.nan,
            "snr_db_median": snr,
            "peak_contrast_median": _median([float(r["peak_contrast"]) for r in mine]),
            "far_power_median": _median([float(r["far_power_over_noise"]) for r in mine]),
            "edge_power_median": _median([float(r["edge_power_over_noise"]) for r in mine]),
            "tau_c_fit_median": _median([float(r["tau_c_fit_s"]) for r in mine]),
            "sigma_median_mps": _median([float(r["drift_sigma_mps"]) for r in formed]),
            "other_echoes_formed": sum(bool(r["formed"]) and not r["facet_echo"] for r in sel),
        })
    return out


def _fmt(x: Any, spec: str = ".2f") -> str:
    if isinstance(x, float):
        return "-" if math.isnan(x) else format(x, spec)
    return str(x)


def _table(summary: list[dict[str, Any]], columns: list[tuple[str, str, str]]) -> list[str]:
    lines = ["| " + " | ".join(c[0] for c in columns) + " |", "|" + "---|" * len(columns)]
    for row in summary:
        lines.append("| " + " | ".join(_fmt(row[c[1]], c[2]) for c in columns) + " |")
    return lines


def _report(rows: list[dict[str, Any]], summary: list[dict[str, Any]], s: Settings,
            failed: list[str]) -> tuple[str, dict[str, Any]]:
    from tss.studies.wp5_campaign import EXPECTED

    lines = ["# wp5_route_one: route one at signal level", "",
             f"{len(rows)} rows; seeds per case {s.seeds}.  Definitions in the module "
             "`tss.studies.wp5_route_one`; every row in `wp5_route_one_echoes.csv.gz`."]
    metrics: dict[str, Any] = {"rows": len(rows), "jobs_failed": len(failed), "seeds": s.seeds}
    scenes = [r for r in rows if r["block"] == "scenes" and r["facet_echo"]]
    if scenes:
        lines += ["", "## Noise-free scenes", "",
                  "| group | estimator | scenes | formed | within 1 m/s | max abs error m/s | classed extended |",
                  "|---|---|---|---|---|---|---|"]
        for group in dict.fromkeys(r["group"] for r in scenes):
            for est in s.estimators:
                g = [r for r in scenes if r["group"] == group and r["estimator"] == est]
                err = [abs(r["err_mps"]) for r in g if r["formed"]]
                within = sum(e < 1.0 for e in err)
                ext = sum(r["echo_class"] == "extended" for r in g)
                lines.append(f"| {group} | {est} | {len(g)} | {len(err)} | {within} | "
                             f"{_fmt(max(err) if err else math.nan)} | {ext} |")
                metrics[f"{group}_{est}_formed"] = len(err)
                metrics[f"{group}_{est}_within_1mps"] = within
            tilted = [r for r in scenes if r["group"] == group and r["formed"] and r["bias_mps"]
                      and math.isfinite(r["bias_mps"])]
            if group == "B_tilted" and tilted:
                for est in s.estimators:
                    q = [r["err_mps"] / r["bias_mps"] for r in tilted if r["estimator"] == est]
                    if q:
                        lines.append(f"| {group} | {est}: error / predicted bias {min(q):.3f} to {max(q):.3f} "
                                     "| | | | | |")
    columns = [("scene", "scene_id", ""), ("level", "level_db", ".0f"), ("tau_c s", "tau_c_s", "g"),
               ("estimator", "estimator", ""), ("n", "n", "d"), ("detected", "detected", "d"),
               ("formed", "formed", "d"), ("RMS m/s", "rms_mps", ".2f"), ("robust", "robust_mps", ".2f"),
               ("bound", "bound_mps", ".2f"), ("outliers", "outliers", ".1%"), ("SNR dB", "snr_db_median", ".1f")]
    mc = [r for r in summary if r["block"] == "monte_carlo"]
    if mc:
        budget = []
        for r in mc:
            if r["group"] == "E_budget" and r["level_db"] == 10.0 and r["tau_c_s"] in ("none", 0.5):
                f = 5.0 if "5MHz" in r["scene_id"] and "10MHz" not in r["scene_id"] else 10.0
                total = EXPECTED["budget_total_mps"][(f, 500.0)][0 if r["tau_c_s"] == "none" else 1]
                budget.append({**r, "budget_mps": total})
        if budget:
            lines += ["", "## The budget points at 10 dB beside TV Table c3r1budget", "",
                      "The budget total includes the ephemeris and range-scale terms, which are not simulated.",
                      ""]
            lines += _table(budget, [("scene", "scene_id", ""), ("tau_c s", "tau_c_s", "g"),
                                     ("estimator", "estimator", ""), ("formed", "formed", "d"),
                                     ("n", "n", "d"), ("RMS m/s", "rms_mps", ".2f"),
                                     ("TV budget total m/s", "budget_mps", ".1f")])
        lines += ["", "## Drift-mode SNR per sweep at the link budget (TV assumption 1)", "",
                  "| scene | SNR per sweep, dB (median) |", "|---|---|"]
        for sid in dict.fromkeys(r["scene_id"] for r in mc):
            v = [r["snr_db_median"] for r in mc if r["scene_id"] == sid and r["level_db"] == "link"
                 and r["tau_c_s"] == "none"]
            if v:
                lines.append(f"| {sid} | {_fmt(float(v[0]), '.1f')} |")
                metrics[f"link_snr_db_{sid}"] = v[0]
        lines += ["", "## Monte Carlo", ""] + _table(mc, columns)
    osc = [r for r in summary if r["block"] == "oscillator"]
    if osc:
        lines += ["", "## Oscillator at the INS-26 goal (the Monte Carlo above has the threshold)", ""]
        lines += _table(osc, columns)
    halo = [r for r in summary if r["block"] == "halo"]
    if halo:
        lines += ["", "## Halo: availability and error against strength", ""]
        lines += _table(halo, [("halo", "halo", ""), ("level", "level_db", ".0f"), ("estimator", "estimator", ""),
                               ("n", "n", "d"), ("formed", "formed", "d"), ("RMS m/s", "rms_mps", ".2f"),
                               ("robust", "robust_mps", ".2f"), ("contrast", "peak_contrast_median", ".1f"),
                               ("far power", "far_power_median", ".2f"),
                               ("edge power", "edge_power_median", ".2f")])
    orbit = [r for r in summary if r["block"] == "orbit"]
    if orbit:
        lines += ["", "## Orbit: the prior's scale height", ""]
        lines += _table(orbit, [("scene", "scene_id", ""), ("level", "level_db", ".0f"),
                                ("H km", "h_assumed_km", ".0f"), ("estimator", "estimator", ""), ("n", "n", "d"),
                                ("formed", "formed", "d"), ("RMS m/s", "rms_mps", ".2f"),
                                ("bias", "bias_mps", ".2f")])
    if failed:
        lines += ["", "## Failed jobs", ""] + [f"- {f}" for f in failed]
    return "\n".join(lines) + "\n", metrics


# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------


def _write_json(path: Path, value: Any) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value), encoding="utf-8")
    os.replace(tmp, path)


def _worker_prepare(args: tuple[str, Settings, Path]) -> str:
    key, s, path = args
    if not path.exists():
        tmp = path.with_suffix(".tmp")
        tmp.write_bytes(pickle.dumps(prepare(key, s)))
        os.replace(tmp, path)
    return key


def _cache_path(folder: Path, key: str) -> Path:
    return folder / (key.replace(":", "__") + ".pkl")


def _worker_job(args: tuple[Job, Settings, Path, Path]) -> tuple[str, float, str]:
    job, s, cache_dir, part = args
    t0 = time.time()
    try:
        cache = {key: pickle.loads(_cache_path(cache_dir, key).read_bytes()) for key in needed_caches([job], s)}
        _write_json(part, run_job(job, s, cache))
    except Exception as exc:  # noqa: BLE001 - one job's failure must not lose the others; it is reported
        return job.job_id, time.time() - t0, f"{type(exc).__name__}: {exc}"
    return job.job_id, time.time() - t0, ""


def _map(pool: ProcessPoolExecutor | None, fn: Any, args: list[Any]):
    if pool is None:
        for a in args:
            yield fn(a)
        return
    for future in as_completed([pool.submit(fn, a) for a in args]):
        yield future.result()


def _check_parts(out_dir: Path, commit: str, s: Settings) -> None:
    """Parts may be reused only from the same commit and settings."""
    meta = out_dir / "parts_meta.json"
    wanted = {"commit": commit, "settings": json.loads(json.dumps(
        {k: v for k, v in dataclasses.asdict(s).items() if k != "processes"}, default=str))}
    if meta.exists():
        found = json.loads(meta.read_text(encoding="utf-8"))
        if found != wanted:
            raise RuntimeError(f"{out_dir} holds parts from another commit or other settings; "
                               "choose a new --out or delete its parts, cache and parts_meta.json")
    else:
        _write_json(meta, wanted)


def _file_name(job: Job) -> str:
    return f"{zlib.crc32(job.job_id.encode()):08x}_{job.block}_{job.seed}.json"


@register_study(
    "wp5_route_one",
    description="Route one at signal level: C3 scenes with noise, fading, oscillator, leakage and halo",
    work_packages=["WP5"],
    claims=["C3"],
    status="implemented",
    requirements=["SCI-15", "INS-15", "INS-26"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    from tss.provenance import git_state

    s = settings(config)
    out_dir = Path(out_dir)
    parts, cache_dir = out_dir / "parts", out_dir / "cache"
    for folder in (parts, cache_dir):
        folder.mkdir(parents=True, exist_ok=True)
    _check_parts(out_dir, str(git_state().get("commit")), s)
    work = jobs(s)
    manifest.seeds.update({"study": s.seed, "seeds_per_case": s.seeds})
    todo = [j for j in work if not (parts / _file_name(j)).exists()]
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ.setdefault(name, "1")  # one thread per worker process
    pool = ProcessPoolExecutor(s.processes) if s.processes > 1 else None
    started = time.time()
    failed: list[str] = []
    try:
        keys = needed_caches(todo, s)
        for key in _map(pool, _worker_prepare, [(k, s, _cache_path(cache_dir, k)) for k in keys]):
            print(f"prepared {key} [{time.time() - started:.0f} s]", flush=True)
        # Interleave the blocks, so every block's cost is known early and the time left is honest.
        order = sorted(range(len(todo)), key=lambda i: (todo[i].seed, todo[i].block, todo[i].job_id))
        args = [(todo[i], s, cache_dir, parts / _file_name(todo[i])) for i in order]
        spent: dict[str, list[float]] = {}
        left = {b: sum(j.block == b for j in todo) for b in BLOCKS}
        clock = time.time()
        for n, (job_id, seconds, error) in enumerate(_map(pool, _worker_job, args), start=1):
            block = job_id.split("__", 1)[0]
            spent.setdefault(block, []).append(seconds)
            left[block] -= 1
            if error:
                failed.append(f"{job_id}: {error}")
                print(f"FAILED {job_id}: {error}", flush=True)
            mean_all = float(np.mean([x for v in spent.values() for x in v]))
            remaining = sum(left[b] * (float(np.mean(spent[b])) if b in spent else mean_all) for b in BLOCKS)
            eta = remaining / max(s.processes, 1)
            if n % 20 == 0 or n == len(args) or n <= 3 or error:
                _write_json(out_dir / "progress.json", {
                    "done": n, "of": len(args), "failed": len(failed), "elapsed_s": time.time() - clock,
                    "remaining_s_if_all_processes_busy": eta,
                    "remaining_s_by_rate": (time.time() - clock) / n * (len(args) - n),
                    "mean_job_s": {b: float(np.mean(v)) for b, v in spent.items()}})
    finally:
        if pool is not None:
            pool.shutdown(cancel_futures=True)

    rows: list[dict[str, Any]] = []
    for job in work:
        part = parts / _file_name(job)
        if part.exists():
            rows += json.loads(part.read_text(encoding="utf-8"))
    summary = summarise(rows)
    text, metrics = _report(rows, summary, s, failed)
    echo_path = out_dir / "wp5_route_one_echoes.csv.gz"
    with gzip.open(echo_path, "wt", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(ROW_FIELDS))
        writer.writeheader()
        writer.writerows(rows)
    summary_path = out_dir / "wp5_route_one_summary.csv"
    with open(summary_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary[0]) if summary else ["block"])
        writer.writeheader()
        writer.writerows(summary)
    json_path = out_dir / "wp5_route_one.json"
    json_path.write_text(json.dumps({"settings": dataclasses.asdict(s), "metrics": metrics, "failed": failed,
                                     "jobs": len(work), "wall_s_this_run": time.time() - started},
                                    indent=2, default=str), encoding="utf-8")
    return StudyResult(outputs=[echo_path, summary_path, json_path], summary=text, metrics=metrics,
                       exit_code=1 if failed else 0)
