"""Route one of C3: drift of a compact scatterer from one transmit and receive pair.

Layer:        study
WP:           WP5
Claims:       C3
Status:       implemented
Requirements: SCI-15, INS-09, INS-15, INS-26

This study holds the numbers of Technical Volume Chapter 8 that describe route
one, so that the chapter can be regenerated rather than retyped.  It is the
analytic model and a Monte Carlo of one step of it.  It is not the estimator of
WP5, which has to work on simulated echoes and is unstarted.

The model
---------
A compact scatterer is fixed while the spacecraft passes over it at ``v_s``.
Inside the plasma its Doppler history near overflight is

    f_D(t) = 2 f mu_sat v_drift / c  -  (2 f / c) (v_s^2 / P') (t - t0)

where ``mu_sat`` is the refractive index at the spacecraft and ``P'`` is the
group path to the scatterer, which is the apparent range the sounder measures.
Both follow from the phase path (Budden 1961, eq. 10.14) and from Snell's law;
the derivation is in :mod:`tss.l1_propagation.doppler` and in Section 8.2 of
the Technical Volume.  The history crosses zero at

    t_zc = t0 + mu_sat v_drift P' / v_s^2

so the Doppler history alone fixes only ``mu_sat v_drift + (v_s^2 / P') t0``.
Drift and the overflight time ``t0`` are not separable from it.

What separates them is the echo power.  Monostatic specular reflection from a
facet of along-track size ``D`` is strongest when the spacecraft lies on the
facet normal, and drift moves the facet without rotating it, so the power peaks
at the geometric ``t0`` whatever the drift.  The field varies as
``sinc(2 pi D s / lambda)`` with ``s = v_s (t - t0) / P'``, and the power has a
full width at half maximum

    W = 0.443 lambda P' / (D v_s).

The Monte Carlo fits that envelope to noisy, fading power samples and reports
how well its centre is found.  The drift error from that term alone is
``(v_s^2 / (mu_sat P')) sigma_t0``.

The reference case
------------------
Chapter 8 quotes one plasma case: a plasma frequency at the spacecraft of
1.5 MHz, with ranges given as the apparent range ``P'``.  That puts 3 MHz at
the floor Requirements Rev G sizes the Doppler mode at, twice the plasma
frequency at the spacecraft, where ``mu_sat`` is 0.87.  The magnetic field is
neglected, as it is everywhere the analytic tier is used.

What is assumed and not shown
-----------------------------
The facet is an ideal uniformly illuminated strip, so its pattern is a sinc.
The signal-to-noise ratio is 10 dB per sweep at the peak of the envelope, which
INS-10 states for the ionogram mode and which the drift mode has still to be
shown to reach.  Fading is a complex Gaussian process with an exponential
correlation, which is a choice and not a measurement; the decorrelation time of
the topside is not known.  The fit is told the form of the envelope.  Each of
these is kind to the result.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tss.constants import C_LIGHT, doppler_range_bias
from tss.l1_propagation.doppler import (
    doppler_shift_hz,
    refractive_index_at_platform,
)
from tss.studies import StudyResult, register_study

__all__ = [
    "REFERENCE_F_N_SAT_HZ",
    "PLATFORM_SPEED_MPS",
    "INS26_THRESHOLD",
    "INS26_GOAL",
    "mu_sat",
    "doppler_slope_hz_per_s",
    "doppler_history_hz",
    "zero_crossing_shift_s",
    "drift_from_zero_crossing_mps",
    "envelope_width_s",
    "envelope_power",
    "overflight_time_tolerance_s",
    "known_t0_drift_sigma_mps",
    "OverflightTrial",
    "simulate_overflight_power",
    "fit_overflight_time",
    "overflight_time_sigma_s",
    "doppler_precision_mps",
    "oscillator_term_mps",
    "ephemeris_term_mps",
    "BudgetRow",
    "route_one_budget",
    "range_bias_table_m",
    "tilt_doppler_table",
]

#: Plasma frequency at the spacecraft in the reference case of TV Chapter 8.
REFERENCE_F_N_SAT_HZ = 1.5e6

#: Orbital speed the chapter uses.
PLATFORM_SPEED_MPS = 7600.0

#: INS-26, Requirements Rev G as amended 2026-10-06.  One tenth of the SCI-15
#: accuracy at mu_sat = sqrt(3)/2, kept unrounded so the budget reproduces 10.0.
INS26_THRESHOLD = 2.0 * 10.0 * np.sqrt(0.75) / C_LIGHT
INS26_GOAL = 3.2e-9

#: First zero of d/dx [sinc^2] half-power point: sinc^2(x) = 1/2 at this x.
_HALF_POWER_X = 1.3915573782515


def mu_sat(f_hz: float | np.ndarray, f_n_sat_hz: float = REFERENCE_F_N_SAT_HZ) -> np.ndarray:
    """Refractive index at the spacecraft in the reference case."""
    return refractive_index_at_platform(f_hz, f_n_sat_hz)


def doppler_slope_hz_per_s(
    f_hz: float | np.ndarray,
    p_group_m: float | np.ndarray,
    v_s_mps: float = PLATFORM_SPEED_MPS,
) -> np.ndarray:
    """Magnitude of the Doppler slope of a fixed compact scatterer, Hz/s.

    ``(2 f / c) v_s^2 / P'``.  It carries no refractive index: the group path
    takes the place of the geometric depth and that is the whole of the change.
    """
    f = np.asarray(f_hz, dtype=float)
    return 2.0 * f * float(v_s_mps) ** 2 / (C_LIGHT * np.asarray(p_group_m, dtype=float))


def doppler_history_hz(
    t_s: np.ndarray,
    f_hz: float,
    p_group_m: float,
    *,
    v_drift_mps: float = 0.0,
    t0_s: float = 0.0,
    mu: float = 1.0,
    v_s_mps: float = PLATFORM_SPEED_MPS,
) -> np.ndarray:
    """Linearised Doppler history of a compact scatterer near overflight, Hz."""
    t = np.asarray(t_s, dtype=float)
    slope = float(doppler_slope_hz_per_s(f_hz, p_group_m, v_s_mps))
    return float(doppler_shift_hz(v_drift_mps, f_hz, mu_sat=mu)) - slope * (t - float(t0_s))


def zero_crossing_shift_s(
    v_drift_mps: float | np.ndarray,
    p_group_m: float | np.ndarray,
    mu: float | np.ndarray = 1.0,
    v_s_mps: float = PLATFORM_SPEED_MPS,
) -> np.ndarray:
    """How far drift moves the zero crossing from the overflight time, s.

    ``mu_sat v_drift P' / v_s^2``.  A drift and an error of this size in ``t0``
    give the same Doppler history.
    """
    return (
        np.asarray(mu, dtype=float)
        * np.asarray(v_drift_mps, dtype=float)
        * np.asarray(p_group_m, dtype=float)
        / float(v_s_mps) ** 2
    )


def drift_from_zero_crossing_mps(
    t_zc_s: float | np.ndarray,
    t0_s: float | np.ndarray,
    p_group_m: float | np.ndarray,
    mu: float | np.ndarray = 1.0,
    v_s_mps: float = PLATFORM_SPEED_MPS,
) -> np.ndarray:
    """Invert :func:`zero_crossing_shift_s`: the route-one drift estimate."""
    dt = np.asarray(t_zc_s, dtype=float) - np.asarray(t0_s, dtype=float)
    return float(v_s_mps) ** 2 * dt / (
        np.asarray(mu, dtype=float) * np.asarray(p_group_m, dtype=float)
    )


def envelope_width_s(
    f_hz: float | np.ndarray,
    p_group_m: float | np.ndarray,
    facet_m: float | np.ndarray,
    v_s_mps: float = PLATFORM_SPEED_MPS,
) -> np.ndarray:
    """Full width at half maximum of the echo power envelope, s."""
    lam = C_LIGHT / np.asarray(f_hz, dtype=float)
    return (
        2.0 * _HALF_POWER_X / (2.0 * np.pi)
        * lam
        * np.asarray(p_group_m, dtype=float)
        / (np.asarray(facet_m, dtype=float) * float(v_s_mps))
    )


def envelope_power(t_s: np.ndarray, t0_s: float, width_s: float) -> np.ndarray:
    """Unit-peak ``sinc^2`` envelope with half-power width ``width_s``."""
    x = 2.0 * _HALF_POWER_X * (np.asarray(t_s, dtype=float) - float(t0_s)) / float(width_s)
    return np.sinc(x / np.pi) ** 2


def overflight_time_tolerance_s(
    sigma_v_mps: float | np.ndarray,
    p_group_m: float | np.ndarray,
    mu: float | np.ndarray = 1.0,
    v_s_mps: float = PLATFORM_SPEED_MPS,
) -> np.ndarray:
    """Error in ``t0`` that costs ``sigma_v_mps`` of drift."""
    return zero_crossing_shift_s(sigma_v_mps, p_group_m, mu, v_s_mps)


def known_t0_drift_sigma_mps(
    f_hz: float, sigma_f_hz: float, n_samples: int, mu: float = 1.0
) -> float:
    """Drift error from ``n`` independent Doppler samples when ``t0`` is known."""
    return float(C_LIGHT * sigma_f_hz / (2.0 * f_hz * mu) / np.sqrt(n_samples))


# ---------------------------------------------------------------------------
# The overflight-time Monte Carlo
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class OverflightTrial:
    """One realisation of the echo power across an overflight."""

    t_s: np.ndarray
    power: np.ndarray
    t0_s: float
    width_s: float


def simulate_overflight_power(
    width_s: float,
    rng: np.random.Generator,
    *,
    snr_db: float = 10.0,
    sample_rate_hz: float = 45.0,
    fading_tau_s: float | None = None,
    half_window_widths: float = 2.5,
) -> OverflightTrial:
    """Echo power per sweep across one overflight.

    The window runs ``half_window_widths`` envelope widths either side of its
    centre and the true overflight time is offset from that centre by up to one
    sample, so the fit is not handed a centre that sits on the grid.  Each sweep
    is ``sqrt(rho) * envelope * g + n`` with unit-variance complex Gaussian
    noise ``n``.  Without fading ``g`` is 1.  With fading it is a unit-variance
    complex Gaussian process whose correlation is ``exp(-|dt| / tau)``, which
    makes the echo amplitude Rayleigh.
    """
    dt = 1.0 / float(sample_rate_hz)
    half = float(half_window_widths) * float(width_s)
    n = int(np.floor(half / dt))
    t = np.arange(-n, n + 1, dtype=float) * dt
    t0 = float(rng.uniform(-dt, dt))
    rho = 10.0 ** (float(snr_db) / 10.0)
    field = np.sqrt(rho * envelope_power(t, t0, width_s))

    if fading_tau_s is None:
        g = np.ones(t.size, dtype=complex)
    else:
        r = float(np.exp(-dt / float(fading_tau_s)))
        w = (rng.standard_normal(t.size) + 1j * rng.standard_normal(t.size)) / np.sqrt(2.0)
        g = np.empty(t.size, dtype=complex)
        g[0] = w[0]
        innov = np.sqrt(1.0 - r * r)
        for k in range(1, t.size):
            g[k] = r * g[k - 1] + innov * w[k]

    noise = (rng.standard_normal(t.size) + 1j * rng.standard_normal(t.size)) / np.sqrt(2.0)
    power = np.abs(field * g + noise) ** 2
    return OverflightTrial(t_s=t, power=power, t0_s=t0, width_s=float(width_s))


def fit_overflight_time(t_s: np.ndarray, power: np.ndarray) -> tuple[float, float]:
    """Least-squares centre and width of a ``sinc^2`` envelope on a noise floor.

    Four parameters are fitted: centre, width, amplitude and floor.  The start
    comes from the data alone, a smoothed peak for the centre and the second
    moment about it for the width, so the fit is not told the answer.
    """
    from scipy.optimize import least_squares

    t = np.asarray(t_s, dtype=float)
    p = np.asarray(power, dtype=float)
    span = float(t[-1] - t[0])
    dt = float(t[1] - t[0])

    floor0 = float(np.median(p))
    k = max(3, int(round(0.08 * t.size)) | 1)
    smooth = np.convolve(p, np.ones(k) / k, mode="same")
    c0 = float(t[int(np.argmax(smooth))])
    excess = np.clip(smooth - floor0, 0.0, None)
    if excess.sum() > 0.0:
        var = float(np.sum(excess * (t - c0) ** 2) / excess.sum())
        w0 = float(np.clip(2.0 * np.sqrt(var), 4.0 * dt, 0.6 * span))
    else:
        w0 = 0.2 * span
    a0 = max(float(smooth.max() - floor0), 1e-3)

    def resid(q: np.ndarray) -> np.ndarray:
        return q[2] * envelope_power(t, q[0], q[1]) + q[3] - p

    sol = least_squares(
        resid,
        x0=[c0, w0, a0, floor0],
        bounds=([t[0], 2.0 * dt, 0.0, 0.0], [t[-1], span, np.inf, np.inf]),
    )
    return float(sol.x[0]), float(sol.x[1])


def overflight_time_sigma_s(
    width_s: float,
    *,
    n_trials: int = 400,
    seed: int = 0,
    snr_db: float = 10.0,
    sample_rate_hz: float = 45.0,
    fading_tau_s: float | None = None,
) -> float:
    """Robust standard deviation of the fitted overflight time, s.

    1.4826 times the median absolute deviation of the centre error over the
    trials.  The robust form is used because a fade can put the fitted centre
    on a sidelobe, and a handful of such trials would otherwise set the number.
    """
    rng = np.random.default_rng(seed)
    err = np.empty(int(n_trials), dtype=float)
    for i in range(int(n_trials)):
        trial = simulate_overflight_power(
            width_s, rng, snr_db=snr_db, sample_rate_hz=sample_rate_hz, fading_tau_s=fading_tau_s
        )
        centre, _ = fit_overflight_time(trial.t_s, trial.power)
        err[i] = centre - trial.t0_s
    return float(1.4826 * np.median(np.abs(err - np.median(err))))


# ---------------------------------------------------------------------------
# The budget
# ---------------------------------------------------------------------------


def doppler_precision_mps(
    f_hz: float,
    mu: float,
    *,
    ft_product_hz_s: float = 1.5e6,
    snr_db: float = 10.0,
    repetition_hz: float = 45.0,
) -> float:
    """Drift error from locating the zero crossing, m/s.

    The frequency of a tone of duration ``T`` is known to
    ``sqrt(6) / (2 pi T sqrt(rho f_rep T))``, and ``T`` is taken from the
    product ``f T``.
    """
    t = float(ft_product_hz_s) / float(f_hz)
    rho_int = 10.0 ** (float(snr_db) / 10.0) * float(repetition_hz) * t
    sigma_f = np.sqrt(6.0) / (2.0 * np.pi * t * np.sqrt(rho_int))
    return float(C_LIGHT * sigma_f / (2.0 * float(f_hz) * float(mu)))


def oscillator_term_mps(mu: float, sigma_y: float = INS26_THRESHOLD) -> float:
    """Drift error from a fractional frequency instability, ``c sigma_y / 2 mu``."""
    return float(C_LIGHT * float(sigma_y) / (2.0 * float(mu)))


def ephemeris_term_mps(
    p_group_m: float,
    mu: float,
    *,
    sigma_velocity_mps: float = 0.1,
    sigma_position_m: float = 10.0,
    v_s_mps: float = PLATFORM_SPEED_MPS,
) -> float:
    """Drift error from ephemeris velocity and position knowledge, m/s."""
    timing = float(sigma_position_m) / float(v_s_mps)
    from_position = float(v_s_mps) ** 2 * timing / (float(mu) * float(p_group_m))
    return float(np.hypot(float(sigma_velocity_mps), from_position))


@dataclass(frozen=True)
class BudgetRow:
    """One row of the route-one drift error budget.  Velocities in m/s."""

    f_hz: float
    p_group_m: float
    mu_sat: float
    width_s: float
    sigma_t0_s: float
    overflight_mps: float
    doppler_mps: float
    oscillator_mps: float
    ephemeris_mps: float
    range_scale_mps: float
    total_mps: float


def route_one_budget(
    f_hz: float,
    p_group_m: float,
    *,
    facet_m: float = 2.0e3,
    fading_tau_s: float | None = None,
    f_n_sat_hz: float = REFERENCE_F_N_SAT_HZ,
    sigma_y: float = INS26_THRESHOLD,
    range_scale_fraction: float = 0.05,
    drift_for_range_scale_mps: float = 50.0,
    n_trials: int = 400,
    seed: int = 0,
    v_s_mps: float = PLATFORM_SPEED_MPS,
) -> BudgetRow:
    """Evaluate the six-term budget for one frequency and apparent range.

    The range-and-Doppler coupling of the sixth term enters only through the
    range scale and is not carried separately.
    """
    mu = float(mu_sat(f_hz, f_n_sat_hz))
    width = float(envelope_width_s(f_hz, p_group_m, facet_m, v_s_mps))
    sigma_t0 = overflight_time_sigma_s(
        width, n_trials=n_trials, seed=seed, fading_tau_s=fading_tau_s
    )
    overflight = float(v_s_mps) ** 2 * sigma_t0 / (mu * float(p_group_m))
    doppler = doppler_precision_mps(f_hz, mu)
    osc = oscillator_term_mps(mu, sigma_y)
    eph = ephemeris_term_mps(p_group_m, mu, v_s_mps=v_s_mps)
    scale = float(range_scale_fraction) * float(drift_for_range_scale_mps)
    total = float(np.sqrt(overflight**2 + doppler**2 + osc**2 + eph**2 + scale**2))
    return BudgetRow(
        f_hz=float(f_hz),
        p_group_m=float(p_group_m),
        mu_sat=mu,
        width_s=width,
        sigma_t0_s=sigma_t0,
        overflight_mps=overflight,
        doppler_mps=doppler,
        oscillator_mps=osc,
        ephemeris_mps=eph,
        range_scale_mps=scale,
        total_mps=total,
    )


# ---------------------------------------------------------------------------
# The chapter's small tables
# ---------------------------------------------------------------------------


def a1_f_over_k_s(f_hz: float | np.ndarray, *, dwell_s: float = 0.0268) -> np.ndarray:
    """``f / K_f`` for the A1 stepped sweep, s.

    Below the 3 MHz crossover a step spans ``0.05 f`` in one dwell, so the ratio
    is ``dwell / 0.05``.  Above it the step is 150 kHz and the rate is constant.
    """
    f = np.asarray(f_hz, dtype=float)
    k_peak = 1.5e5 / float(dwell_s)
    return np.where(f < 3.0e6, float(dwell_s) / 0.05, f / k_peak)


def range_bias_table_m(
    velocities_mps: np.ndarray,
    f_hz: np.ndarray,
    f_n_sat_hz: float = REFERENCE_F_N_SAT_HZ,
) -> np.ndarray:
    """Apparent range bias ``mu_sat v f / K_f`` for the A1 sweep, m.

    Rows are velocities and columns are frequencies.  NaN where the frequency
    is at or below the plasma frequency at the spacecraft.
    """
    v = np.asarray(velocities_mps, dtype=float)[:, None]
    f = np.asarray(f_hz, dtype=float)[None, :]
    free = doppler_range_bias(v, 1.0, 1.0 / a1_f_over_k_s(f))
    return mu_sat(f, f_n_sat_hz) * free


def tilt_doppler_table(
    tilts_deg: np.ndarray,
    f_hz: np.ndarray,
    f_n_sat_hz: float = REFERENCE_F_N_SAT_HZ,
    v_s_mps: float = PLATFORM_SPEED_MPS,
) -> tuple[np.ndarray, np.ndarray]:
    """Apparent radial velocity and Doppler shift from a tilted layer.

    Returns ``(v_r, f_D)`` with ``v_r`` one value per tilt and ``f_D`` an array
    of tilts by frequencies.
    """
    v_r = float(v_s_mps) * np.sin(np.radians(np.asarray(tilts_deg, dtype=float)))
    f = np.asarray(f_hz, dtype=float)
    f_d = np.array([doppler_shift_hz(v, f, mu_sat=mu_sat(f, f_n_sat_hz)) for v in v_r])
    return v_r, f_d


# ---------------------------------------------------------------------------
# The registered study
# ---------------------------------------------------------------------------

#: Facet sizes of the Monte Carlo figure, km.  3.9 km is the Fresnel radius at
#: 5 MHz and 500 km in free space, above which a facet is not compact.
FACETS_KM = (0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 3.9)

#: Fading cases: label and correlation time, s.
FADING_CASES = (("none", None), ("tau_0p1", 0.1), ("tau_0p5", 0.5))

#: Rows of the budget table: frequency in MHz and apparent range in km.
BUDGET_CASES = ((3.0, 500.0), (3.0, 1000.0), (5.0, 500.0), (5.0, 1000.0), (10.0, 500.0), (10.0, 1000.0))


@register_study(
    "route_one",
    description="Route one of C3: the overflight-time Monte Carlo and the drift error budget",
    work_packages=["WP5"],
    claims=["C3"],
    status="implemented",
    requirements=["SCI-15", "INS-09", "INS-15", "INS-26"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    opts = dict(config.get("route_one", {}))
    n_trials = int(opts.get("n_trials", 1000))
    seed = int(opts.get("seed", 20261006))
    f_n_sat = float(opts.get("f_n_sat_hz", REFERENCE_F_N_SAT_HZ))

    f_ref, p_ref = 5.0e6, 5.0e5
    mu_ref = float(mu_sat(f_ref, f_n_sat))

    mc_rows = []
    for d_km in FACETS_KM:
        width = float(envelope_width_s(f_ref, p_ref, d_km * 1e3))
        for j, (_label, tau) in enumerate(FADING_CASES):
            sig = overflight_time_sigma_s(
                width, n_trials=n_trials, seed=seed + 100 * j + int(round(100 * d_km)), fading_tau_s=tau
            )
            drift = PLATFORM_SPEED_MPS**2 * sig / (mu_ref * p_ref)
            mc_rows.append((d_km, j, width, sig, drift, sig / width))
    mc_path = out_dir / "route_one_overflight_mc.csv"
    np.savetxt(
        mc_path, np.array(mc_rows), delimiter=",", comments="",
        header="facet_km,fading_case,width_s,sigma_t0_s,drift_error_mps,kappa",
    )

    budget_rows = []
    for i, (f_mhz, p_km) in enumerate(BUDGET_CASES):
        for j, tau in enumerate((None, 0.5)):
            row = route_one_budget(
                f_mhz * 1e6, p_km * 1e3, fading_tau_s=tau, f_n_sat_hz=f_n_sat,
                n_trials=n_trials, seed=seed + 1000 + 10 * i + j,
            )
            budget_rows.append(
                (f_mhz, p_km, j, row.mu_sat, row.width_s, row.sigma_t0_s, row.overflight_mps,
                 row.doppler_mps, row.oscillator_mps, row.ephemeris_mps, row.range_scale_mps,
                 row.total_mps)
            )
    budget_path = out_dir / "route_one_budget.csv"
    np.savetxt(
        budget_path, np.array(budget_rows), delimiter=",", comments="",
        header=("f_mhz,p_group_km,fading,mu_sat,width_s,sigma_t0_s,overflight_mps,doppler_mps,"
                "oscillator_mps,ephemeris_mps,range_scale_mps,total_mps"),
    )

    b = np.array(budget_rows)
    clear, faded = b[b[:, 2] == 0], b[b[:, 2] == 1]
    mc = np.array(mc_rows)
    metrics = {
        "f_n_sat_hz": f_n_sat,
        "n_trials": n_trials,
        "seed": seed,
        "mu_sat_3mhz": float(mu_sat(3.0e6, f_n_sat)),
        "mu_sat_5mhz": mu_ref,
        "mu_sat_10mhz": float(mu_sat(1.0e7, f_n_sat)),
        "zero_crossing_shift_s_25mps_5mhz_500km": float(zero_crossing_shift_s(25.0, p_ref, mu_ref)),
        "total_mps_no_fading_min": float(clear[:, 11].min()),
        "total_mps_no_fading_max": float(clear[:, 11].max()),
        "total_mps_fading_min": float(faded[:, 11].min()),
        "total_mps_fading_max": float(faded[:, 11].max()),
        "drift_error_mps_no_fading_min": float(mc[mc[:, 1] == 0][:, 4].min()),
        "drift_error_mps_no_fading_max": float(mc[mc[:, 1] == 0][:, 4].max()),
        "drift_error_mps_tau_0p5_max": float(mc[mc[:, 1] == 2][:, 4].max()),
    }

    summary = f"""
Route one measures the drift of a compact scatterer from the zero crossing of
its Doppler history and the peak of its echo power.  The reference case is a
plasma frequency at the spacecraft of {f_n_sat / 1e6:.1f} MHz, so the refractive index
there is {metrics["mu_sat_3mhz"]:.2f} at 3 MHz, {mu_ref:.2f} at 5 MHz and
{metrics["mu_sat_10mhz"]:.2f} at 10 MHz.  Ranges are apparent ranges.

For a 2 km facet at 10 dB per sweep and 45 Hz the total drift error is
{metrics["total_mps_no_fading_min"]:.1f} to {metrics["total_mps_no_fading_max"]:.1f} m/s without fading and
{metrics["total_mps_fading_min"]:.1f} to {metrics["total_mps_fading_max"]:.1f} m/s with fading at a correlation time of
0.5 s, over 3 to 10 MHz and 500 to 1000 km.  The oscillator allocation of
INS-26 dominates the first figure and the overflight time dominates the second.

The envelope is an ideal sinc, the signal-to-noise ratio is assumed, and the
fading model is a choice.  {n_trials} trials per case, seed {seed}.
""".strip()

    return StudyResult(outputs=[mc_path, budget_path], summary=summary, metrics=metrics)
