"""Layer: L3
WP: WP5
Claims: C3
Status: partial
Requirements: SCI-15, INS-15

The drift tier (``l3/drift``): route one at signal level on a drift-mode
sequence (WP5 Phase 3; ``docs/drift-mode.md`` sections 2 to 5).

WHAT IT DOES
------------
:meth:`DriftProcessing.process_sequence` takes the :class:`ReceivedSequence` of
``Pipeline.run_sequence`` and returns a :class:`DriftProduct` with the
contract's names.  One receive channel (channel 0).  Steps:

1.  **Detect.**  Each record is transformed over the samples after the longest
    delay searched (``max_range_m``), the power spectra are summed over the
    sequence, and up to ``max_echoes`` peaks above ``detect_db`` (over the
    median) become echoes, nearest first.  The beat of each is refined by a
    parabola through the zero-padded peak.
2.  **Gate.**  For each echo the first group delay of every record holds the
    previous sweep's echo at another beat; L2 leaves it empty and this tier
    drops it (decision (f), 2026-11-25): the gate starts ``guard_samples``
    after the echo's delay.  The fraction dropped is reported
    (``blanked_fraction``) and is a loss of SNR per sweep.
3.  **One complex value per sweep**, the gated record demodulated at the echo's
    beat with time measured from the gate's centre ``t_ref``, stamped with the
    time the echo was reflected, ``t_n + t_ref - tau / 2``.  Its phase is the
    echo's carrier phase then, at the frequency transmitted for it,
    ``f_eff = f_start + K (t_ref - tau)`` (D8 condition 2): velocities are
    converted at ``f_eff``, not at the sub-band's centre.
4.  **Platform Doppler.**  ``f_D,plat = (2 f_eff / c) mu_sat v_down`` for a
    local-down wave normal is removed by a phase ramp (section 2, step 1).
5.  **Route one** (:func:`estimate_route_one`, reused): Doppler from phase
    steps, the power envelope's peak ``t_pk``, the zero crossing ``t_zc`` with
    the ephemeris slope, and ``drift = v_eff^2 (t_zc - t_pk) / (mu_sat P')``,
    with its refusals (extended, no peak).

THE EFFECTIVE SPEED ON AN ORBIT
-------------------------------
Near overflight the phase path's second derivative is ``v^2 / P'_flat`` with
``P'_flat = int (r_p / r)^2 dz / mu`` (the group-path integral of the
Earth-flattened medium about the platform's radius), not ``v^2 / P'``.  So
route one uses ``v_eff^2 = v^2 P' / P'_flat``, which is ``v^2 r_f / r_p`` in
free space (4.3 % under ``v^2`` for a facet 300 km below a 600 km orbit).
``P'`` is measured; the facet's true depth, which ``P'_flat`` needs, comes from
the N(h) prior by inverting ``P' = int_0^D dz / mu`` (decision (d)).  The prior
is an explicit input (:class:`OnboardPrior`); ``Pipeline.run_sequence`` has no
ionogram chain to supply it yet, so a test or study passes one and its source
is in ``flags["prior_source"]`` (``"truth"`` when the truth profile was used).
On the straight track (platform states without velocity) the geometry is flat
and ``v_eff`` is ``v_s_mps``.

WHAT IS NOT DONE YET
--------------------
The diffuse halo, ``tau_c_s``, ``doppler_spread_hz``, the field-aligned and
several-frequency products, and ``drift_sigma_mps`` (NaN until the
Cramer-Rao harness, Phase 3).  Every facet product carries ``TILT_UNFLAGGED``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from tss.constants import C_LIGHT
from tss.interfaces import DriftProduct, ReceivedSequence
from tss.l3_signal.novel.route_one_estimator import estimate_route_one
from tss.platform.orbits import RE
from tss.registry import register

__all__ = ["OnboardPrior", "DriftProcessing", "effective_speed_ratio"]

#: ``f_N^2 / N_e``, Hz^2 m^3.
_FN2_PER_NE = 80.6163


@dataclass(frozen=True)
class OnboardPrior:
    """The N(h) profile the drift tier works from, and where it came from.

    ``alt_m``, ``ne``: the profile (m, m^-3), interpolated in ``log N_e``.
    ``source``: ``"truth"``, ``"onboard"`` (the ionogram chain's retrieval)
    or ``"configured"``; reported in the product's flags.
    """

    alt_m: np.ndarray
    ne: np.ndarray
    source: str = "configured"
    extra: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def exponential(cls, *, f_n_sat_hz: float = 1.5e6, scale_height_m: float = 200e3,
                    ref_alt_m: float = 600e3, source: str = "configured") -> OnboardPrior:
        """The reference topside of TV Ch. 8 about ``ref_alt_m``, on a fine grid."""
        alt = np.linspace(ref_alt_m - 900e3, ref_alt_m + 200e3, 4401)
        fn2 = float(f_n_sat_hz) ** 2 * np.exp((ref_alt_m - alt) / float(scale_height_m))
        return cls(alt_m=alt, ne=fn2 / _FN2_PER_NE, source=source,
                   extra={"f_n_sat_hz": float(f_n_sat_hz), "scale_height_m": float(scale_height_m),
                          "ref_alt_m": float(ref_alt_m)})

    def mu_squared(self, f_hz: float, alt_m: np.ndarray) -> np.ndarray:
        alt = np.asarray(alt_m, dtype=float)
        order = np.argsort(self.alt_m)
        ln = np.interp(alt, np.asarray(self.alt_m)[order], np.log(np.clip(np.asarray(self.ne)[order], 1.0, None)))
        return 1.0 - _FN2_PER_NE * np.exp(ln) / float(f_hz) ** 2

    def mu_at(self, f_hz: float, alt_m: float) -> float:
        m2 = float(self.mu_squared(f_hz, np.array([alt_m]))[0])
        return math.sqrt(m2) if m2 > 0.0 else math.nan


def _group_integrals(prior: OnboardPrior, f_hz: float, alt_p: float, depth: float, r_p: float,
                     n: int = 400) -> tuple[float, float]:
    """``int_0^D dz / mu`` and ``int_0^D (r_p / r)^2 dz / mu`` below ``alt_p``."""
    x, w = np.polynomial.legendre.leggauss(n)
    z = 0.5 * depth * (1.0 + x)
    hw = 0.5 * depth * w
    m2 = prior.mu_squared(f_hz, alt_p - z)
    if np.any(m2 <= 0.0):
        return math.nan, math.nan
    inv = 1.0 / np.sqrt(m2)
    ratio = r_p / (r_p - z)
    return float(np.sum(hw * inv)), float(np.sum(hw * ratio * ratio * inv))


def effective_speed_ratio(prior: OnboardPrior, f_hz: float, alt_p: float, group_path_m: float,
                          r_p: float) -> tuple[float, float]:
    """``(v_eff / v, depth)`` for an echo of group path ``P'`` (module docstring).

    The depth is found from ``P' = int_0^D dz / mu`` by bisection.
    """
    target = float(group_path_m)
    lo, hi = 0.0, target  # group path >= geometric depth
    if not math.isfinite(_group_integrals(prior, f_hz, alt_p, hi, r_p)[0]):
        # The prior's reflection level lies above: search below it.
        while hi - lo > 1.0 and not math.isfinite(_group_integrals(prior, f_hz, alt_p, hi, r_p)[0]):
            hi = 0.5 * (lo + hi)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        g, _ = _group_integrals(prior, f_hz, alt_p, mid, r_p)
        if g < target:
            lo = mid
        else:
            hi = mid
    depth = 0.5 * (lo + hi)
    g, g_flat = _group_integrals(prior, f_hz, alt_p, depth, r_p)
    return math.sqrt(g / g_flat), depth


@register("l3", "drift")
class DriftProcessing:
    """Route one at signal level on a drift-mode sequence.  See the module docstring."""

    layer = "l3"
    name = "drift"

    def __init__(
        self,
        *,
        prior: OnboardPrior | None = None,
        max_range_m: float = 1050e3,
        min_range_m: float = 20e3,
        max_echoes: int = 2,
        detect_db: float = 13.0,
        dynamic_db: float = 25.0,
        guard_samples: int = 3,
        zero_pad: int = 8,
        v_s_mps: float = 7600.0,
        earth_radius_m: float = RE,
        slope_tolerance: float = 0.2,
        min_contrast: float = 4.0,
        use_effective_speed: bool = True,
    ) -> None:
        self.prior = prior
        self.max_range_m = float(max_range_m)
        self.min_range_m = float(min_range_m)
        self.max_echoes = int(max_echoes)
        self.detect_db = float(detect_db)
        self.dynamic_db = float(dynamic_db)
        self.guard_samples = int(guard_samples)
        self.zero_pad = int(zero_pad)
        self.v_s_mps = float(v_s_mps)
        self.earth_radius_m = float(earth_radius_m)
        self.slope_tolerance = float(slope_tolerance)
        self.min_contrast = float(min_contrast)
        self.use_effective_speed = bool(use_effective_speed)

    def describe(self) -> Mapping[str, Any]:
        return {"name": self.name, "layer": self.layer, "max_range_m": self.max_range_m,
                "min_range_m": self.min_range_m, "max_echoes": self.max_echoes, "detect_db": self.detect_db, "dynamic_db": self.dynamic_db,
                "guard_samples": self.guard_samples, "zero_pad": self.zero_pad, "v_s_mps": self.v_s_mps,
                "earth_radius_m": self.earth_radius_m, "slope_tolerance": self.slope_tolerance,
                "min_contrast": self.min_contrast, "use_effective_speed": self.use_effective_speed,
                "prior_source": None if self.prior is None else self.prior.source}

    # -- steps 1 to 3 ----------------------------------------------------------------

    def _detect(self, x: np.ndarray, fs: float, k: float) -> list[float]:
        """Beats of the echoes, Hz, nearest first (step 1)."""
        n = x.shape[1]
        i0 = min(n - 8, int(math.ceil(2.0 * self.max_range_m / C_LIGHT * fs)))
        seg = x[:, i0:]
        m = seg.shape[1]
        nfft = int(2 ** math.ceil(math.log2(m * self.zero_pad)))
        win = np.hanning(m)
        p = np.sum(np.abs(np.fft.fft(seg * win, n=nfft, axis=1)) ** 2, axis=0)
        fr = np.fft.fftfreq(nfft, 1.0 / fs)
        lo, hi = k * 2.0 * self.min_range_m / C_LIGHT, k * 2.0 * self.max_range_m / C_LIGHT
        ok = np.nonzero((fr > lo) & (fr < hi))[0]
        if ok.size < 3:
            return []
        floor = float(np.median(p[ok]))
        thresh = floor * 10.0 ** (self.detect_db / 10.0)
        cell = max(3, int(round(4 * nfft / m)))  # four resolution cells either side
        found: list[float] = []
        p_top = float(p[ok].max())
        pp = p.copy()
        for _ in range(self.max_echoes):
            j = ok[int(np.argmax(pp[ok]))]
            if pp[j] < thresh or (found and pp[j] < p_top * 10.0 ** (-self.dynamic_db / 10.0)):
                break
            a, b, c = np.log(p[j - 1]), np.log(p[j]), np.log(p[j + 1])
            d = 0.5 * (a - c) / (a - 2.0 * b + c) if (a - 2.0 * b + c) != 0 else 0.0
            found.append(float(fr[j] + d * fs / nfft))
            pp[max(0, j - cell): j + cell + 1] = 0.0
        return sorted(found)

    # -- the tier ---------------------------------------------------------------------

    def process_sequence(self, sequence: ReceivedSequence, *, prior: OnboardPrior | None = None) -> DriftProduct:
        prior = prior if prior is not None else self.prior
        recs = list(sequence.records)
        if not recs:
            raise ValueError("an empty sequence")
        sweep = recs[0].sweep
        fs = float(recs[0].fs_hz)
        k = float(sweep.k_f_hz_per_s)
        f0 = float(sweep.f_start_hz)
        n_min = min(r.samples.shape[1] for r in recs)
        x = np.stack([np.asarray(r.samples[0][:n_min], dtype=complex) for r in recs])
        t_n = np.asarray(sequence.t_s, dtype=float)
        tt = np.arange(n_min) / fs
        beats = self._detect(x, fs, k)
        flags: dict[str, Any] = {
            "prior_source": None if prior is None else prior.source,
            "channel": 0, "n_sweeps": len(recs), "detected": len(beats),
            "not_done": ["diffuse_halo", "tau_c_s", "doppler_spread_hz", "drift_sigma_mps",
                         "field_aligned", "several_frequency"],
        }
        ids, rows = [], []
        for i, fb in enumerate(beats):
            ids.append(f"e{i + 1}")
            rows.append(self._echo(x, tt, t_n, fs, k, f0, fb, sequence.platforms, prior))
        names = sorted({key for r in rows for key in r if key != "flags"})
        products: dict[str, np.ndarray] = {}
        for key in names:
            vals = [r.get(key, math.nan) for r in rows]
            products[key] = np.array(vals, dtype=object if any(isinstance(v, str) for v in vals) else float)
        flags["echo_flags"] = [r["flags"] for r in rows]
        return DriftProduct(t_s=t_n, echo_id=ids, products=products, flags=flags,
                            source=f"l3/drift ({len(ids)} echoes)")

    def _echo(self, x, tt, t_n, fs, k, f0, fb, platforms: Sequence[Any], prior):
        n = tt.size
        tau0 = fb / k
        i0 = min(n - 8, int(math.ceil(tau0 * fs)) + self.guard_samples)
        g = slice(i0, n)
        t_ref = float(np.mean(tt[g]))
        osc = np.exp(-2j * math.pi * fb * (tt[g] - t_ref))
        z = x[:, g] @ osc  # one complex value per sweep, phase +2 pi f tau
        f_eff = f0 + k * (t_ref - tau0)
        # The geometry time: received at t_ref, reflected half a delay earlier
        # (L2 module docstring, geometry time and receive time).
        times = t_n + t_ref - 0.5 * tau0
        echo_flags: list[str] = ["TILT_UNFLAGGED"]
        p0 = platforms[0] if len(platforms) else None
        alt_p = float(getattr(p0, "alt_m", math.nan))
        mu_sat = prior.mu_at(f_eff, alt_p) if prior is not None and math.isfinite(alt_p) else math.nan
        if not math.isfinite(mu_sat):
            echo_flags.append("NO_PRIOR" if prior is None else "EVANESCENT")
            return {"range_m": C_LIGHT * tau0 / 2.0, "f_eff_hz": f_eff, "blanked_fraction": i0 / n,
                    "echo_class": "unresolved", "flags": echo_flags}
        # Step 4: the platform term for a local-down wave normal.
        v_down = np.array([float(np.asarray(p.velocity_ned_mps)[2]) if getattr(p, "velocity_ned_mps", None)
                           is not None else 0.0 for p in platforms])
        fd_plat = 2.0 * f_eff * mu_sat * v_down / C_LIGHT
        ramp = 2.0 * math.pi * np.concatenate([[0.0], np.cumsum(0.5 * (fd_plat[1:] + fd_plat[:-1]) * np.diff(times))])
        field_e = np.conj(z) * np.exp(-1j * ramp)  # plate convention exp(-j 2 k0 P): closing -> +Doppler
        power = np.abs(field_e) ** 2
        ipk = int(np.argmax(power))
        # Apparent range at the peak, debiased for the Doppler (beat = K tau - f_D).
        step = field_e[1:] * np.conj(field_e[:-1])
        fd_local = float(np.angle(step[min(ipk, step.size - 1)]) / (2.0 * math.pi * np.diff(times).mean()))
        p_group = C_LIGHT * (fb + fd_local) / k / 2.0
        geometry = "flat"
        v = self.v_s_mps
        ratio = 1.0
        depth = math.nan
        p_pk = platforms[ipk] if len(platforms) > ipk else p0
        vel = getattr(p_pk, "velocity_ned_mps", None)
        if vel is not None:
            vv = np.asarray(vel, dtype=float)
            v = float(math.hypot(vv[0], vv[1]))
            geometry = "orbit"
            if self.use_effective_speed:
                r_p = self.earth_radius_m + float(p_pk.alt_m)
                ratio, depth = effective_speed_ratio(prior, f_eff, float(p_pk.alt_m), p_group, r_p)
            else:
                echo_flags.append("FLAT_SPEED")
        v_eff = v * ratio
        est = estimate_route_one(times, field_e, f_hz=f_eff, group_path_m=p_group, mu_sat=mu_sat,
                                 v_s_mps=v_eff, slope_tolerance=self.slope_tolerance,
                                 min_contrast=self.min_contrast)
        if not est.separable:
            echo_flags.append("DECLINED: " + est.reason)
        echo_class = "compact_facet" if est.separable else (
            "extended" if abs(est.slope_ratio - 1.0) > self.slope_tolerance else "unresolved")
        return {
            "range_m": p_group, "range_beat_m": C_LIGHT * tau0 / 2.0, "f_eff_hz": f_eff, "mu_sat": mu_sat,
            "blanked_fraction": i0 / n, "t_ref_s": t_ref,
            "platform_doppler_hz": float(np.mean(fd_plat)), "los_rate_mps": est.los_rate_mps,
            "drift_mps": est.drift_mps, "drift_sigma_mps": math.nan, "t_pk_s": est.t_power_peak_s,
            "t_zc_s": est.t_zero_crossing_s, "slope_ratio": est.slope_ratio, "peak_contrast": est.peak_contrast,
            "along_track_m": v_eff * est.t_power_peak_s if est.separable else math.nan,
            "v_mps": v, "v_eff_mps": v_eff, "facet_depth_m": depth, "geometry": geometry,
            "power_peak_db": 10.0 * math.log10(float(power.max())) if power.max() > 0 else -math.inf,
            "echo_class": echo_class, "flags": echo_flags,
        }
