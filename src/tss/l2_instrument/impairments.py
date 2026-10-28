"""Layer: L2
WP: WP3
Claims: C4, C3
Status: implemented
Requirements: INS-21, INS-22, INS-24, INS-25, INS-26

Transmitter impairments, applied to the samples (WP3 Phase 2).

WHAT THIS MODULE IS FOR
-----------------------
Every impairment of the transmitted chirp reaches the receiver twice: on the
leakage, a few nanoseconds late, and on each echo, milliseconds late.  The
de-chirp reference subtracts whatever part of the impairment it shares with the
transmitted signal, delayed by the path.  Rev G section 5.9 derives what
survives for one sinusoidal frequency error; this module puts the real set into
the samples so the derivation can be tested and then gone beyond:

* **Reference-oscillator noise.**  The fractional frequency noise ``y(t)`` of
  the reference, with the power-law spectrum its Allan deviation and noise type
  imply (:class:`tss.l2_instrument.oscillator.Oscillator`), multiplied up to the
  instantaneous carrier: the phase error is ``2 pi int f_c(t) y(t) dt``.
* **Synthesiser phase-noise floor.**  Additive phase noise, flat at
  ``phase_noise_floor_dbc_hz`` (single-sideband), not scaled by the carrier.
* **AM noise.**  Flat at ``am_noise_dbc_hz`` (single-sideband) on the radiated
  amplitude.  The de-chirp mixer's reference port is driven into saturation and
  does not carry the transmit amplifier's AM, so AM is never range-correlated:
  it reaches the residual whole, whatever the delay (Committee Review 5).
* **Discrete spurs** (:class:`Spur`): phase-modulation spurs from the
  synthesiser, which the reference shares under INS-21, and amplitude spurs,
  which it does not.  INS-22 is written on these.
* **Sweep nonlinearity.**  ``SweepSchedule.nonlinearity_ppm`` on the transmitted
  law, per step (INS-25), when ``sweep_nonlinearity`` is on.  Until Phase 2 the
  L2 tier always used the ideal law.
* **Doppler.**  ``EchoSet.doppler_hz``, which L1 computes with the refractive
  index at the spacecraft (B4), applied to the echoes when ``doppler`` is on.
  An echo closing at ``f_d`` beats at ``f_R - f_d`` on an up-sweep and
  ``-(f_R + f_d)`` on a down-sweep (:mod:`tss.l1_propagation.doppler`).

``single_reference`` is INS-21.  When it is on, the de-chirp reference carries
the same oscillator noise, the same synthesiser floor, the same PM spurs and the
same sweep law as the transmitted chirp, and only their lag-``tau`` difference
survives.  When it is off, the reference is an independent source of the same
specification: its noise and spur phases are independent draws, and its sweep
law is ideal.  The difference between the two is the range-correlation credit,
which the tests hold to the prediction.

HOW THE NOISE IS MADE AND DELAYED
---------------------------------
Each noise process is synthesised once per record, in the frequency domain, on
a grid that extends ``max_delay_s`` before the record so that the longest echo
has a past.  A process synthesised this way is periodic on its grid, so a delay
by any ``tau``, integer number of samples or not, is exact as a phase ramp on
its spectrum.  The oscillator term is the running integral of a delayed
``f_c y``, which is exact to within a constant phase of order
``2 pi f_c y tau`` (below 1e-4 rad).  One :class:`ImpairmentRealization` serves
the echoes and the leakage of the same record, as one oscillator would.

Each noise process is band-limited to the digitiser's Nyquist rate.  Noise
above it is taken as removed by the anti-alias filter.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from dataclasses import replace as _replace
from typing import Any

import numpy as np
from scipy import fft as _fft
from scipy import special as _special

from tss import sweep_law
from tss.interfaces import SweepSchedule
from tss.l2_instrument.oscillator import FM_NOISE_TYPES, Oscillator

__all__ = [
    "Spur",
    "TransmitImpairments",
    "ImpairmentRealization",
    "power_law_h",
    "pm_spur_beta",
    "am_spur_index",
    "dechirped_pm_density_dbc_hz",
    "oscillator_phase_density_dbc_hz",
    "local_dechirp_phase_rad",
]


# ----------------------------------------------------------------------------
# The specification
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class Spur:
    """A discrete spur on the radiated chirp, ``level_dbc`` per sideband.

    ``kind`` is ``"pm"`` (phase modulation from the synthesiser: shared with
    the reference under INS-21) or ``"am"`` (amplitude: never shared).
    """

    offset_hz: float
    level_dbc: float
    kind: str = "pm"
    phase_rad: float = 0.0

    def __post_init__(self) -> None:
        if not self.offset_hz > 0.0:
            raise ValueError("a spur needs a positive offset")
        if self.kind not in ("pm", "am"):
            raise ValueError("spur kind must be 'pm' or 'am'")
        if not self.level_dbc < 0.0:
            raise ValueError("a spur's level must be below the carrier")


def pm_spur_beta(level_dbc: float) -> float:
    """Phase-modulation index whose first sideband is ``level_dbc``, exactly.

    Solves ``J1(beta) / J0(beta) = 10^(level / 20)``; for a small spur
    ``beta = 2 * 10^(level / 20)``.
    """
    from scipy.optimize import brentq

    target = 10.0 ** (float(level_dbc) / 20.0)
    return float(brentq(lambda b: _special.j1(b) / _special.j0(b) - target, 0.0, 2.0))


def am_spur_index(level_dbc: float) -> float:
    """Amplitude-modulation index whose sidebands are each ``level_dbc``: ``m = 2 * 10^(L/20)``."""
    return 2.0 * 10.0 ** (float(level_dbc) / 20.0)


def power_law_h(oscillator: Oscillator, tau_s: float = 1.0) -> tuple[float, int]:
    """``(h_a, a)`` of ``S_y(f) = h_a f^a`` (one-sided) from the Allan deviation at ``tau_s``.

    White FM: ``sigma^2 = h_0 / (2 tau)``.  Flicker FM: ``sigma^2 = 2 ln2 h_-1``.
    Random-walk FM: ``sigma^2 = (2 pi^2 / 3) h_-2 tau``.
    """
    a = int(FM_NOISE_TYPES[oscillator.noise_type])
    s2 = float(oscillator.sigma_y(float(tau_s))) ** 2
    tau = float(tau_s)
    if a == 0:
        return 2.0 * tau * s2, 0
    if a == -1:
        return s2 / (2.0 * math.log(2.0)), -1
    return 3.0 * s2 / (2.0 * math.pi**2 * tau), -2


@dataclass(frozen=True)
class TransmitImpairments:
    """What is wrong with the transmitted chirp, each part switchable.

    Parameters
    ----------
    oscillator
        The reference oscillator whose fractional-frequency noise is applied,
        scaled to the carrier, or ``None`` for none.
    phase_noise_floor_dbc_hz
        Additive synthesiser phase-noise floor, single-sideband, or ``None``.
    am_noise_dbc_hz
        Flat AM noise of the radiated signal, single-sideband, or ``None``.
    spurs
        :class:`Spur` terms.
    sweep_nonlinearity
        Apply ``SweepSchedule.nonlinearity_ppm`` to the transmitted law.
    single_reference
        INS-21.  On: the reference shares the oscillator noise, the
        synthesiser floor, the PM spurs and the sweep law.  Off: it is an
        independent source of the same specification with an ideal law.
    doppler
        Apply ``EchoSet.doppler_hz`` to the echoes (no effect on leakage).
    """

    oscillator: Oscillator | None = None
    phase_noise_floor_dbc_hz: float | None = None
    am_noise_dbc_hz: float | None = None
    spurs: tuple[Spur, ...] = ()
    sweep_nonlinearity: bool = False
    single_reference: bool = True
    doppler: bool = False
    label: str = field(default="", compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "spurs", tuple(self.spurs))

    @property
    def needs_realization(self) -> bool:
        """True when any random or modulating term is on."""
        return (
            self.oscillator is not None
            or self.phase_noise_floor_dbc_hz is not None
            or self.am_noise_dbc_hz is not None
            or bool(self.spurs)
        )

    @property
    def any(self) -> bool:
        return self.needs_realization or self.sweep_nonlinearity or self.doppler

    @classmethod
    def none(cls) -> TransmitImpairments:
        return cls()

    @classmethod
    def from_config(
        cls, cfg: Mapping[str, Any] | None, oscillator: Oscillator | None = None
    ) -> TransmitImpairments:
        """From a ``transmit_impairments`` block.

        Keys: ``oscillator_noise`` (bool; uses ``oscillator``),
        ``phase_noise_floor_dbc_hz``, ``am_noise_dbc_hz``, ``spurs`` (a list of
        mappings with ``offset_hz``, ``level_dbc``, ``kind``),
        ``sweep_nonlinearity``, ``single_reference``, ``doppler``.  An absent or
        empty block is no impairment, which is what every run before Phase 2
        had.
        """
        c = dict(cfg or {})
        osc = oscillator if bool(c.get("oscillator_noise", False)) else None
        if osc is not None and "single_reference" not in c:
            c["single_reference"] = osc.single_reference

        def _opt(key: str) -> float | None:
            v = c.get(key)
            return None if v is None else float(v)

        return cls(
            oscillator=osc,
            phase_noise_floor_dbc_hz=_opt("phase_noise_floor_dbc_hz"),
            am_noise_dbc_hz=_opt("am_noise_dbc_hz"),
            spurs=tuple(Spur(**dict(s)) for s in (c.get("spurs") or ())),
            sweep_nonlinearity=bool(c.get("sweep_nonlinearity", False)),
            single_reference=bool(c.get("single_reference", True)),
            doppler=bool(c.get("doppler", False)),
            label=str(c.get("label", "")),
        )

    def describe(self) -> Mapping[str, Any]:
        return {
            "oscillator": None if self.oscillator is None else dict(self.oscillator.describe()),
            "phase_noise_floor_dbc_hz": self.phase_noise_floor_dbc_hz,
            "am_noise_dbc_hz": self.am_noise_dbc_hz,
            "spurs": [asdict(s) for s in self.spurs],
            "sweep_nonlinearity": self.sweep_nonlinearity,
            "single_reference": self.single_reference,
            "doppler": self.doppler,
            "label": self.label,
        }

    def realize(
        self,
        sweep: SweepSchedule,
        fs_hz: float,
        n_samples: int,
        rng: np.random.Generator,
        *,
        t0_s: float = 0.0,
        max_delay_s: float = 0.0,
    ) -> ImpairmentRealization:
        """One draw of every random term, for a record of ``n_samples`` from ``t0_s``."""
        return ImpairmentRealization(
            self, sweep, fs_hz, n_samples, rng, t0_s=t0_s, max_delay_s=max_delay_s
        )


# ----------------------------------------------------------------------------
# One draw
# ----------------------------------------------------------------------------


def _carrier_hz(sweep: SweepSchedule, t: np.ndarray) -> np.ndarray:
    tt = np.clip(t, 0.0, None)
    if sweep_law.is_stepped(sweep):
        return sweep_law.frequency_hz(sweep, tt)
    d = 1.0 if int(sweep.direction) >= 0 else -1.0
    f0 = float(sweep.f_start_hz) if d > 0 else float(sweep.f_stop_hz)
    return f0 + d * float(sweep.k_f_hz_per_s) * tt


def _synth_spectrum(s_one_sided: np.ndarray, m: int, fs: float, rng: np.random.Generator) -> np.ndarray:
    """rfft coefficients of a real process of one-sided PSD ``s`` on ``m`` points."""
    df = fs / m
    k = s_one_sided.size
    z = rng.normal(size=k) + 1j * rng.normal(size=k)
    spec = m * np.sqrt(np.clip(s_one_sided, 0.0, None) * df / 4.0) * z
    spec[0] = 0.0
    if m % 2 == 0:
        spec[-1] = 0.0
    return spec


class ImpairmentRealization:
    """One draw of the impairments on one record, shared by echoes and leakage.

    The record is ``n`` samples at ``fs`` from ``t0``.  Methods return arrays on
    that record:

    * :meth:`tx_phase` (``tau``): the transmit phase error at ``t - tau``;
    * :meth:`ref_phase`: the de-chirp reference's phase error at ``t``;
    * :meth:`tx_amplitude` (``tau``): the radiated amplitude factor at
      ``t - tau``, ``1 + a``.
    """

    def __init__(
        self,
        imp: TransmitImpairments,
        sweep: SweepSchedule,
        fs_hz: float,
        n: int,
        rng: np.random.Generator,
        *,
        t0_s: float = 0.0,
        max_delay_s: float = 0.0,
    ) -> None:
        self.imp = imp
        self.sweep = sweep
        self.fs = float(fs_hz)
        self.n = int(n)
        self.t0 = float(t0_s)
        self.pad = int(math.ceil(max(float(max_delay_s), 0.0) * self.fs)) + 16
        self.m = int(_fft.next_fast_len(self.n + self.pad, real=True))
        self.t_ext = self.t0 + (np.arange(self.m) - self.pad) / self.fs
        self.t = self.t_ext[self.pad : self.pad + self.n]
        self.f = np.fft.rfftfreq(self.m, 1.0 / self.fs)
        fpos = self.f.copy()
        fpos[0] = np.inf
        self._fc = _carrier_hz(sweep, self.t_ext)

        def _osc_spec() -> np.ndarray | None:
            if imp.oscillator is None:
                return None
            h, a = power_law_h(imp.oscillator)
            s_y = h * fpos ** float(a)
            s_y[0] = 0.0
            y = _fft.irfft(_synth_spectrum(s_y, self.m, self.fs, rng), n=self.m)
            return _fft.rfft(self._fc * y)

        def _flat_spec(l_dbc: float | None) -> np.ndarray | None:
            if l_dbc is None:
                return None
            s = np.full(self.f.shape, 2.0 * 10.0 ** (float(l_dbc) / 10.0))
            return _synth_spectrum(s, self.m, self.fs, rng)

        # Transmit chain draws.
        self._g_tx = _osc_spec()
        self._pm_tx = _flat_spec(imp.phase_noise_floor_dbc_hz)
        self._am = _flat_spec(imp.am_noise_dbc_hz)
        self._spur_phase_tx = [float(s.phase_rad) for s in imp.spurs]
        # The reference: shared under INS-21, independent otherwise.
        if imp.single_reference:
            self._g_ref, self._pm_ref = self._g_tx, self._pm_tx
            self._spur_phase_ref = list(self._spur_phase_tx)
        else:
            self._g_ref = _osc_spec()
            self._pm_ref = _flat_spec(imp.phase_noise_floor_dbc_hz)
            self._spur_phase_ref = [float(rng.uniform(0.0, 2.0 * np.pi)) for _ in imp.spurs]
        self._betas = [pm_spur_beta(s.level_dbc) if s.kind == "pm" else 0.0 for s in imp.spurs]
        self._m_am = [am_spur_index(s.level_dbc) if s.kind == "am" else 0.0 for s in imp.spurs]

    # -- helpers ---------------------------------------------------------------

    def _delayed(self, spec: np.ndarray, tau: float) -> np.ndarray:
        if tau == 0.0:
            return _fft.irfft(spec, n=self.m)
        return _fft.irfft(spec * np.exp(-2j * np.pi * self.f * float(tau)), n=self.m)

    def _cycles(self, spec_g: np.ndarray, tau: float) -> np.ndarray:
        g = self._delayed(spec_g, tau)
        return 2.0 * np.pi * np.cumsum(g) / self.fs

    def _record(self, x: np.ndarray) -> np.ndarray:
        return x[self.pad : self.pad + self.n]

    # -- the record ------------------------------------------------------------

    def tx_phase(self, tau_s: float = 0.0) -> np.ndarray:
        """Transmit phase error at ``t - tau``, radians, on the record."""
        tau = float(tau_s)
        out = np.zeros(self.n)
        if self._g_tx is not None:
            out += self._record(self._cycles(self._g_tx, tau))
        if self._pm_tx is not None:
            out += self._record(self._delayed(self._pm_tx, tau))
        tt = self.t - tau
        for s, beta, ph in zip(self.imp.spurs, self._betas, self._spur_phase_tx, strict=True):
            if beta:
                out += beta * np.sin(2.0 * np.pi * s.offset_hz * tt + ph)
        return out

    def ref_phase(self) -> np.ndarray:
        """The de-chirp reference's phase error at ``t``, radians, on the record."""
        if getattr(self, "_ref_cache", None) is not None:
            return self._ref_cache
        self._ref_cache = self._ref_phase()
        return self._ref_cache

    def _ref_phase(self) -> np.ndarray:
        if self.imp.single_reference:
            return self.tx_phase(0.0)
        out = np.zeros(self.n)
        if self._g_ref is not None:
            out += self._record(self._cycles(self._g_ref, 0.0))
        if self._pm_ref is not None:
            out += self._record(self._delayed(self._pm_ref, 0.0))
        for s, beta, ph in zip(self.imp.spurs, self._betas, self._spur_phase_ref, strict=True):
            if beta:
                out += beta * np.sin(2.0 * np.pi * s.offset_hz * self.t + ph)
        return out

    def tx_amplitude(self, tau_s: float = 0.0) -> np.ndarray:
        """Radiated amplitude factor at ``t - tau``, ``1 + a``, on the record."""
        out = np.ones(self.n)
        tau = float(tau_s)
        if self._am is not None:
            out += self._record(self._delayed(self._am, tau))
        tt = self.t - tau
        for s, m_am, ph in zip(self.imp.spurs, self._m_am, self._spur_phase_tx, strict=True):
            if m_am:
                out += m_am * np.cos(2.0 * np.pi * s.offset_hz * tt + ph)
        return out

    def residual_phase(self, tau_s: float) -> np.ndarray:
        """``e_ref(t) - e_tx(t - tau)``: what a path of delay ``tau`` leaves."""
        return self.ref_phase() - self.tx_phase(tau_s)


# ----------------------------------------------------------------------------
# The de-chirped phase of the law, nonlinearity included
# ----------------------------------------------------------------------------


def _local_cycles(sweep: SweepSchedule, i: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Cycles accumulated from the start of step ``i`` to local time ``u``."""
    dwell = float(sweep.step_dwell_s)
    d = 1.0 if int(sweep.direction) >= 0 else -1.0
    start = sweep_law.frequency_hz(_replace(sweep, nonlinearity_ppm=0.0), i * dwell)
    k = sweep_law.rate_hz_per_s(sweep, i * dwell)
    c = start * u + d * k * u**2 / 2.0
    eps = float(sweep.nonlinearity_ppm) * 1.0e-6
    if eps != 0.0:
        s = u / dwell
        c = c + d * (k * dwell) * eps * 4.0 * dwell * (s**2 / 2.0 - s**3 / 3.0)
    return c


def local_dechirp_phase_rad(
    sweep: SweepSchedule, t_s: np.ndarray, tau_s: np.ndarray | float
) -> np.ndarray:
    """``phi(t) - phi(t - tau)`` of a stepped law, nonlinearity included, stably.

    Inside a step it is the difference of the step's local cycles, which are at
    most a few times 10^5 and keep full precision; across a step boundary it is
    the difference of :func:`tss.sweep_law.phase_rad`.
    """
    t = np.asarray(t_s, dtype=float)
    tau = np.broadcast_to(np.asarray(tau_s, dtype=float), t.shape)
    te = t - tau
    dwell = float(sweep.step_dwell_s)
    i = sweep_law.step_index(sweep, t)
    same = (sweep_law.step_index(sweep, te) == i) & (te >= 0.0)
    u = t - i * dwell
    out = 2.0 * np.pi * (_local_cycles(sweep, i, u) - _local_cycles(sweep, i, u - tau))
    if not np.all(same):
        sel = ~same
        out = np.array(out, dtype=float)
        out[sel] = sweep_law.phase_rad(sweep, t[sel]) - sweep_law.phase_rad(sweep, te[sel])
    return out


# ----------------------------------------------------------------------------
# Predictions
# ----------------------------------------------------------------------------


def dechirped_pm_density_dbc_hz(
    l_dbc_hz: float, f_hz: np.ndarray | float, tau_s: float, *, single_reference: bool = True
) -> np.ndarray:
    """Residual phase-noise density after de-chirp, per sideband, dBc/Hz.

    A phase noise of single-sideband density ``L`` on the transmitter leaves
    ``L 4 sin^2(pi f tau)`` with a shared reference, and ``2 L`` with an
    independent one of the same specification.
    """
    f = np.asarray(f_hz, dtype=float)
    if single_reference:
        return float(l_dbc_hz) + 10.0 * np.log10(4.0 * np.sin(np.pi * f * float(tau_s)) ** 2)
    return np.full(f.shape, float(l_dbc_hz) + 10.0 * math.log10(2.0))


def oscillator_phase_density_dbc_hz(
    oscillator: Oscillator, carrier_hz: float, f_hz: np.ndarray | float
) -> np.ndarray:
    """``L(f) = (f_c / f)^2 S_y(f) / 2`` of the reference multiplied to ``f_c``, dBc/Hz."""
    h, a = power_law_h(oscillator)
    f = np.asarray(f_hz, dtype=float)
    return 10.0 * np.log10((float(carrier_hz) / f) ** 2 * h * f ** float(a) / 2.0)
