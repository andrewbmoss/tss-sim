"""Layer: L2
WP: WP3
Claims: C4
Status: implemented
Requirements: INS-11, INS-23, INS-24

The receive chain ahead of the digitiser (WP3 Phase 3): front-end compression
and the leakage high-pass.

WHERE THE LEAKAGE GOES
----------------------
The leakage arrives at the receiver input tens of decibels above the noise and
a hundred above the weakest echo.  Two things stand between it and the
digitiser, and this module models both:

1. **The front end** (LNA if any, and the de-chirp mixer) sees the composite of
   leakage, echoes and noise at radio frequency.  Its compression is a
   memoryless function of the composite's envelope.  De-chirping multiplies by
   a constant-envelope reference, which translates the spectrum and leaves the
   envelope unchanged, so the bandpass-equivalent model of the front end
   (AM/AM and AM/PM applied to the complex envelope) can be applied to the
   de-chirped baseband directly.  ``tests/unit/test_receive_chain.py`` checks
   that against a passband simulation of a real memoryless nonlinearity.
   :class:`FrontEnd` is a Rapp limiter specified by its input 1 dB compression
   point, with optional AM/PM.

2. **The leakage high-pass** (INS-11, Rev G section 5.9) sits after the
   de-chirp, in I and Q, and removes the leakage carrier, which de-chirps to a
   fraction of a hertz.  Its corner is programmable and follows the step:
   ``f_c = corner_fraction * 2 K_i R_min / c``, a fixed fraction of the beat of
   the nearest echo INS-04 must hold (100 km), so every step keeps its whole
   echo band (decision of 2026-10-28).  At the A1 settings and the default
   fraction of one half the corner runs from 205 Hz at 0.33 MHz to 1.87 kHz
   above 3 MHz, and a 100 km echo loses 0.017 dB in a fourth-order Butterworth.

   The filter is a cascade of state-variable sections (two integrators per
   biquad, one for an odd order), integrated by the trapezoidal rule.  Its
   states are the integrator outputs, which are the capacitor voltages of a
   real programmable filter, so they carry across a change of corner as they
   would in hardware.  Within a step the filter is time-invariant and is run
   in its modal form with :func:`scipy.signal.lfilter`.

WHAT IS SCORED HERE AND WHAT IS NOT
-----------------------------------
INS-11 is scored against the noise-referenced rule (decision of 2026-10-28,
G5b): the rejection must bring the leakage to 10 dB below the converter's full
scale, which INS-11 places 60 dB above the noise in the digitiser band, so
``required = L/N_fs - 60 + 10``.  Rev G's rule, ``L/E - 60 + 10``, is reported
beside it.  Settling is reported per step against the step's own echo
transition time.  The residual the high-pass cannot remove, the leakage that
spurs, phase noise and AM spread into the echo band, is INS-24's and is scored
in Phase 4.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import signal as _signal

from tss import sweep_law
from tss.constants import C_LIGHT
from tss.interfaces import SweepSchedule

__all__ = [
    "FrontEnd",
    "LeakageHighpass",
    "butterworth_qs",
    "rapp_saturation_amplitude",
    "small_signal_gain_db",
    "compression_image_dbc",
]


# ----------------------------------------------------------------------------
# Front end
# ----------------------------------------------------------------------------


def rapp_saturation_amplitude(input_p1db_w: float, p: float) -> float:
    """Saturation amplitude of a Rapp limiter whose gain is 1 dB down at ``input_p1db_w``.

    ``g(r) = r / (1 + (r/A)^(2p))^(1/(2p))``; the gain ``g/r`` is 1 dB down where
    ``(r/A)^(2p) = 10^(0.1 p) - 1``.
    """
    r1 = math.sqrt(float(input_p1db_w))
    return r1 / (10.0 ** (0.1 * float(p)) - 1.0) ** (1.0 / (2.0 * float(p)))


@dataclass(frozen=True)
class FrontEnd:
    """Bandpass-equivalent compression of the composite at the receiver input.

    Parameters
    ----------
    input_p1db_dbm
        Input 1 dB compression point.  ``None`` is a linear front end.  A
        level-7 diode mixer is about +1 dBm; WP9-2 asks for 10 dB above the
        leakage presented to it.
    rapp_p
        Rapp smoothness.  2 is a soft limiter, large values a hard one.
    am_pm_deg_at_p1db
        Phase conversion at the 1 dB point, degrees; scales as
        ``2 (r/r1)^2 / (1 + (r/r1)^2)`` of it.
    am_am
        Optional replacement AM/AM function of the envelope amplitude, in
        root-watts, for example a describing function.  Overrides the Rapp form.
    """

    input_p1db_dbm: float | None = None
    rapp_p: float = 2.0
    am_pm_deg_at_p1db: float = 0.0
    am_am: Any = None

    @property
    def linear(self) -> bool:
        return self.input_p1db_dbm is None and self.am_am is None

    @property
    def input_p1db_w(self) -> float:
        return 10.0 ** ((float(self.input_p1db_dbm) - 30.0) / 10.0)

    def gain(self, r: np.ndarray) -> np.ndarray:
        """Complex gain at envelope amplitude ``r`` (root-watts)."""
        r = np.asarray(r, dtype=float)
        if self.am_am is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                mag = np.where(r > 0.0, np.asarray(self.am_am(r), dtype=float) / np.where(r > 0, r, 1.0), 1.0)
        elif self.input_p1db_dbm is None:
            return np.ones(r.shape, dtype=complex)
        else:
            a_sat = rapp_saturation_amplitude(self.input_p1db_w, self.rapp_p)
            mag = (1.0 + (r / a_sat) ** (2.0 * self.rapp_p)) ** (-1.0 / (2.0 * self.rapp_p))
        phase = np.zeros(r.shape)
        if self.am_pm_deg_at_p1db and self.input_p1db_dbm is not None:
            u = (r / math.sqrt(self.input_p1db_w)) ** 2
            phase = math.radians(self.am_pm_deg_at_p1db) * 2.0 * u / (1.0 + u)
        return mag * np.exp(1j * phase)

    def apply(self, x: np.ndarray) -> np.ndarray:
        """The composite through the front end.  Unity small-signal gain."""
        x = np.asarray(x, dtype=complex)
        if self.linear:
            return x.copy()
        return x * self.gain(np.abs(x))

    def describe(self) -> Mapping[str, Any]:
        return {
            "input_p1db_dbm": self.input_p1db_dbm,
            "rapp_p": self.rapp_p,
            "am_pm_deg_at_p1db": self.am_pm_deg_at_p1db,
            "custom_am_am": self.am_am is not None,
            "model": "linear" if self.linear else "bandpass-equivalent envelope (Rapp AM/AM, AM/PM)",
        }

    @classmethod
    def from_config(cls, cfg: Mapping[str, Any] | None) -> FrontEnd:
        c = dict(cfg or {})
        p1 = c.get("input_p1db_dbm")
        return cls(
            input_p1db_dbm=None if p1 is None else float(p1),
            rapp_p=float(c.get("rapp_smoothness", 2.0)),
            am_pm_deg_at_p1db=float(c.get("am_pm_deg_at_p1db", 0.0)),
        )


def _linearised(front_end: FrontEnd, carrier_power_w: float, eps: float = 1e-5) -> tuple[complex, complex]:
    """``(a, b)`` of a weak signal ``d`` beside a carrier of amplitude ``R``: out = a d + b conj(d).

    With ``F(r)`` the complex output for a real input ``r``, a weak signal sees
    the slope ``F'(R)`` on its in-phase part and the chord ``F(R)/R`` on its
    quadrature part, so ``a = (F' + F/R) / 2`` and ``b = (F' - F/R) / 2``.
    """
    r = math.sqrt(float(carrier_power_w))

    def big_f(x: float) -> complex:
        return complex(front_end.gain(np.array([x]))[0]) * x

    slope = (big_f(r * (1.0 + eps)) - big_f(r * (1.0 - eps))) / (2.0 * eps * r)
    chord = big_f(r) / r
    return 0.5 * (slope + chord), 0.5 * (slope - chord)


def small_signal_gain_db(front_end: FrontEnd, carrier_power_w: float) -> float:
    """Gain of a weak signal beside a strong carrier, dB: ``|a|^2``.

    This is the desensitisation of an echo by the leakage.  The rest of the
    weak signal's power goes into its image (:func:`compression_image_dbc`).
    """
    if front_end.linear or carrier_power_w <= 0.0:
        return 0.0
    a, _ = _linearised(front_end, carrier_power_w)
    return 10.0 * math.log10(abs(a) ** 2)


def compression_image_dbc(front_end: FrontEnd, carrier_power_w: float) -> float:
    """The image a weak signal acquires beside a compressing carrier, dBc: ``|b/a|^2``.

    After de-chirping the carrier is the leakage near zero beat, so the image of
    an echo at beat ``f_b`` appears at ``-f_b``: a ghost on the other side of
    the spectrum, at this level below the echo.
    """
    if front_end.linear or carrier_power_w <= 0.0:
        return float("-inf")
    a, b = _linearised(front_end, carrier_power_w)
    return 20.0 * math.log10(max(abs(b), 1e-300) / abs(a))


# ----------------------------------------------------------------------------
# The leakage high-pass
# ----------------------------------------------------------------------------


def butterworth_qs(order: int) -> list[float]:
    """Q of each second-order section of a Butterworth filter of ``order``."""
    n = int(order)
    return [1.0 / (2.0 * math.sin((2 * k - 1) * math.pi / (2 * n))) for k in range(1, n // 2 + 1)]


def _continuous_ss(order: int, wc: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """State space of the high-pass cascade in integrator (physical) states."""
    blocks = []
    if order % 2:
        blocks.append((np.array([[-wc]]), np.array([[wc]]), np.array([[-1.0]]), 1.0))
    for q in butterworth_qs(order):
        a = wc * np.array([[-1.0 / q, -1.0], [1.0, 0.0]])
        b = wc * np.array([[1.0], [0.0]])
        c = np.array([[-1.0 / q, -1.0]])
        blocks.append((a, b, c, 1.0))
    a_t, b_t, c_t, d_t = blocks[0]
    for a2, b2, c2, d2 in blocks[1:]:
        n1, n2 = a_t.shape[0], a2.shape[0]
        a_n = np.zeros((n1 + n2, n1 + n2))
        a_n[:n1, :n1] = a_t
        a_n[n1:, :n1] = b2 @ c_t
        a_n[n1:, n1:] = a2
        b_n = np.vstack([b_t, b2 * d_t])
        c_n = np.hstack([d2 * c_t, c2])
        a_t, b_t, c_t, d_t = a_n, b_n, c_n, d2 * d_t
    return a_t, b_t, c_t, float(d_t)


@dataclass(frozen=True)
class LeakageHighpass:
    """Programmable high-pass after the de-chirp, in I and Q.

    Parameters
    ----------
    order
        Butterworth order.  The prototype board (TV BOM A18) assumes four.
    corner_fraction
        Corner as a fraction of each step's beat at ``r_min_m``.
    r_min_m
        Nearest echo range the step must keep, INS-04's 100 km.
    fixed_corner_hz
        A fixed corner instead of the per-step rule, for comparison.
    """

    order: int = 4
    corner_fraction: float = 0.5
    r_min_m: float = 1.0e5
    fixed_corner_hz: float | None = None

    def __post_init__(self) -> None:
        if self.order < 1:
            raise ValueError("order must be at least 1")
        if not (self.corner_fraction > 0.0 and self.r_min_m > 0.0):
            raise ValueError("corner_fraction and r_min_m must be positive")

    def corner_hz(self, rate_hz_per_s: np.ndarray | float) -> np.ndarray:
        """Corner for a step of the given sweep rate."""
        k = np.abs(np.asarray(rate_hz_per_s, dtype=float))
        if self.fixed_corner_hz is not None:
            return np.full(k.shape, float(self.fixed_corner_hz))
        return float(self.corner_fraction) * 2.0 * k * float(self.r_min_m) / C_LIGHT

    def response(self, f_hz: np.ndarray | float, corner_hz: float) -> np.ndarray:
        """Complex analog response at ``f_hz`` for a corner ``corner_hz``."""
        b, a = _signal.butter(self.order, 2.0 * np.pi * float(corner_hz), btype="high", analog=True)
        _, h = _signal.freqs(b, a, worN=2.0 * np.pi * np.atleast_1d(np.abs(np.asarray(f_hz, dtype=float))))
        return h

    def response_db(self, f_hz: np.ndarray | float, corner_hz: float) -> np.ndarray:
        return 20.0 * np.log10(np.abs(self.response(f_hz, corner_hz)))

    def _discrete(self, corner_hz: float, fs: float):
        """Trapezoidal integration of the physical states, in modal form.

        ``x[n+1] = P x[n] + Q (u[n] + u[n+1])``, ``y[n] = C x[n] + D u[n]``, with
        ``P = (I - hA/2)^-1 (I + hA/2)`` and ``Q = (I - hA/2)^-1 hB/2``.  The state
        is the integrator outputs themselves, so it can be carried across a
        change of corner.  The corner is prewarped.
        """
        wc = 2.0 * fs * math.tan(math.pi * float(corner_hz) / fs)
        a, b, c, d = _continuous_ss(self.order, wc)
        h = 1.0 / fs
        eye = np.eye(a.shape[0])
        m = np.linalg.inv(eye - 0.5 * h * a)
        p = m @ (eye + 0.5 * h * a)
        q = (m @ (0.5 * h * b)).ravel()
        lam, v = np.linalg.eig(p)
        vinv = np.linalg.inv(v)
        return lam, v, vinv, vinv @ q, (c @ v).ravel(), float(d)

    def segments(self, sweep: SweepSchedule, t0_s: float, n: int, fs: float) -> list[tuple[int, int, float]]:
        """``(start, stop, corner)`` for each run of constant corner in the record."""
        t = float(t0_s) + np.arange(int(n)) / float(fs)
        if sweep_law.is_stepped(sweep):
            idx = sweep_law.step_index(sweep, t)
            rates = sweep_law.rate_hz_per_s(sweep, t)
        else:
            idx = np.zeros(t.shape, dtype=int)
            rates = np.full(t.shape, abs(float(sweep.k_f_hz_per_s)))
        cuts = np.flatnonzero(np.diff(idx)) + 1
        starts = np.concatenate([[0], cuts])
        stops = np.concatenate([cuts, [t.size]])
        return [(int(s), int(e), float(self.corner_hz(rates[s]))) for s, e in zip(starts, stops, strict=True)]

    def apply(
        self,
        x: np.ndarray,
        sweep: SweepSchedule,
        fs_hz: float,
        *,
        t0_s: float = 0.0,
        state: np.ndarray | None = None,
        return_state: bool = False,
    ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        """Filter a de-chirped record, corner switched at each step boundary.

        I and Q are filtered by identical real filters, which for a real filter
        is the same as filtering the complex stream.  ``state`` is the physical
        state vector (complex, I + jQ) to start from; zero by default.
        """
        x = np.asarray(x, dtype=complex)
        fs = float(fs_hz)
        n_state = self.order
        xs = np.zeros(n_state, dtype=complex) if state is None else np.asarray(state, dtype=complex).copy()
        y = np.empty_like(x)
        cache: dict[float, Any] = {}
        for s, e, fc in self.segments(sweep, t0_s, x.size, fs):
            key = round(fc, 9)
            if key not in cache:
                cache[key] = self._discrete(fc, fs)
            lam, v, vinv, qm, cm, dd = cache[key]
            u = x[s:e]
            u_next = x[s + 1 : e + 1] if e < x.size else np.concatenate([x[s + 1 : e], x[-1:]])
            vsum = u + u_next
            z0 = vinv @ xs
            out = dd * u
            z_end = np.empty(n_state, dtype=complex)
            for j in range(n_state):
                # z_j[n+1] = lam z_j[n] + q_j (u[n] + u[n+1]); lfilter returns z_j[n+1].
                w = _signal.lfilter([qm[j]], [1.0, -lam[j]], vsum, zi=np.array([lam[j] * z0[j]]))[0]
                zj = np.concatenate([[z0[j]], w[:-1]])
                out = out + cm[j] * zj
                z_end[j] = w[-1]
            y[s:e] = out
            xs = v @ z_end
        if return_state:
            return y, xs
        return y

    def steady_state(
        self, u0: complex, sweep: SweepSchedule, fs_hz: float, *, t0_s: float = 0.0
    ) -> np.ndarray:
        """The state after a long run on a constant input ``u0``, at the corner in force at ``t0_s``.

        ``x* = (I - P)^-1 2 Q u0``, solved in modal form.  It is the state a
        filter that has been running on the leakage carrier is in when the
        record starts: the transmitter runs between sweeps (WP3 Phase 4
        operating assumption, 2026-10-28), so the de-chirped leakage at the
        start of a sweep is a carrier the filter has already absorbed.  The
        output from it is ``D u0 + C x* = 0`` for a high-pass, which is the
        check.
        """
        fs = float(fs_hz)
        seg = self.segments(sweep, t0_s, 1, fs)[0]
        lam, v, _vinv, qm, _cm, _dd = self._discrete(seg[2], fs)
        z = 2.0 * qm * complex(u0) / (1.0 - lam)
        return v @ z

    def describe(self) -> Mapping[str, Any]:
        return {
            "order": self.order,
            "corner_fraction": self.corner_fraction,
            "r_min_m": self.r_min_m,
            "fixed_corner_hz": self.fixed_corner_hz,
            "realisation": "Butterworth, state-variable sections, trapezoidal; states carried across corner changes",
        }

    @classmethod
    def from_config(cls, cfg: Mapping[str, Any] | None) -> LeakageHighpass:
        c = dict(cfg or {})
        fixed = c.get("fixed_corner_hz")
        return cls(
            order=int(c.get("order", 4)),
            corner_fraction=float(c.get("corner_fraction", 0.5)),
            r_min_m=float(c.get("r_min_m", 1.0e5)),
            fixed_corner_hz=None if fixed is None else float(fixed),
        )
