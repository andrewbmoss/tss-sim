"""Layer: L2
WP: WP3
Claims: C4
Status: partial
Requirements: INS-21, INS-23, INS-24

Transmit-to-receive leakage, as de-chirped complex baseband (WP3 Phase 1).

WHY THIS IS THE HARDEST PART OF THE INSTRUMENT
----------------------------------------------
The sounder transmits continuously while it receives, so the transmitter is the
strongest signal at the receiver input at all times.  The echo arrives after
free-space spreading over twice the range and after ionospheric absorption,
while the leakage takes a path of a metre to a few tens of metres.  At INS-23's
40 dB isolation the leakage lies 97 dB above the free-space echo from INS-04's
1050 km at 5 MHz, and 120.5 dB above it at 1.26 MHz, where antenna efficiency
is low (:mod:`tss.l2_instrument.reference_echo`; Rev G section 5.9 as amended
2026-10-27).  The rest has to come from the de-chirp.

De-chirping is what makes it survivable.  Leakage at a delay of a few
nanoseconds de-chirps to a beat of a fraction of a hertz, while an echo at
600 km de-chirps to kilohertz, so the two separate in the transform even though
they do not separate in time.  INS-24 is a requirement on what is left after
that separation.

WHAT THIS MODULE PRODUCES
-------------------------
:func:`leakage_signal` returns the leakage at the receiver input in the same
units and with the same phase convention as the L2 baseline's echoes: complex
samples in root-watts, de-chirped against the reference, with the phase
``phi_ref(t) - phi_rx(t)``.  For a coupling path of delay ``tau`` and complex
gain ``H(f)`` the received leakage carries the transmitted phase at ``t - tau``
plus the excess phase of ``H``, so

    x(t) = sqrt(P_tx) conj(H(f_e)) exp(j [Dphi(t) + e_ref(t) - e_tx(t - tau)])

with ``f_e = f(t - tau)`` the frequency that is arriving, ``Dphi`` the
de-chirped phase of the ideal law and ``e_tx`` and ``e_ref`` the phase errors
of the transmitted chirp and of the de-chirp reference.  Under INS-21 they are
the same source, ``e_ref = e_tx``, and the residual error is the lag-``tau``
difference of one error; with an independent reference it is not, and the
range-correlation credit is lost.  Phase 1 carries the single-tone frequency
error of Rev G section 5.9 (:class:`SingleToneFM`) so that the analytic forms
can be tested; WP3 Phase 2 replaces it with the oscillator, spurs, AM noise and
sweep nonlinearity.

``Dphi`` is computed as the integral of the transmitted frequency over
``[t - tau, t]`` in closed form inside a step rather than as the difference of
two accumulated phases, which at 4 s into a sweep are 10^8 radians apart from
the answer and would cost the residual seven digits.  Where ``t - tau`` falls
in the previous step the exact difference of :func:`tss.sweep_law.phase_rad`
is used, so the step transition is in the samples as it is for the echoes.

THE ISOLATION PROFILE
---------------------
WP3 decision G1 (2026-10-27): the scored isolation profiles are parametric
brackets, with one profile anchored by a moment-method model of the antennas
(``data/reference/wp3_isolation_mom``).  :class:`LeakageModel` composes

* a base isolation, flat (``isolation_db``) or tabulated against frequency
  (:func:`mom_anchor_model` loads the moment-method anchor);
* a sinusoidal ripple in decibels, ``isolation_ripple_db`` peak to peak with
  period ``ripple_period_hz`` in frequency, which the sweep turns into
  amplitude modulation at ``K / ripple_period_hz``;
* any number of :class:`Resonance` terms, each a parallel resonant path that
  lowers the isolation by ``depth_db`` at ``f0_hz`` and turns the coupling phase
  through it;
* a retuning transient: when the match retunes at a step boundary the coupling
  departs by ``retune_transient_db`` and settles with ``retune_time_constant_s``.

The coupling delay is ``coupling_delay_s``, with an optional slow random
jitter (``delay_jitter_s`` rms, first-order Gauss-Markov with bandwidth
``jitter_bandwidth_hz``) for a moving deployed structure.  WP3 decision G2
(2026-10-27) brackets the path at 1, 3, 10 and 30 m (:data:`COUPLING_PATHS_M`);
the baseline stays at 3 ns (Open Item 13).

WHAT IS NOT HERE YET
--------------------
The leakage is not yet added to the baseline's samples, and the receiver's
compression and the leakage high-pass are not modelled: both are WP3 Phase 3.
:func:`dechirped_residual` and :func:`score_ins23_ins24` (WP3 Phase 4) measure
the residual in each step's echo band and score INS-23 and INS-24.  The
companion canceller ``tss.l3_signal.leakage_suppress`` is WP5; Phase 4 reports
an ideal-canceller bound instead.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from dataclasses import replace as _replace
from typing import Any

import numpy as np
from scipy import signal as _signal
from scipy import special as _special

from tss import sweep_law
from tss.constants import C_LIGHT
from tss.interfaces import SweepSchedule
from tss.l2_instrument.impairments import (
    ImpairmentRealization,
    TransmitImpairments,
    local_dechirp_phase_rad,
)

__all__ = [
    "BASELINE_ISOLATION_DB",
    "BASELINE_COUPLING_DELAY_S",
    "ISOLATION_BRACKETS_DB",
    "COUPLING_PATHS_M",
    "Resonance",
    "SingleToneFM",
    "LeakageModel",
    "isolation_profile_db",
    "coupling_response",
    "dechirp_phase_rad",
    "leakage_signal",
    "leakage_over_echo_db",
    "range_correlation_suppression_db",
    "fm_spur_dbc",
    "dechirped_fm_sideband_dbc",
    "small_angle_sideband_dbc",
    "ripple_sideband_dbc",
    "ripple_modulation_hz",
    "jitter_phase_rms_rad",
    "MOM_ANCHOR_RELPATH",
    "MOM_BUS_RELPATH",
    "coupled_resonator_isolation_db",
    "mom_anchor_table",
    "mom_anchor_model",
    "dechirped_residual",
    "score_ins23_ins24",
    "INS24_MARGIN_DB",
    "INS23_THRESHOLD_DB",
]

#: INS-23 threshold, dB.  A requirement, not a measurement.
BASELINE_ISOLATION_DB = 40.0

#: Coupling delay of the baseline configuration, seconds.  Open Item 13: the
#: deployed geometry is not fixed.  WP3 decision G2 keeps 3 ns as the baseline
#: and brackets the path at :data:`COUPLING_PATHS_M`.
BASELINE_COUPLING_DELAY_S = 3.0e-9

#: WP3 decision G1(a): the scored flat isolation brackets, dB.
ISOLATION_BRACKETS_DB = (30.0, 40.0, 50.0, 60.0)

#: WP3 decision G2: the coupling-path bracket, m (3.3 to 100 ns).
COUPLING_PATHS_M = (1.0, 3.0, 10.0, 30.0)

_LN10_OVER_20 = math.log(10.0) / 20.0


# ----------------------------------------------------------------------------
# The coupling path
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class Resonance:
    """A parallel resonant coupling path.

    The coupling gain is multiplied by ``1 + (g - 1) / (1 + j x)`` with
    ``g = 10^(depth_db / 20)`` and ``x = 2 q (f - f0) / f0``.  At ``f0`` the
    isolation falls by ``depth_db``; far from it the factor is one.  Across the
    resonance the coupling phase turns, which is what a structural resonance
    does to the leakage and what a flat-isolation model omits.
    """

    f0_hz: float
    q: float
    depth_db: float

    def __post_init__(self) -> None:
        if not (self.f0_hz > 0.0 and self.q > 0.0):
            raise ValueError("a resonance needs f0_hz > 0 and q > 0")
        if not math.isfinite(self.depth_db):
            raise ValueError("depth_db must be finite")

    def response(self, f_hz: np.ndarray | float) -> np.ndarray:
        f = np.asarray(f_hz, dtype=float)
        g = 10.0 ** (float(self.depth_db) / 20.0)
        x = 2.0 * float(self.q) * (f - float(self.f0_hz)) / float(self.f0_hz)
        return 1.0 + (g - 1.0) / (1.0 + 1j * x)


@dataclass(frozen=True)
class SingleToneFM:
    """Transmit frequency error ``deviation_hz cos(2 pi rate_hz t + phase_rad)``.

    The error of Rev G section 5.9.  Its phase is
    ``(deviation_hz / rate_hz) sin(2 pi rate_hz t + phase_rad)``.
    """

    deviation_hz: float
    rate_hz: float
    phase_rad: float = 0.0

    def __post_init__(self) -> None:
        if not self.rate_hz > 0.0:
            raise ValueError("rate_hz must be positive")

    def __call__(self, t_s: np.ndarray) -> np.ndarray:
        t = np.asarray(t_s, dtype=float)
        beta = float(self.deviation_hz) / float(self.rate_hz)
        return beta * np.sin(2.0 * np.pi * float(self.rate_hz) * t + float(self.phase_rad))


@dataclass(frozen=True)
class LeakageModel:
    """The coupling path from transmitter to receiver.

    Parameters
    ----------
    isolation_db
        Base isolation, positive decibels, when no table is given.  Default
        :data:`BASELINE_ISOLATION_DB`, the INS-23 threshold.
    coupling_delay_s
        One-way delay of the coupling path.  Default
        :data:`BASELINE_COUPLING_DELAY_S` (Open Item 13).
    isolation_ripple_db, ripple_period_hz, ripple_phase_rad
        A sinusoidal departure of the isolation from its base, in decibels,
        peak to peak, with the given period in frequency.  A ripple needs a
        period.
    resonances
        :class:`Resonance` terms applied on top of the base and the ripple.
    table_f_hz, table_isolation_db
        A tabulated base isolation, interpolated linearly in decibels and held
        at the end values outside the table.  Replaces ``isolation_db``.
    delay_jitter_s, jitter_bandwidth_hz
        RMS and bandwidth of a first-order Gauss-Markov variation of the
        coupling delay.  Zero jitter draws nothing from the generator.
    retune_transient_db, retune_time_constant_s
        Departure of the coupling at each step boundary of a stepped sweep, in
        decibels (positive is less isolation), decaying exponentially.
    label
        Free text carried into :meth:`describe`.
    """

    isolation_db: float = BASELINE_ISOLATION_DB
    coupling_delay_s: float = BASELINE_COUPLING_DELAY_S
    isolation_ripple_db: float = 0.0
    ripple_period_hz: float | None = None
    ripple_phase_rad: float = 0.0
    resonances: tuple[Resonance, ...] = ()
    table_f_hz: tuple[float, ...] | None = None
    table_isolation_db: tuple[float, ...] | None = None
    delay_jitter_s: float = 0.0
    jitter_bandwidth_hz: float = 1.0
    retune_transient_db: float = 0.0
    retune_time_constant_s: float = 1.0e-3
    label: str = ""

    def __post_init__(self) -> None:
        if not math.isfinite(self.isolation_db):
            raise ValueError("isolation_db must be finite")
        if not self.coupling_delay_s >= 0.0:
            raise ValueError("coupling_delay_s must not be negative")
        if self.isolation_ripple_db < 0.0:
            raise ValueError("isolation_ripple_db is peak to peak and must not be negative")
        if self.isolation_ripple_db > 0.0 and not (self.ripple_period_hz and self.ripple_period_hz > 0.0):
            raise ValueError("an isolation ripple needs a positive ripple_period_hz")
        if (self.table_f_hz is None) != (self.table_isolation_db is None):
            raise ValueError("table_f_hz and table_isolation_db go together")
        if self.table_f_hz is not None:
            f = np.asarray(self.table_f_hz, dtype=float)
            iso = np.asarray(self.table_isolation_db, dtype=float)
            if f.shape != iso.shape or f.size < 2:
                raise ValueError("an isolation table needs two or more parallel entries")
            if np.any(np.diff(f) <= 0.0):
                raise ValueError("table_f_hz must be strictly increasing")
            if not np.all(np.isfinite(iso)):
                raise ValueError("table_isolation_db must be finite")
            object.__setattr__(self, "table_f_hz", tuple(float(v) for v in f))
            object.__setattr__(self, "table_isolation_db", tuple(float(v) for v in iso))
        object.__setattr__(self, "resonances", tuple(self.resonances))
        if self.delay_jitter_s < 0.0 or self.jitter_bandwidth_hz <= 0.0:
            raise ValueError("delay_jitter_s must be >= 0 and jitter_bandwidth_hz > 0")
        if self.retune_time_constant_s <= 0.0:
            raise ValueError("retune_time_constant_s must be positive")

    # -- construction ------------------------------------------------------

    @classmethod
    def from_coupling_path(cls, path_m: float, **kwargs: Any) -> LeakageModel:
        """A model whose delay is a free-space path of ``path_m`` metres."""
        return cls(coupling_delay_s=float(path_m) / C_LIGHT, **kwargs)

    @classmethod
    def from_table(
        cls, f_hz: Sequence[float], isolation_db: Sequence[float], **kwargs: Any
    ) -> LeakageModel:
        """A model on a tabulated isolation, for a modelled or measured profile."""
        return cls(table_f_hz=tuple(f_hz), table_isolation_db=tuple(isolation_db), **kwargs)

    @classmethod
    def from_config(cls, isolation: Mapping[str, Any], **overrides: Any) -> LeakageModel:
        """From the ``isolation`` block of an instrument configuration.

        Reads ``transmit_receive_db`` and ``coupling_delay_s`` (the keys of
        ``config/instrument/baseline_6u.yaml``) and the optional profile keys of
        this class under their own names.  Unknown keys are ignored.
        """
        cfg = dict(isolation or {})
        kw: dict[str, Any] = {
            "isolation_db": float(cfg.get("transmit_receive_db", BASELINE_ISOLATION_DB)),
            "coupling_delay_s": float(cfg.get("coupling_delay_s", BASELINE_COUPLING_DELAY_S)),
        }
        for key in (
            "isolation_ripple_db", "ripple_period_hz", "ripple_phase_rad", "delay_jitter_s",
            "jitter_bandwidth_hz", "retune_transient_db", "retune_time_constant_s",
        ):
            if cfg.get(key) is not None:
                kw[key] = float(cfg[key])
        if cfg.get("resonances"):
            kw["resonances"] = tuple(Resonance(**dict(r)) for r in cfg["resonances"])
        kw.update(overrides)
        return cls(**kw)

    # -- properties --------------------------------------------------------

    @property
    def coupling_path_m(self) -> float:
        return float(self.coupling_delay_s) * C_LIGHT

    def describe(self) -> Mapping[str, Any]:
        out = asdict(self)
        out["resonances"] = [asdict(r) for r in self.resonances]
        out["coupling_path_m"] = self.coupling_path_m
        out["profile"] = "tabulated" if self.table_f_hz is not None else "flat"
        if self.table_f_hz is not None:
            out["table_f_hz"] = [self.table_f_hz[0], self.table_f_hz[-1], len(self.table_f_hz)]
            out["table_isolation_db"] = [min(self.table_isolation_db), max(self.table_isolation_db)]
        return out


def _base_isolation_db(f: np.ndarray, model: LeakageModel) -> np.ndarray:
    if model.table_f_hz is None:
        return np.full(f.shape, float(model.isolation_db))
    return np.interp(f, np.asarray(model.table_f_hz), np.asarray(model.table_isolation_db))


def _ripple_db(f: np.ndarray, model: LeakageModel) -> np.ndarray:
    if model.isolation_ripple_db == 0.0:
        return np.zeros(f.shape)
    arg = 2.0 * np.pi * f / float(model.ripple_period_hz) + float(model.ripple_phase_rad)
    return 0.5 * float(model.isolation_ripple_db) * np.sin(arg)


def coupling_response(f_hz: np.ndarray | float, model: LeakageModel) -> np.ndarray:
    """Complex voltage gain of the coupling path, less its pure delay.

    ``|H|^2`` is ``10^(-isolation/10)``; the phase is the excess phase of the
    resonances.  The delay is carried separately, by the transmitted law at
    ``t - tau``.
    """
    f = np.asarray(f_hz, dtype=float)
    mag = 10.0 ** (-(_base_isolation_db(f, model) + _ripple_db(f, model)) / 20.0)
    h = mag.astype(complex)
    for r in model.resonances:
        h = h * r.response(f)
    return h


def isolation_profile_db(f_hz: np.ndarray | float, model: LeakageModel) -> np.ndarray:
    """Isolation against frequency, positive decibels, parallel to ``f_hz``.

    INS-23 is scored on this array, and its departures from flat are what put
    leakage sidebands into the echo band, so INS-23 and INS-24 share it.
    """
    h = coupling_response(f_hz, model)
    return -20.0 * np.log10(np.abs(h))


# ----------------------------------------------------------------------------
# The de-chirped phase of a short delay
# ----------------------------------------------------------------------------


def dechirp_phase_rad(
    sweep: SweepSchedule, t_s: np.ndarray, tau_s: np.ndarray | float
) -> np.ndarray:
    """``phi(t) - phi(t - tau)`` for the transmitted law, computed stably.

    Inside a step it is the difference of the step's local cycles
    (:func:`tss.l2_instrument.impairments.local_dechirp_phase_rad`), which keeps
    full precision for a nanosecond delay with or without the law's
    nonlinearity; for a uniform sweep it is the closed form
    ``2 pi [f(t - tau) tau + d K tau^2 / 2]``.  Where ``t - tau`` lies in the
    previous step it is the difference of :func:`tss.sweep_law.phase_rad`.
    """
    t = np.asarray(t_s, dtype=float)
    tau = np.broadcast_to(np.asarray(tau_s, dtype=float), t.shape)
    te = t - tau
    if sweep_law.is_stepped(sweep):
        # Local cycles inside a step, nonlinearity included (WP3 Phase 2).
        return local_dechirp_phase_rad(sweep, t, tau)
    d = 1.0 if int(sweep.direction) >= 0 else -1.0
    k = float(sweep.k_f_hz_per_s)
    f0 = float(sweep.f_start_hz) if d > 0 else float(sweep.f_stop_hz)
    if float(sweep.nonlinearity_ppm) != 0.0:
        raise NotImplementedError(
            "a nonlinear uniform sweep is WP3 Phase 2; the stepped law carries the nonlinearity"
        )
    return 2.0 * np.pi * ((f0 + d * k * te) * tau + d * k * tau**2 / 2.0)


def _tx_frequency_hz(sweep: SweepSchedule, t: np.ndarray) -> np.ndarray:
    if sweep_law.is_stepped(sweep):
        return sweep_law.frequency_hz(sweep, t)
    d = 1.0 if int(sweep.direction) >= 0 else -1.0
    f0 = float(sweep.f_start_hz) if d > 0 else float(sweep.f_stop_hz)
    return f0 + d * float(sweep.k_f_hz_per_s) * t


def _sweep_duration_s(sweep: SweepSchedule) -> float:
    if sweep_law.is_stepped(sweep):
        return sweep_law.active_duration_s(sweep)
    return abs(float(sweep.f_stop_hz) - float(sweep.f_start_hz)) / abs(float(sweep.k_f_hz_per_s))


def _delay_trace(
    model: LeakageModel, n: int, fs_hz: float, rng: np.random.Generator
) -> np.ndarray:
    tau0 = float(model.coupling_delay_s)
    sigma = float(model.delay_jitter_s)
    if sigma == 0.0:
        return np.full(n, tau0)
    a = math.exp(-2.0 * math.pi * float(model.jitter_bandwidth_hz) / float(fs_hz))
    w = rng.normal(0.0, 1.0, n)
    x = _signal.lfilter([math.sqrt(1.0 - a * a)], [1.0, -a], w, zi=[a * rng.normal()])[0]
    return tau0 + sigma * x


# ----------------------------------------------------------------------------
# The leakage
# ----------------------------------------------------------------------------


def leakage_signal(
    sweep: SweepSchedule,
    fs_hz: float,
    model: LeakageModel,
    rng: np.random.Generator,
    *,
    n_samples: int | None = None,
    t0_s: float = 0.0,
    transmit_power_w: float = 1.0,
    tx_chain_efficiency_db: float = 0.0,
    tx_phase_error: Callable[[np.ndarray], np.ndarray] | None = None,
    reference_phase_error: str | Callable[[np.ndarray], np.ndarray] | None = "common",
    impairments: TransmitImpairments | None = None,
    realization: ImpairmentRealization | None = None,
    primed: bool = False,
) -> np.ndarray:
    """The leakage at the receiver input, de-chirped, as complex root-watts.

    Samples are at ``t0_s + n / fs_hz``.  ``n_samples`` defaults to the rest of
    the sweep from ``t0_s``.  Before ``tau`` nothing has been transmitted and the
    sample is zero, as for the echoes.

    ``tx_phase_error`` is a function of time giving the phase error of the
    transmitted chirp, radians (for example :class:`SingleToneFM`).
    ``reference_phase_error`` is ``"common"`` for INS-21 (the de-chirp
    reference carries the same error), ``None`` for an independent reference
    with no error of its own, or a function of time.

    ``primed`` (WP3 Phase 4 default, decision of 2026-10-28) says the
    transmitter was already running when the record starts, so there is no
    switch-on at ``tau``.

    ``impairments`` (WP3 Phase 2) applies the transmitter's oscillator noise,
    phase-noise floor, AM noise, spurs and sweep nonlinearity, with INS-21 as
    its ``single_reference``.  ``realization`` is the draw to use, so that the
    leakage and the echoes of one record share one oscillator; without it a
    draw is made from ``rng`` for this record.  With ``impairments`` given, the
    law's ``nonlinearity_ppm`` is used only if ``impairments.sweep_nonlinearity``
    is on; without it the sweep is used as given.

    It is meant to be added to the echoes before the receiver is applied, not
    after, because whether the front end compresses is decided by the sum
    (WP3 Phase 3).  Randomness comes from ``rng`` only, and only when the model
    has delay jitter.
    """
    fs = float(fs_hz)
    if not fs > 0.0:
        raise ValueError("fs_hz must be positive")
    if transmit_power_w <= 0.0:
        raise ValueError("transmit_power_w must be positive")
    if n_samples is None:
        n_samples = int(np.floor((_sweep_duration_s(sweep) - float(t0_s)) * fs))
    n = int(n_samples)
    if n < 1:
        raise ValueError("no samples to generate")
    t = float(t0_s) + np.arange(n, dtype=float) / fs

    tau = _delay_trace(model, n, fs, rng)
    te = t - tau
    # Primed: the transmitter was already running, so leakage is present from
    # the first sample (the law extended back before the sweep starts).
    live = np.ones(t.shape, dtype=bool) if primed else te >= 0.0
    law = sweep
    if impairments is not None and not impairments.sweep_nonlinearity:
        law = _replace(sweep, nonlinearity_ppm=0.0)
    dphi = dechirp_phase_rad(law, t, tau)
    if (
        impairments is not None
        and impairments.sweep_nonlinearity
        and not impairments.single_reference
        and float(sweep.nonlinearity_ppm) != 0.0
    ):
        # An independent reference sweeps the ideal law; the transmitter does not.
        ideal = _replace(sweep, nonlinearity_ppm=0.0)
        dphi = dphi + (sweep_law.phase_rad(ideal, t) - sweep_law.phase_rad(sweep, t))
    f_e = _tx_frequency_hz(law, np.clip(te, 0.0, None))

    phase = dphi
    if tx_phase_error is not None:
        e_tx_delayed = np.asarray(tx_phase_error(te), dtype=float)
        if isinstance(reference_phase_error, str):
            if reference_phase_error != "common":
                raise ValueError("reference_phase_error must be 'common', None or a callable")
            e_ref = np.asarray(tx_phase_error(t), dtype=float)
        elif reference_phase_error is None:
            e_ref = np.zeros_like(t)
        else:
            e_ref = np.asarray(reference_phase_error(t), dtype=float)
        phase = phase + (e_ref - e_tx_delayed)
    elif callable(reference_phase_error):
        phase = phase + np.asarray(reference_phase_error(t), dtype=float)

    h = coupling_response(f_e, model)
    if model.retune_transient_db != 0.0 and sweep_law.is_stepped(sweep):
        dwell = float(sweep.step_dwell_s)
        u = te - sweep_law.step_index(sweep, np.clip(te, 0.0, None)) * dwell
        excess_db = float(model.retune_transient_db) * np.exp(-u / float(model.retune_time_constant_s))
        h = h * 10.0 ** (excess_db / 20.0)

    amp: np.ndarray | float = 1.0
    if impairments is not None and impairments.needs_realization:
        tau0 = float(model.coupling_delay_s)
        if realization is None:
            realization = impairments.realize(
                sweep, fs, n, rng, t0_s=float(t0_s),
                max_delay_s=tau0 + 10.0 * float(model.delay_jitter_s),
            )
        elif realization.n != n or abs(realization.t0 - float(t0_s)) > 0.5 / fs or realization.fs != fs:
            raise ValueError("the realization is for a different record")
        phase = phase + realization.residual_phase(tau0)
        amp = realization.tx_amplitude(tau0)

    p_tx = float(transmit_power_w) * 10.0 ** (float(tx_chain_efficiency_db) / 10.0)
    x = math.sqrt(p_tx) * amp * np.conj(h) * np.exp(1j * phase)
    return np.where(live, x, 0.0 + 0.0j)


def leakage_over_echo_db(
    instrument: Any,
    sweep: SweepSchedule,
    model: LeakageModel,
    *,
    f_hz: np.ndarray | None = None,
    r_max_m: float | None = None,
) -> dict[str, np.ndarray]:
    """Leakage over the reference echoes at each step, with this model's profile.

    The echoes are those of :func:`tss.l2_instrument.reference_echo.reference_levels`
    (the free-space echo at R_max and the weakest echo INS-10 admits); the
    isolation is :func:`isolation_profile_db` at each frequency rather than one
    number.  With a flat model this reproduces the Phase 0 reference exactly.
    """
    from tss.l2_instrument import reference_echo

    lev = reference_echo.reference_levels(
        instrument, sweep, f_hz=f_hz, r_max_m=r_max_m, isolation_db=0.0
    )
    iso = isolation_profile_db(lev["f_hz"], model)
    return {
        "f_hz": lev["f_hz"],
        "isolation_db": iso,
        "leakage_dbw": lev["leakage_dbw"] - iso,
        "leakage_over_rmax_echo_db": lev["leakage_over_rmax_echo_db"] - iso,
        "leakage_over_weakest_db": lev["leakage_over_weakest_db"] - iso,
        "leakage_over_noise_db": lev["leakage_over_noise_db"] - iso,
        "beat_min_hz": lev["beat_min_hz"],
        "beat_max_hz": lev["beat_max_hz"],
        "step_rate_hz_per_s": lev["step_rate_hz_per_s"],
    }


# ----------------------------------------------------------------------------
# The analytic forms the model is tested against (Rev G section 5.9)
# ----------------------------------------------------------------------------


def range_correlation_suppression_db(f_m_hz: np.ndarray | float, tau_s: np.ndarray | float) -> np.ndarray:
    """``10 log10(4 sin^2(pi f_m tau))``: the suppression of a transmit spur."""
    x = np.pi * np.asarray(f_m_hz, dtype=float) * np.asarray(tau_s, dtype=float)
    return 10.0 * np.log10(4.0 * np.sin(x) ** 2)


def fm_spur_dbc(deviation_hz: float, rate_hz: float) -> float:
    """First FM sideband of the transmitted chirp relative to its carrier, exact."""
    beta = float(deviation_hz) / float(rate_hz)
    return float(20.0 * np.log10(abs(_special.j1(beta) / _special.j0(beta))))


def dechirped_fm_sideband_dbc(deviation_hz: float, rate_hz: float, tau_s: float) -> float:
    """That sideband after de-chirping against a common reference, exact.

    The residual phase deviation is ``(2 deviation / rate) sin(pi rate tau)``.
    """
    beta = 2.0 * float(deviation_hz) / float(rate_hz) * math.sin(math.pi * float(rate_hz) * float(tau_s))
    return float(20.0 * np.log10(abs(_special.j1(beta) / _special.j0(beta))))


def small_angle_sideband_dbc(deviation_hz: float, tau_s: float) -> float:
    """Rev G's form for ``f_m tau << 1``: ``(pi deviation tau)^2``."""
    return float(20.0 * math.log10(math.pi * float(deviation_hz) * float(tau_s)))


def ripple_sideband_dbc(ripple_pp_db: float) -> float:
    """First AM sideband of a sinusoidal isolation ripple, relative to the carrier.

    The amplitude is ``exp(-z sin theta)`` with ``z = (ln 10 / 20) ripple / 2``,
    whose Fourier series has ``I_n(z)`` coefficients, so the first sideband is
    ``I1(z) / I0(z)`` exactly.
    """
    z = _LN10_OVER_20 * 0.5 * float(ripple_pp_db)
    return float(20.0 * np.log10(_special.i1(z) / _special.i0(z)))


def ripple_modulation_hz(rate_hz_per_s: float, ripple_period_hz: float) -> float:
    """The AM rate a frequency ripple becomes under a sweep: ``K / period``."""
    return abs(float(rate_hz_per_s)) / float(ripple_period_hz)


def jitter_phase_rms_rad(f_hz: float, delay_jitter_s: float) -> float:
    """RMS de-chirped phase from delay jitter: ``2 pi f sigma_tau``."""
    return 2.0 * math.pi * float(f_hz) * float(delay_jitter_s)


# ----------------------------------------------------------------------------
# The moment-method anchor (WP3 decision G1(b))
# ----------------------------------------------------------------------------

#: Where the anchoring runs are kept, relative to the repository root.  Produced
#: by ``scripts/matlab/wp3_isolation_mom.m`` (free space) and
#: ``scripts/matlab/wp3_isolation_mom_bus.m`` (on the 6U bus), MATLAB Antenna
#: Toolbox; see the README beside them.
MOM_ANCHOR_RELPATH = "data/reference/wp3_isolation_mom/isolation_mom.csv"
MOM_BUS_RELPATH = "data/reference/wp3_isolation_mom/isolation_mom_bus.csv"

#: (geometry, match) -> (file, column).  The symmetric bus mounts
#: (crossed_ends_centred, crossed_colocated) sit at the numerical floor of the
#: solver and are not offered as profiles.
_MOM_COLUMNS = {
    ("bus_crossed_ends_offset", "simultaneous"): ("bus", "iso_crossed_ends_offset_db"),
    ("bus_crossed_ends_offset", "self"): ("bus", "iso_crossed_ends_offset_self_db"),
    ("bus_parallel_ends", "simultaneous"): ("bus", "iso_parallel_ends_db"),
    ("bus_parallel_ends", "self"): ("bus", "iso_parallel_ends_self_db"),
    ("crossed_offset", "simultaneous"): ("free", "iso_crossed_offset_db"),
    ("crossed_offset", "self"): ("free", "iso_crossed_offset_self_db"),
    ("parallel", "simultaneous"): ("free", "iso_parallel_db"),
    ("parallel", "self"): ("free", "iso_parallel_self_db"),
}


def mom_anchor_table(which: str = "free") -> dict[str, np.ndarray]:
    """An anchoring run's columns, keyed by the CSV header.

    ``which`` is ``"free"`` (two dipoles in free space) or ``"bus"`` (on the 6U
    bus, with the coupling coefficient of each mount).
    """
    import csv

    from tss.config import repo_root

    rel = {"free": MOM_ANCHOR_RELPATH, "bus": MOM_BUS_RELPATH}[str(which)]
    path = repo_root() / rel
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{path} is empty")
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def mom_anchor_model(
    geometry: str = "bus_crossed_ends_offset", match: str = "simultaneous", **kwargs: Any
) -> LeakageModel:
    """A tabulated :class:`LeakageModel` on a moment-method isolation.

    ``geometry``:

    * ``"bus_crossed_ends_offset"`` (the default anchor): on the 6U bus,
      transmit dipole along x 5 cm below it, receive along y 5 cm above it, the
      receive feed displaced 3 cm and 8 cm from the axis;
    * ``"bus_parallel_ends"``: both along x, the worst case;
    * ``"crossed_offset"`` and ``"parallel"``: the free-space runs.

    ``match`` is ``"simultaneous"`` (both ports conjugate-matched with the other
    antenna present: the least isolation the pair allows) or ``"self"``.  All
    include the matching-coil loss at Q = 150.
    """
    try:
        which, column = _MOM_COLUMNS[(str(geometry), str(match))]
    except KeyError:
        raise ValueError(f"no anchor column for geometry={geometry!r}, match={match!r}") from None
    tab = mom_anchor_table(which)
    kwargs.setdefault("label", f"MoM anchor: {geometry}, {match} match")
    return LeakageModel.from_table(tab["f_hz"], tab[column], **kwargs)


def coupled_resonator_isolation_db(k: np.ndarray | float, q: float) -> np.ndarray:
    """Isolation of two identical resonated antennas coupled by ``k``, exact.

    With each antenna's coil loss ``R = |X| / Q`` dominating and the mutual
    impedance ``j k |X|``, both ports conjugate-matched together give
    ``1 / MAG = kappa + sqrt(kappa^2 - 1)``, ``kappa = 1 + 2 / (k Q)^2``, which is
    ``-20 log10(k Q / 2)`` for ``k Q << 1``.  This is why the isolation is flat
    at the bottom of the band, where the coil loss dominates, rather than
    growing as the efficiency falls.
    """
    kq = np.asarray(k, dtype=float) * float(q)
    kappa = 1.0 + 2.0 / kq**2
    return 10.0 * np.log10(kappa + np.sqrt(kappa**2 - 1.0))

# ----------------------------------------------------------------------------
# Phase 4: the residual and its score
# ----------------------------------------------------------------------------

#: INS-24: the residual shall lie this far below the weakest admitted echo.
INS24_MARGIN_DB = 10.0
#: INS-23 threshold, dB.
INS23_THRESHOLD_DB = 40.0


def _step_order(sweep: SweepSchedule, n_steps: int) -> np.ndarray:
    """Ascending-frequency index of each step in transmission order."""
    idx = np.arange(n_steps)
    return idx if int(sweep.direction) >= 0 else idx[::-1]


def dechirped_residual(
    residual: np.ndarray,
    sweep: SweepSchedule,
    fs_hz: float,
    *,
    levels: Mapping[str, np.ndarray],
    t0_s: float = 0.0,
    sub_bands: int = 4,
    min_part_bins: int = 64,
    false_alarm: float = 1.0e-3,
    confirm_false_alarm: float = 1.0e-2,
) -> dict[str, np.ndarray]:
    """The residual leakage in each step's echo band, as an L3 cell would see it.

    ``residual`` is the leakage's contribution to the de-chirped record after
    the analog chain, in root-watts referred to the receiver input (the L2
    baseline passes the leakage alone through the front end and the
    high-pass).  ``levels`` is
    :func:`tss.l2_instrument.reference_echo.reference_levels` for the same
    instrument and sweep, in ascending frequency.

    Each step is transformed over its whole dwell with the Hann window L3 uses.
    The level in bin ``k`` is ``|X_k|^2 / (sum w)^2``: the power of a tone that
    would give that bin, so it compares directly with an echo; for a noise-like
    residual it is the density times the cell bandwidth, which is how
    ``noise_cell_dbw`` is defined.  The echo band is ``2 K_i R / c`` from 100 km
    to R_max, on the side the echoes fall (positive for an up-sweep).  In it:

    * ``pedestal_dbw``: the largest of the means over up to ``sub_bands``
      equal parts of the band, each at least ``min_part_bins`` wide, so a
      pedestal that is not flat is read where it is highest (the largest of
      four means of 64 exponential bins reads a flat pedestal about 1 dB high,
      of 240 bins about 0.4 dB);
      ``pedestal_mean_dbw`` is the mean over the whole band;
    * ``line_dbw``: the strongest line, NaN when there is none (below);
    * ``level_dbw``: the larger of the two, which is what INS-24 is scored on.

    A line is a bin that stands above its part's pedestal (the part's mean
    without its largest bin) by more than ``ln(N / false_alarm)``, which a
    noise-like pedestal of ``N`` exponentially distributed bins exceeds in a
    step with probability ``false_alarm``, and that is confirmed by the same
    frequency standing ``ln(1 / confirm_false_alarm)`` above the pedestal in
    a neighbouring step.  A transmitter line sits at the same de-chirped
    offset in every step; a noise peak does not.  Without the confirmation one
    step in a sweep or so would carry a false line some 10 dB above its
    pedestal, and the worst step is what INS-24 is scored on.  A line is found
    only if it stands 11 to 12 dB above the pedestal (the 100 to 1000 bins of
    an A1 echo band); deterministic spurs are checked against INS-22
    analytically (Phase 2) and need not be found here.

    Returned arrays are in ascending frequency.
    """
    x = np.asarray(residual, dtype=complex)
    fs = float(fs_hz)
    dwell = float(sweep.step_dwell_s)
    n_steps = sweep_law.step_count(sweep)
    order = _step_order(sweep, n_steps)
    d = 1.0 if int(sweep.direction) >= 0 else -1.0
    keys = ("pedestal_dbw", "pedestal_mean_dbw", "line_dbw", "line_offset_hz", "level_dbw")
    out = {k: np.full(n_steps, np.nan) for k in keys}
    bands: list[tuple[np.ndarray, np.ndarray, np.ndarray] | None] = [None] * n_steps
    for i in range(n_steps):                       # transmission order
        lo = int(math.ceil((i * dwell - float(t0_s)) * fs))
        hi = int(math.floor(((i + 1) * dwell - float(t0_s)) * fs))
        lo, hi = max(lo, 0), min(hi, x.size)
        if hi - lo < 16:
            continue
        seg = x[lo:hi]
        w = np.hanning(seg.size)
        spec = np.abs(np.fft.fft(seg * w)) ** 2 / np.sum(w) ** 2
        f = np.fft.fftfreq(seg.size, 1.0 / fs)
        j = order[i]
        b_lo, b_hi = float(levels["beat_min_hz"][j]), float(levels["beat_max_hz"][j])
        sel = np.flatnonzero((d * f >= b_lo) & (d * f <= min(b_hi, fs / 2.0)))
        if sel.size == 0:
            continue
        sel = sel[np.argsort(d * f[sel])]
        vals = spec[sel]
        n_parts = max(1, min(int(sub_bands), vals.size // max(int(min_part_bins), 1)))
        parts = np.array_split(np.arange(vals.size), n_parts)
        ref = np.empty(vals.size)
        for q in parts:
            v = vals[q]
            ref[q] = (np.sum(v) - np.max(v)) / max(v.size - 1, 1)
        bands[i] = (f[sel], vals, vals / np.maximum(ref, 1e-300))
        ped = max(float(np.mean(vals[q])) for q in parts)
        out["pedestal_dbw"][j] = 10.0 * np.log10(ped + 1e-300)
        out["pedestal_mean_dbw"][j] = 10.0 * np.log10(float(np.mean(vals)) + 1e-300)
        out["level_dbw"][j] = out["pedestal_dbw"][j]

    t_confirm = math.log(1.0 / float(confirm_false_alarm))
    for i, band in enumerate(bands):
        if band is None:
            continue
        fb, vals, excess = band
        thresh = math.log(vals.size / float(false_alarm))
        best = -np.inf
        best_f = np.nan
        for m in np.flatnonzero(excess > thresh):
            confirmed = False
            for nb in (i - 1, i + 1):
                if 0 <= nb < n_steps and bands[nb] is not None:
                    fn, _vn, en = bands[nb]
                    k = int(np.argmin(np.abs(fn - fb[m])))
                    if abs(fn[k] - fb[m]) <= fs / max(vals.size, 1) and en[max(k - 1, 0) : k + 2].max() > t_confirm:
                        confirmed = True
            if confirmed and vals[m] > best:
                best, best_f = float(vals[m]), float(fb[m])
        if np.isfinite(best):
            j = order[i]
            out["line_dbw"][j] = 10.0 * np.log10(best + 1e-300)
            out["line_offset_hz"][j] = best_f
            out["level_dbw"][j] = max(out["level_dbw"][j], out["line_dbw"][j])
    out["f_hz"] = np.asarray(levels["f_hz"], dtype=float)
    return out


def score_ins23_ins24(
    f_hz: np.ndarray,
    model: LeakageModel,
    residual: Mapping[str, np.ndarray] | None = None,
    *,
    levels: Mapping[str, np.ndarray] | None = None,
    canceller_residual: Mapping[str, np.ndarray] | None = None,
) -> dict[str, Any]:
    """Score the isolation and the residual together.

    INS-23 is the least isolation of ``model`` over ``f_hz``.  INS-24 (decision
    G5a, 2026-10-28) is the margin of the residual's level (the larger of its
    pedestal and its strongest line, :func:`dechirped_residual`) below the
    weakest echo INS-10 admits, ``weakest_echo_dbw - level_dbw``, at its worst
    step;
    10 dB is the threshold.  The margin against the free-space echo from R_max
    (Rev G's reading, which INS-22 uses) is reported beside it, and so is the
    margin with an ideal canceller ahead of the receiver when
    ``canceller_residual`` is given.  Steps whose residual could not be
    measured are left out and counted.
    """
    f = np.asarray(f_hz, dtype=float)
    iso = isolation_profile_db(f, model)
    i_iso = int(np.argmin(iso))
    out: dict[str, Any] = {
        "isolation_db_min": float(iso[i_iso]),
        "isolation_db_min_freq_hz": float(f[i_iso]),
        "ins23_met": bool(iso[i_iso] >= INS23_THRESHOLD_DB),
    }
    if residual is None or levels is None:
        return out

    def _margins(res: Mapping[str, np.ndarray]) -> dict[str, Any]:
        pk = np.asarray(res["level_dbw"], dtype=float)
        weak = np.asarray(levels["weakest_echo_dbw"], dtype=float) - pk
        rmax = np.asarray(levels["echo_rmax_dbw"], dtype=float) - pk
        ok = np.isfinite(weak)
        if not np.any(ok):
            return {"margin_db": None}
        j = int(np.flatnonzero(ok)[np.argmin(weak[ok])])
        return {
            "margin_db": float(weak[j]),
            "worst_freq_hz": float(np.asarray(levels["f_hz"])[j]),
            "margin_db_rmax": float(np.nanmin(rmax)),
            "steps_failing": int(np.count_nonzero(weak[ok] < INS24_MARGIN_DB)),
            "steps_scored": int(np.count_nonzero(ok)),
            "per_step_margin_db": weak,
            "per_step_margin_db_rmax": rmax,
        }

    m = _margins(residual)
    out.update({
        "ins24_margin_db": m["margin_db"],
        "ins24_worst_freq_hz": m.get("worst_freq_hz"),
        "ins24_margin_db_rmax": m.get("margin_db_rmax"),
        "ins24_steps_failing": m.get("steps_failing"),
        "ins24_steps_scored": m.get("steps_scored"),
        "ins24_per_step_margin_db": m.get("per_step_margin_db"),
        "ins24_per_step_margin_db_rmax": m.get("per_step_margin_db_rmax"),
    })
    if canceller_residual is not None:
        c = _margins(canceller_residual)
        out["ins24_margin_db_canceller"] = c["margin_db"]
        out["ins24_steps_failing_canceller"] = c.get("steps_failing")
    return out
