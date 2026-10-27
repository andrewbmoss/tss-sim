"""Layer: L2
WP: WP3
Claims: C1, C4
Status: partial
Requirements: INS-11, INS-12, INS-14, INS-21, INS-26

The baseline L2 tier: echoes in, complex baseband samples out.

This is the tier registered as ``l2/baseline`` and named by
``config/default.yaml``.  It composes the other modules in this package rather
than computing anything itself, which is the point: the antenna model, the noise
reference, the sweep law, the receive chain and the oscillator are each testable
on their own, and this module is the wiring.

WHAT IT PRODUCES
----------------
De-chirped complex baseband, one record per sweep.  The instrument de-chirps in
hardware against the same reference that generated the transmission, which is
INS-21, so the samples that leave the receiver are already at baseband and an
echo appears as a tone rather than as a delayed sweep.  Generating the passband
signal and then mixing it down in software would be a more literal simulation
and would cost a factor of several hundred in sample rate for no change in the
answer, so it is not done.

The de-chirped phase of an echo at two-way group delay ``tau`` is exact rather
than approximated:

    phi_b(t) = 2 pi [ f0 tau + d K tau t - d K tau^2 / 2 ]

from subtracting the transmit phase law of :mod:`tss.l2_instrument.waveform` at
``t - tau`` from the same law at ``t``.  The middle term is the beat frequency,
``d K tau``, which is negative for a down-sweep; carrying the sign rather than
the magnitude is what lets an up/down pair be processed as a pair, which C3
needs.  The last term is a constant phase and is the piece a first-order
treatment drops.

DISPERSION WITHIN A DWELL
-------------------------
The delay of an ionospheric echo depends on frequency, and the sweep moves
through frequency while a step is being collected.  The delay therefore changes
during the dwell: over a 5 per cent step in an exponential topside of scale
height ``H`` the virtual range changes by about ``0.1 H / mu_sat``, which at
3 MHz is ten to forty of the 1 km range cells, and near foF2 it is far more.  A
tier that gives each step one delay produces a tone where the instrument would
see a slide, and every resolution and sensitivity figure taken from it is the
undispersed case.

With ``in_dwell_dispersion`` on, which is the default, the rows of each
magnetoionic mode are joined into a trace ``tau_g(f)``, piecewise linear through
the rows and extended linearly past its ends, and the echo is built from that
trace rather than from one delay per row.  The received sample at time ``t``
left the antenna at ``t_e``, where

    t = t_e + tau_g(f(t_e))

and its phase is that of the transmitted law at ``t_e`` advanced along the phase
path, so the de-chirped phase is

    phi_b(t) = phi(t) - phi(t_e) - 2 pi f_e tau_g(f_e) + 2 pi int^f_e tau_g df

with ``f_e = f(t_e)``.  This is the stationary-phase result for a signal whose
frequency changes slowly against the dispersion.  Its derivative is
``2 pi (f(t) - f_e)``, the beat against the frequency that is arriving, and for
a constant delay the last two terms cancel and it reduces to the expression
above.  ``t_e`` is found by fixed-point iteration, which converges quickly
because ``d tau / d t`` is of order one per cent.

A mode whose rows are not single valued in frequency, as a ray-traced tier with
several rays per frequency would give, has no unique trace.  Its rows keep one
delay each and the count is reported in ``flags["rows_without_trace"]``.

WHEN AN ECHO IS PRESENT
-----------------------
An echo at sounding frequency ``f_k`` is only there while the sweep is near
``f_k``.  The tier gives each row a rectangular window centred on the time the
sweep crosses that row's frequency, of width equal to the frequency spacing of
the echo set divided by the sweep rate.  That is the dwell the sounding actually
gets at that frequency, and it is why a fine frequency step costs dwell.

The window is rectangular and not tapered.  A taper would put a more realistic
skirt on each tone, and it would also be a processing choice made at L2 where it
belongs at L3, so the taper is left for the transform to apply.

ABSORPTION IS NOT APPLIED HERE
------------------------------
``EchoSet.amplitude_db`` from the L1 analytic tier already includes the
absorption integral; ``EchoSet.absorption_db`` is carried beside it for
traceability, not as a second loss to apply.  This tier therefore uses
``amplitude_db`` and records ``absorption_db`` in the flags.  Applying both
would double-count the absorption, which at the top of the band near foF2 is
several decibels.

THE SAMPLE RATE
---------------
The digitiser rate is a property of the instrument, not of the echoes.  It is
``sample_rate_hz`` when the configuration fixes one, and otherwise the INS-06
rate ``4 K_f R / c`` at the sweep's fastest rate, with ``R`` the design range
(INS-04's 1050 km by default), so that the unambiguous range ``c fs / (4 K_f)``
equals the design range.  Until 2026-10-17 the rate was 2.5 times the largest
beat among the simulated echoes, which made the INS-04 score a statement about
the simulated ionosphere (always 1.25 times its deepest echo) and sampled a
record with no usable echo at 2.5 Hz.  Echoes beyond the unambiguous range fold
back into the record, as they would in the instrument, because the analytic
phase is evaluated at the sample times; they are counted in
``flags["rows_beyond_unambiguous_range"]``.

WHAT IS NOT IN IT YET
---------------------
Transmit-to-receive leakage.  :mod:`tss.l2_instrument.leakage` models it (WP3
Phase 1) but this tier does not yet add it (Phase 3), so the samples this tier
produces have no leakage in them and the receiver therefore never clips on leakage.  That is the single
largest gap between this tier and the instrument, it is recorded in every
``ReceivedSignal`` this tier returns under the flag ``leakage_modelled``, and no
verdict on INS-11 or INS-24 may be taken from a run made without it.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace as _replace
from typing import Any

import numpy as np

from tss import sweep_law
from tss.constants import C_LIGHT, K_BOLTZ
from tss.interfaces import EchoSet, ReceivedSignal, SweepSchedule
from tss.l2_instrument.antenna import ShortDipole, antenna_noise_temperature_k
from tss.l2_instrument.noise import (
    is_extrapolated,
    noise_density_dbm_per_hz,
    sky_temperature,
)
from tss.l2_instrument.oscillator import Oscillator
from tss.l2_instrument.receiver import Receiver
from tss.l2_instrument.waveform import (
    sweep_bandwidth_hz,
    sweep_duration_s,
    validate_schedule,
)
from tss.registry import register
from tss.units import db, undb, w_to_dbm

__all__ = [
    "BaselineInstrument",
    "DEFAULT_DESIGN_RANGE_M",
    "DEFAULT_POLARIZATION_LOSS_DB",
    "DelayTrace",
]

#: Design unambiguous range, m: the INS-04 threshold.  With no configured
#: sample rate the digitiser runs at the INS-06 rate that reaches it.
DEFAULT_DESIGN_RANGE_M = 1.05e6

#: Relative guard on the INS-06 rate, so that the unambiguous range computed
#: back from it is not a rounding error short of the design range.
_RATE_GUARD = 1.0e-9

#: Fixed-point passes used to find the emission time of a dispersed echo.  The
#: contraction factor is ``|d tau_g / d f| K``, about 0.005 for the A1 sweep in
#: a 150 km scale height, so six passes are far more than double precision needs.
_DISPERSION_ITERATIONS = 6


class DelayTrace:
    """Group delay against frequency along one echo trace, and its integral.

    Piecewise linear through the knots and extended linearly beyond them with
    the slope of the end segment, so that the delay and its frequency integral
    are exact for the model and continuous everywhere.  A trace needs at least
    two knots at distinct frequencies.
    """

    def __init__(self, freq_hz: np.ndarray, delay_s: np.ndarray) -> None:
        f = np.asarray(freq_hz, dtype=float)
        tau = np.asarray(delay_s, dtype=float)
        order = np.argsort(f)
        f, tau = f[order], tau[order]
        if f.size < 2 or np.any(np.diff(f) <= 0.0):
            raise ValueError("a delay trace needs two or more distinct frequencies")
        self.freq_hz = f
        self.tau_s = tau
        self.slope_s_per_hz = np.diff(tau) / np.diff(f)
        # Cycles of phase path at each knot, int tau df, anchored so that a
        # constant delay gives f * tau exactly.
        self.cycles = np.concatenate(
            [[f[0] * tau[0]], f[0] * tau[0] + np.cumsum(0.5 * (tau[:-1] + tau[1:]) * np.diff(f))]
        )

    def _segment(self, f_hz: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
        f = np.asarray(f_hz, dtype=float)
        j = np.clip(np.searchsorted(self.freq_hz, f, side="right") - 1, 0, self.freq_hz.size - 2)
        return j, f - self.freq_hz[j]

    def delay_s(self, f_hz: np.ndarray | float) -> np.ndarray:
        """Group delay at ``f_hz``."""
        j, d = self._segment(f_hz)
        return self.tau_s[j] + self.slope_s_per_hz[j] * d

    def phase_cycles(self, f_hz: np.ndarray | float) -> np.ndarray:
        """``int tau_g df`` up to ``f_hz``, in cycles."""
        j, d = self._segment(f_hz)
        return self.cycles[j] + self.tau_s[j] * d + 0.5 * self.slope_s_per_hz[j] * d**2


def _traces_by_mode(
    freq: np.ndarray, delay: np.ndarray, mode: np.ndarray | None, usable: np.ndarray
) -> dict[str, DelayTrace]:
    """One trace per mode, for the modes whose usable rows define one."""
    if mode is None or np.shape(mode) != np.shape(freq):
        labels = np.array([""] * freq.size, dtype=object)
    else:
        labels = np.array([str(m) for m in np.asarray(mode, dtype=object).ravel()], dtype=object)
    out: dict[str, DelayTrace] = {}
    for label in sorted(set(labels[usable].tolist())):
        sel = usable & (labels == label)
        f = freq[sel]
        if f.size < 2 or np.unique(f).size != f.size:
            continue
        out[label] = DelayTrace(f, delay[sel])
    return out

#: Polarization loss, dB.  -3 dB for the O and X split and -3 dB for an
#: elliptically polarized wave received on a linear antenna.  The same figure
#: the link budget uses as its default.  It is not in the instrument
#: configuration because it is a property of the propagation and the antenna
#: together rather than of the instrument alone.
DEFAULT_POLARIZATION_LOSS_DB = -6.0


@register("l2", "baseline")
class BaselineInstrument:
    """Compose the L2 modules into an :class:`~tss.interfaces.L2Instrument`.

    Parameters
    ----------
    tx_power_w
        Transmit power at the power amplifier output.  INS-14's threshold is 1 W.
    tx_chain_efficiency_db
        Loss between the amplifier and the antenna terminals, negative, not
        counting the matching network.  The matching-network loss is part of
        the antenna's radiation efficiency and is carried there only (decision
        of 2026-10-17); the baseline configuration states 0 dB here.
    antenna
        A :class:`~tss.l2_instrument.antenna.ShortDipole`, or ``None`` to build
        the default 5 m one.  The same antenna is used for transmit and receive;
        the instrument shares one aperture and a separate receive antenna would
        be a different configuration.
    receiver, oscillator
        The receive chain and the frequency reference.
    stated_efficiency_db
        When given, this flat radiation efficiency replaces the frequency
        dependent one the physical antenna model computes.  The baseline
        configuration does not state one (decision of 2026-10-17): its
        efficiency follows the sounding frequency, matching-coil loss included.
    polarization_loss_db, extra_loss_db
        Applied to every echo.
    sample_rate_hz
        Complex digitiser rate, Hz.  ``None``, the default, runs at the INS-06
        rate for the sweep in use: see :meth:`sample_rate_for` and the module
        docstring, THE SAMPLE RATE.
    design_range_m
        The unambiguous range the INS-06 rate is set to reach.  Default
        :data:`DEFAULT_DESIGN_RANGE_M`, the INS-04 threshold.  Ignored when
        ``sample_rate_hz`` is given.
    max_samples
        Cap on the record length.  A sweep longer than this is generated only up
        to the cap and the truncation is recorded in the flags, together with
        how many echo rows fell outside the generated window.  A silently
        shortened record would read as an ionogram with a missing top end.
    sky_model
        Passed to :func:`tss.l2_instrument.noise.sky_temperature`.
    quantise
        Whether to run the quantiser.  Off is useful for isolating whether an
        effect is quantisation or something else.
    in_dwell_dispersion
        Whether the delay of an echo follows its trace through the dwell, as it
        does in the ionosphere, or is held at one value per row.  On by default.
        Off reproduces the undispersed samples this tier produced before
        2026-10-05 and is there so the two can be compared.
    """

    def __init__(
        self,
        *,
        tx_power_w: float = 1.0,
        tx_chain_efficiency_db: float = 0.0,
        antenna: ShortDipole | None = None,
        receiver: Receiver | None = None,
        oscillator: Oscillator | None = None,
        stated_efficiency_db: float | None = None,
        polarization_loss_db: float = DEFAULT_POLARIZATION_LOSS_DB,
        extra_loss_db: float = 0.0,
        sample_rate_hz: float | None = None,
        design_range_m: float = DEFAULT_DESIGN_RANGE_M,
        max_samples: int = 1 << 20,
        sky_model: str = "franklin",
        quantise: bool = True,
        in_dwell_dispersion: bool = True,
        isolation: Mapping[str, Any] | None = None,
        resources: Mapping[str, Any] | None = None,
        configuration_name: str = "",
        configuration_notes: tuple[str, ...] = (),
    ) -> None:
        if tx_power_w <= 0.0:
            raise ValueError("tx_power_w must be positive")
        if sample_rate_hz is not None and not float(sample_rate_hz) > 0.0:
            raise ValueError("sample_rate_hz must be positive")
        if not float(design_range_m) > 0.0:
            raise ValueError("design_range_m must be positive")
        if max_samples < 16:
            raise ValueError("max_samples must be at least 16")
        self.tx_power_w = float(tx_power_w)
        self.tx_chain_efficiency_db = float(tx_chain_efficiency_db)
        self.antenna = antenna if antenna is not None else ShortDipole()
        self.receiver = receiver if receiver is not None else Receiver()
        self.oscillator = oscillator if oscillator is not None else Oscillator()
        self.stated_efficiency_db = (
            None if stated_efficiency_db is None else float(stated_efficiency_db)
        )
        self.polarization_loss_db = float(polarization_loss_db)
        self.extra_loss_db = float(extra_loss_db)
        self.sample_rate_hz = None if sample_rate_hz is None else float(sample_rate_hz)
        self.design_range_m = float(design_range_m)
        self.max_samples = int(max_samples)
        self.sky_model = str(sky_model)
        self.quantise = bool(quantise)
        self.in_dwell_dispersion = bool(in_dwell_dispersion)
        self.isolation = dict(isolation or {})
        self.resources = dict(resources or {})
        self.configuration_name = str(configuration_name)
        self.configuration_notes = tuple(str(n) for n in configuration_notes)

    # -- configuration -------------------------------------------------------

    @classmethod
    def from_config(cls, cfg: Mapping[str, Any], **overrides: Any) -> BaselineInstrument:
        """Build from the instrument configuration shape.

        The shape is the one in ``config/instrument/baseline_6u.yaml``: blocks
        named ``transmit``, ``antenna``, ``receiver``, ``isolation``,
        ``oscillator`` and ``resources``.  Keys that are not recognised are
        ignored rather than rejected, so that adding a key to the configuration
        for another consumer does not break this one.
        """
        transmit = dict(cfg.get("transmit", {}))
        antenna_cfg = dict(cfg.get("antenna", {}))
        antenna = ShortDipole.from_config(antenna_cfg)
        receiver_cfg = dict(cfg.get("receiver", {}))
        rate = receiver_cfg.get("sample_rate_hz")
        kwargs: dict[str, Any] = dict(
            tx_power_w=float(transmit.get("power_w", 1.0)),
            tx_chain_efficiency_db=float(transmit.get("efficiency_db", 0.0)),
            antenna=antenna,
            receiver=Receiver.from_config(dict(cfg.get("receiver", {}))),
            oscillator=Oscillator.from_config(dict(cfg.get("oscillator", {}))),
            stated_efficiency_db=antenna_cfg.get("radiation_efficiency_db"),
            in_dwell_dispersion=bool(cfg.get("in_dwell_dispersion", True)),
            isolation=dict(cfg.get("isolation", {})),
            resources=dict(cfg.get("resources", {})),
            configuration_name=str(cfg.get("name", "")),
            configuration_notes=tuple(cfg.get("configuration_overrides", ()) or ()),
            sample_rate_hz=None if rate is None else float(rate),
            design_range_m=float(receiver_cfg.get("design_range_m", DEFAULT_DESIGN_RANGE_M)),
        )
        kwargs.update(overrides)
        return cls(**kwargs)

    def sample_rate_for(self, sweep: SweepSchedule) -> float:
        """Complex digitiser rate for a sweep, Hz.

        The configured ``sample_rate_hz`` when there is one.  Otherwise the
        INS-06 rate ``4 K_f R / c`` at the sweep's fastest rate (for a stepped
        sweep ``k_f_hz_per_s`` is the peak step rate), with ``R`` the design
        range, so that ``c fs / (4 K_f)`` is the design range.  It never
        depends on the echoes.
        """
        if self.sample_rate_hz is not None:
            return self.sample_rate_hz
        k = abs(float(sweep.k_f_hz_per_s))
        if k <= 0.0:
            raise ValueError("SweepSchedule.k_f_hz_per_s must be non-zero")
        return 4.0 * k * self.design_range_m / C_LIGHT * (1.0 + _RATE_GUARD)

    # -- antenna terms -------------------------------------------------------

    def efficiency_db(self, f_mhz: np.ndarray | float) -> np.ndarray:
        """Radiation efficiency used for this run, dB.

        The stated figure when the configuration gives one, otherwise the
        physical model.  Which was used is in :meth:`describe`.
        """
        if self.stated_efficiency_db is not None:
            return np.full(np.shape(np.asarray(f_mhz, dtype=float)), self.stated_efficiency_db)
        return np.asarray(self.antenna.efficiency_db(f_mhz), dtype=float)

    def gain_dbi(self, f_mhz: np.ndarray | float) -> np.ndarray:
        """``G = D + eta``.  Directivity and efficiency added, never lumped."""
        return self.efficiency_db(f_mhz) + self.antenna.directivity_dbi

    def antenna_temperature_k(self, f_mhz: np.ndarray | float) -> np.ndarray:
        """``T_a = eta T_sky + (1-eta) T_phys`` at the given frequencies."""
        eta = undb(self.efficiency_db(f_mhz))
        t_sky = np.asarray(sky_temperature(f_mhz, self.sky_model), dtype=float)
        return antenna_noise_temperature_k(eta, t_sky, self.antenna.physical_temperature)

    # -- the contract --------------------------------------------------------

    def receive(
        self,
        echoes: EchoSet,
        sweep: SweepSchedule,
        rng: np.random.Generator,
    ) -> ReceivedSignal:
        """De-chirped baseband for one sweep.  Inputs are never modified."""
        validate_schedule(sweep)
        freq = np.asarray(echoes.freq_hz, dtype=float)
        delay = np.asarray(echoes.group_delay_s, dtype=float)
        amp_db = np.asarray(echoes.amplitude_db, dtype=float)
        if not (freq.shape == delay.shape == amp_db.shape):
            raise ValueError(
                "EchoSet.freq_hz, .group_delay_s and .amplitude_db must be parallel; "
                f"got {freq.shape}, {delay.shape}, {amp_db.shape}"
            )

        flags_in = (
            np.asarray(echoes.flags, dtype=object)
            if echoes.flags is not None
            else np.array([""] * freq.size, dtype=object)
        )
        usable = (
            np.isfinite(freq)
            & np.isfinite(delay)
            & np.isfinite(amp_db)
            & (delay > 0.0)
            & np.array(
                [str(f) not in ("evanescent", "penetrating") for f in flags_in.ravel()[: freq.size]],
                dtype=bool,
            )
        )

        duration_s = sweep_duration_s(sweep)
        direction = int(sweep.direction)
        k = float(sweep.k_f_hz_per_s)
        f0 = float(sweep.f_start_hz) if direction > 0 else float(sweep.f_stop_hz)

        max_delay = float(np.max(delay[usable])) if np.any(usable) else 0.0
        beat_max = k * max_delay
        # The digitiser rate is the instrument's (INS-06), never the echoes'.
        fs_hz = self.sample_rate_for(sweep)
        unambiguous_m = C_LIGHT * fs_hz / (4.0 * k)
        n_full = int(np.floor(duration_s * fs_hz))
        n_samples = int(min(max(n_full, 16), self.max_samples))
        truncated = n_samples < n_full
        modelled_s = n_samples / fs_hz

        t = np.arange(n_samples, dtype=float) / fs_hz
        stepped = sweep_law.is_stepped(sweep)
        if stepped:
            # The A1 stepped law: each step has its own rate, so the beat of an
            # echo depends on which step it was transmitted in and which step
            # it is received in.  The ideal law is used, as it is for the
            # uniform sweep below.
            ideal = _replace(sweep, nonlinearity_ppm=0.0)
            f_tx = sweep_law.frequency_hz(ideal, t)
            phi_rx = sweep_law.phase_rad(ideal, t)

            def tx_frequency(x: np.ndarray) -> np.ndarray:
                return sweep_law.frequency_hz(ideal, x)

            def tx_phase(x: np.ndarray) -> np.ndarray:
                return sweep_law.phase_rad(ideal, x)
        else:
            f_tx = f0 + direction * k * t

            def tx_frequency(x: np.ndarray) -> np.ndarray:
                return f0 + direction * k * x

            def tx_phase(x: np.ndarray) -> np.ndarray:
                return 2.0 * np.pi * (f0 * x + direction * k * x**2 / 2.0)

        traces = (
            _traces_by_mode(freq, delay, getattr(echoes, "mode", None), usable)
            if self.in_dwell_dispersion
            else {}
        )
        mode_in = getattr(echoes, "mode", None)
        if mode_in is None or np.shape(mode_in) != freq.shape:
            row_mode = np.array([""] * freq.size, dtype=object)
        else:
            row_mode = np.array([str(m) for m in np.asarray(mode_in, dtype=object).ravel()], dtype=object)

        # Per-row visit time and dwell window.  The spacing is taken from the
        # echo set itself, because that is the frequency step the sounding
        # actually used, and the dwell that follows from it is the dwell the
        # sounding actually got.
        unique_f = np.unique(freq[usable]) if np.any(usable) else np.array([], dtype=float)
        if stepped:
            # A stepped sweep fixes the dwell: each row is filed in the step
            # that contains its frequency and fills that step.
            dwell_s = float(sweep.step_dwell_s)
        elif unique_f.size > 1:
            df = float(np.median(np.diff(unique_f)))
            dwell_s = max(abs(df) / k, 2.0 / fs_hz)
        else:
            df = sweep_bandwidth_hz(sweep)
            dwell_s = max(abs(df) / k, 2.0 / fs_hz)

        samples = np.zeros(n_samples, dtype=complex)
        rows_used = 0
        rows_outside = 0
        rows_dispersed = 0
        max_delay_change_s = 0.0
        for row in np.nonzero(usable)[0]:
            f_k = float(freq[row])
            if stepped:
                t_k = float(sweep_law.time_at_frequency_s(sweep, f_k))
            else:
                t_k = direction * (f_k - f0) / k
            if not (0.0 <= t_k <= modelled_s):
                rows_outside += 1
                continue
            if stepped:
                i_step = int(sweep_law.step_index(sweep, t_k))
                w_lo, w_hi = i_step * dwell_s, (i_step + 1) * dwell_s
            else:
                w_lo, w_hi = t_k - dwell_s / 2.0, t_k + dwell_s / 2.0
            lo = int(np.floor(max(w_lo, 0.0) * fs_hz))
            hi = int(np.ceil(min(w_hi, modelled_s) * fs_hz))
            if hi <= lo:
                rows_outside += 1
                continue

            tau = float(delay[row])
            f_mhz = f_k / 1.0e6
            # The same antenna transmits and receives, so its gain appears
            # twice.  amp_db already carries the spreading and the absorption.
            gain_dbi = float(self.gain_dbi(f_mhz))
            power_dbw = (
                db(self.tx_power_w)
                + self.tx_chain_efficiency_db
                + gain_dbi
                + float(amp_db[row])
                + self.polarization_loss_db
                + self.extra_loss_db
                + gain_dbi
            )
            amplitude = np.sqrt(undb(power_dbw))

            tt = t[lo:hi]
            trace = traces.get(str(row_mode[row]))
            if trace is not None:
                # The delay follows the trace through the dwell.  Find when each
                # received sample left the antenna, t = t_e + tau_g(f(t_e)), and
                # take the phase of the transmitted law there advanced along the
                # phase path (module docstring, DISPERSION WITHIN A DWELL).  The
                # step transition needs no separate case here either, because
                # f(t_e) is read from the law whichever step t_e falls in.
                t_e = tt - tau
                for _ in range(_DISPERSION_ITERATIONS):
                    t_e = tt - np.clip(trace.delay_s(tx_frequency(t_e)), 0.0, None)
                f_e = tx_frequency(t_e)
                tau_e = tt - t_e
                phase = (
                    tx_phase(tt)
                    - tx_phase(t_e)
                    - 2.0 * np.pi * f_e * tau_e
                    + 2.0 * np.pi * trace.phase_cycles(f_e)
                )
                live = t_e >= 0.0
                samples[lo:hi] += np.where(live, amplitude * np.exp(1j * phase), 0.0)
                rows_used += 1
                rows_dispersed += 1
                if np.any(live):
                    seen = tau_e[live]
                    max_delay_change_s = max(max_delay_change_s, float(seen.max() - seen.min()))
                continue
            if stepped:
                # The exact de-chirped phase, phi(t) - phi(t - tau), from the
                # transmitted law itself.  Inside a step it is the uniform form
                # below with that step's start frequency and rate.  In the first
                # tau of a step the echo left the antenna during the previous
                # step, and the beat ramps from K_(i-1) tau to K_i tau; the
                # difference of the exact phases carries that transition with no
                # separate case.  Nothing had been transmitted before the sweep
                # started, so samples earlier than tau carry no echo.
                phase = phi_rx[lo:hi] - sweep_law.phase_rad(ideal, tt - tau)
                live = tt >= tau
                samples[lo:hi] += np.where(live, amplitude * np.exp(1j * phase), 0.0)
                rows_used += 1
                continue
            else:
                phase = 2.0 * np.pi * (
                    f0 * tau + direction * k * tau * tt - direction * k * tau**2 / 2.0
                )
            samples[lo:hi] += amplitude * np.exp(1j * phase)
            rows_used += 1

        # Echoes beyond the unambiguous range fold back into the record, as they
        # would in the instrument.  They are counted against the rate of the
        # step that sounds them, so a reader can tell a folded trace from a
        # measured one.
        rows_beyond = 0
        if np.any(usable):
            f_u = freq[usable]
            if stepped:
                k_u = np.asarray(
                    sweep_law.rate_hz_per_s(sweep, sweep_law.time_at_frequency_s(sweep, f_u)),
                    dtype=float,
                )
            else:
                k_u = np.full(f_u.shape, k)
            rows_beyond = int(np.count_nonzero(k_u * delay[usable] > fs_hz / 2.0))

        # Noise.  The antenna temperature follows the sweep, because the sky is
        # three orders of magnitude colder at the top of the band than at the
        # bottom and a single band-averaged figure would be wrong at both ends.
        t_ant = self.antenna_temperature_k(np.clip(f_tx, 1.0, None) / 1.0e6)
        t_sys = t_ant + self.receiver.noise_temperature_k
        sigma = np.sqrt(K_BOLTZ * t_sys * fs_hz / 2.0)
        samples = samples + rng.normal(0.0, sigma) + 1j * rng.normal(0.0, sigma)

        # Receive chain.  Gain, then clipping against the instantaneous dynamic
        # range referred to the median noise floor, then quantisation.
        amplified = self.receiver.apply_gain(samples)
        noise_floor_w = float(np.median(K_BOLTZ * t_sys * fs_hz)) * float(
            undb(self.receiver.gain_db)
        )
        clip_level_w = self.receiver.clip_level_w(noise_floor_w)
        clipped, n_clipped = self.receiver.clip(amplified, clip_level_w)
        out = self.receiver.quantise(clipped, clip_level_w) if self.quantise else clipped

        channels = max(int(self.receiver.channels), 1)
        rows_out = [out]
        # The vector configuration of C5 is three channels.  Until
        # tss.l2_instrument.vector_array exists the extra channels carry the same
        # signal with an independent noise realisation, which is enough to
        # exercise the interface and is not enough to test anything about C5.
        # The flag ``vector_channels_are_placeholders`` says so on every record.
        for _ in range(channels - 1):
            extra_noise = rng.normal(0.0, sigma) + 1j * rng.normal(0.0, sigma)
            chain = self.receiver.apply_gain(samples + extra_noise)
            chain, _ = self.receiver.clip(chain, clip_level_w)
            rows_out.append(
                self.receiver.quantise(chain, clip_level_w) if self.quantise else chain
            )
        block = np.vstack(rows_out)

        floor_f_mhz = unique_f / 1.0e6 if unique_f.size else np.array([1.0])
        noise_floor = noise_density_dbm_per_hz(
            self.antenna_temperature_k(floor_f_mhz) + self.receiver.noise_temperature_k
        )

        flags: dict[str, Any] = {
            # The samples are de-chirped baseband (see the module docstring).
            # L3 reads this and must not mix them against the sweep again.
            "dechirped": True,
            "fs_hz": fs_hz,
            "n_samples": n_samples,
            "sweep_duration_s": duration_s,
            "sweep_law": "stepped" if stepped else "uniform",
            # For a stepped sweep the de-chirped phase is the exact difference
            # phi(t) - phi(t - tau), so the transition in the first tau of each
            # step is in the samples (Technical Volume, stepped-rate sweep).
            "step_transition_modelled": bool(stepped),
            # Whether the delay of each echo follows its trace through the dwell
            # (Technical Volume, dispersion within a dwell).  rows_without_trace
            # counts the rows that kept one delay because their mode had no
            # single-valued trace; max_in_dwell_delay_change_s is the largest
            # change of delay any one row saw inside its window.
            "in_dwell_dispersion_modelled": bool(self.in_dwell_dispersion and rows_dispersed > 0),
            "rows_dispersed": rows_dispersed,
            "rows_without_trace": rows_used - rows_dispersed,
            "max_in_dwell_delay_change_s": max_delay_change_s,
            "modelled_duration_s": modelled_s,
            "truncated": bool(truncated),
            "rows_in": int(freq.size),
            "rows_used": rows_used,
            "rows_unusable": int(freq.size - int(np.count_nonzero(usable))),
            "rows_outside_window": rows_outside,
            "dwell_s": dwell_s,
            "max_beat_frequency_hz": beat_max,
            "sample_rate_rule": (
                "configured sample_rate_hz" if self.sample_rate_hz is not None
                else f"INS-06: 4 K_f R / c at the peak rate, R = {self.design_range_m:g} m"
            ),
            "unambiguous_range_m": unambiguous_m,
            "rows_beyond_unambiguous_range": rows_beyond,
            "clip_level_dbm": float(w_to_dbm(clip_level_w)),
            "n_clipped": n_clipped,
            "quantised": bool(self.quantise and self.receiver.n_bits is not None),
            "n_bits": self.receiver.n_bits,
            "sky_model": self.sky_model,
            "noise_floor_f_mhz": floor_f_mhz.tolist(),
            "noise_reference_extrapolated": bool(np.any(is_extrapolated(floor_f_mhz))),
            "leakage_modelled": False,
            "leakage_note": (
                "tss.l2_instrument.leakage models the leakage (WP3 Phase 1) but it "
                "is not yet added to these samples (Phase 3), so the receiver never "
                "clips on leakage and no INS-11 or INS-24 verdict may be taken from "
                "this run."
            ),
            "vector_channels_are_placeholders": channels > 1,
            "absorption_db_carried_not_applied": True,
        }

        return ReceivedSignal(
            samples=block,
            fs_hz=fs_hz,
            sweep=sweep,
            noise_floor_dbm_hz=np.asarray(noise_floor, dtype=float),
            flags=flags,
            source=self.source_label(),
        )

    # -- manifest ------------------------------------------------------------

    def source_label(self) -> str:
        name = self.configuration_name or "baseline"
        return (
            f"l2/baseline({name}, Pt={self.tx_power_w:g} W, "
            f"L={self.antenna.length_m:g} m, NF={self.receiver.noise_figure_db:g} dB)"
        )

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable identity and settings, for the run manifest."""
        return {
            "layer": "l2",
            "name": getattr(self, "name", type(self).__name__),
            "class": type(self).__name__,
            "module": type(self).__module__,
            "configuration": self.configuration_name,
            "tx_power_w": self.tx_power_w,
            "tx_chain_efficiency_db": self.tx_chain_efficiency_db,
            "polarization_loss_db": self.polarization_loss_db,
            "extra_loss_db": self.extra_loss_db,
            "sample_rate_hz": self.sample_rate_hz,
            "design_range_m": self.design_range_m,
            "configuration_notes": list(self.configuration_notes),
            "max_samples": self.max_samples,
            "sky_model": self.sky_model,
            "quantise": self.quantise,
            "in_dwell_dispersion": self.in_dwell_dispersion,
            "efficiency_source": (
                "configuration radiation_efficiency_db"
                if self.stated_efficiency_db is not None
                else "ShortDipole physical model"
            ),
            "stated_efficiency_db": self.stated_efficiency_db,
            "antenna": dict(self.antenna.describe()),
            "receiver": dict(self.receiver.describe()),
            "oscillator": dict(self.oscillator.describe()),
            "isolation": dict(self.isolation),
            "resources": dict(self.resources),
            "leakage_modelled": False,
            "matching_network_loss": (
                "carried in the antenna radiation efficiency only (decision of "
                "2026-10-17); tx_chain_efficiency_db excludes it"
            ),
        }
