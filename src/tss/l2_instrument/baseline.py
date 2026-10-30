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

import math
from collections.abc import Mapping
from dataclasses import replace as _replace
from typing import Any

import numpy as np

from tss import sweep_law
from tss.constants import C_LIGHT, K_BOLTZ
from tss.interfaces import EchoSet, ReceivedSignal, SweepSchedule
from tss.l2_instrument.antenna import ShortDipole, antenna_noise_temperature_k
from tss.l2_instrument.impairments import TransmitImpairments
from tss.l2_instrument.noise import (
    is_extrapolated,
    noise_density_dbm_per_hz,
    sky_temperature,
)
from tss.l2_instrument.oscillator import Oscillator
from tss.l2_instrument.receive_chain import FrontEnd, LeakageHighpass
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

def _leakage_model_from_config(leak: Mapping[str, Any], isolation: Mapping[str, Any]) -> Any:
    """The leakage model a ``leakage`` block asks for.

    ``profile: mom_bus`` (the default when leakage is on, decision of
    2026-10-28) loads the moment-method isolation of the dipoles on the 6U
    bus, ``mom_geometry`` and ``mom_match`` choosing the column; ``flat`` takes
    ``isolation.transmit_receive_db``.  The coupling delay is
    ``isolation.coupling_delay_s`` either way.
    """
    from tss.l2_instrument import leakage as _lk

    profile = str(leak.get("profile", "mom_bus"))
    delay = float(isolation.get("coupling_delay_s", _lk.BASELINE_COUPLING_DELAY_S))
    if profile == "flat":
        return _lk.LeakageModel.from_config(isolation)
    if profile == "mom_bus":
        return _lk.mom_anchor_model(
            str(leak.get("mom_geometry", "bus_crossed_ends_offset")),
            str(leak.get("mom_match", "simultaneous")),
            coupling_delay_s=delay,
        )
    raise ValueError(f"leakage profile must be 'mom_bus' or 'flat', got {profile!r}")


#: Headroom the INS-11 noise rule keeps between the residual leakage and the
#: digitiser's clip level (Requirements Rev G sec. 5.9; decision G5b).
_INS11_HEADROOM_DB = 10.0


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
        impairments: TransmitImpairments | None = None,
        leakage: Any = None,
        leakage_enabled: bool = False,
        front_end: FrontEnd | None = None,
        front_end_enabled: bool | None = None,
        highpass: LeakageHighpass | None = None,
        highpass_enabled: bool | None = None,
        leakage_primed: bool = True,
        canceller_bound_db: float | None = 30.0,
        canceller: Any = None,
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
        self.impairments = impairments if impairments is not None else TransmitImpairments.none()
        # WP3 Phase 3.  The leakage is off unless asked for (decision of
        # 2026-10-28: off until Phase 4 scores it).  The front end and the
        # high-pass follow the leakage switch unless set explicitly.
        self.leakage_enabled = bool(leakage_enabled)
        if self.leakage_enabled and leakage is None:
            from tss.l2_instrument import leakage as _lk

            leakage = _lk.LeakageModel.from_config(self.isolation)
        self.leakage = leakage
        self.front_end = front_end if front_end is not None else FrontEnd()
        self.front_end_enabled = (
            self.leakage_enabled if front_end_enabled is None else bool(front_end_enabled)
        )
        self.highpass = highpass if highpass is not None else LeakageHighpass()
        self.highpass_enabled = (
            self.leakage_enabled if highpass_enabled is None else bool(highpass_enabled)
        )
        # WP3 Phase 4.  The transmitter runs between sweeps (operating
        # assumption of 2026-10-28), so the leakage is already present and the
        # high-pass already settled on it when a record starts.  The canceller
        # bound is the residual with the leakage taken down by that much ahead
        # of the receiver, an ideal canceller reported beside the score (WP5
        # owns the real one).
        self.leakage_primed = bool(leakage_primed)
        self.canceller_bound_db = None if canceller_bound_db is None else float(canceller_bound_db)
        # WP3 Phase 5: an RF canceller ahead of the front end
        # (tss.l2_instrument.leakage.RFCanceller), or None.
        self.canceller = canceller
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
        if "impairments" not in overrides:
            kwargs["impairments"] = TransmitImpairments.from_config(
                cfg.get("transmit_impairments"), oscillator=kwargs["oscillator"]
            )
        leak_cfg = dict(cfg.get("leakage") or {})
        kwargs["leakage_enabled"] = bool(leak_cfg.get("enabled", False))
        if kwargs["leakage_enabled"] and "leakage" not in overrides:
            kwargs["leakage"] = _leakage_model_from_config(leak_cfg, dict(cfg.get("isolation") or {}))
        fe_cfg = dict(cfg.get("front_end") or {})
        kwargs["front_end"] = FrontEnd.from_config(fe_cfg)
        kwargs["front_end_enabled"] = fe_cfg.get("enabled")
        hp_cfg = dict(cfg.get("leakage_highpass") or {})
        kwargs["highpass"] = LeakageHighpass.from_config(hp_cfg)
        kwargs["highpass_enabled"] = hp_cfg.get("enabled")
        kwargs["leakage_primed"] = bool(leak_cfg.get("primed", True))
        cb = leak_cfg.get("canceller_bound_db", 30.0)
        kwargs["canceller_bound_db"] = None if cb is None else float(cb)
        if "canceller" not in overrides:
            from tss.l2_instrument.leakage import RFCanceller

            kwargs["canceller"] = RFCanceller.from_config(cfg.get("canceller"))
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
        imp = self.impairments
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
            # The transmitted law carries the sweep's nonlinearity only when
            # the impairments say so (WP3 Phase 2); an independent reference
            # (INS-21 off) sweeps the ideal law.
            ideal = _replace(sweep, nonlinearity_ppm=0.0)
            law = sweep if imp.sweep_nonlinearity else ideal
            ref_law = law if imp.single_reference else ideal
            f_tx = sweep_law.frequency_hz(law, t)
            phi_rx = sweep_law.phase_rad(ref_law, t)

            def tx_frequency(x: np.ndarray) -> np.ndarray:
                return sweep_law.frequency_hz(law, x)

            def tx_phase(x: np.ndarray) -> np.ndarray:
                return sweep_law.phase_rad(law, x)

            def ref_phase_law(x: np.ndarray) -> np.ndarray:
                return sweep_law.phase_rad(ref_law, x)
        else:
            f_tx = f0 + direction * k * t

            def tx_frequency(x: np.ndarray) -> np.ndarray:
                return f0 + direction * k * x

            def tx_phase(x: np.ndarray) -> np.ndarray:
                return 2.0 * np.pi * (f0 * x + direction * k * x**2 / 2.0)

            ref_phase_law = tx_phase

        # Clutter rows (WP4 Phase 1, tss.l1_propagation.clutter) are synthesised
        # like echoes but never join a mode's dispersion trace: they are not on
        # it, and a trace with two delays at one frequency is not a trace.
        is_clutter = np.array([str(f).startswith("clutter:") for f in flags_in.ravel()[: freq.size]],
                              dtype=bool)
        traces = (
            _traces_by_mode(freq, delay, getattr(echoes, "mode", None), usable & ~is_clutter)
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

        # Transmitter impairments (WP3 Phase 2): one draw for the record, which
        # the leakage of the same record must share (Phase 3).
        realization = None
        ref_noise = None
        if imp.needs_realization:
            realization = imp.realize(
                sweep, fs_hz, n_samples, rng, max_delay_s=1.2 * max_delay + 1.0e-3
            )
            ref_noise = realization.ref_phase()
        dopp_in = getattr(echoes, "doppler_hz", None)
        dopp = (
            np.asarray(dopp_in, dtype=float)
            if imp.doppler and dopp_in is not None and np.shape(dopp_in) == freq.shape
            else None
        )
        rows_doppler = 0

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
            imp_phase: np.ndarray | float = 0.0
            imp_amp: np.ndarray | float = 1.0
            if realization is not None:
                imp_phase = ref_noise[lo:hi] - realization.tx_phase_window(tau, lo, hi)
                imp_amp = realization.tx_amplitude_window(tau, lo, hi)
            if dopp is not None and np.isfinite(dopp[row]):
                # A closing echo beats at f_R - f_d on an up-sweep and
                # -(f_R + f_d) on a down-sweep (tss.l1_propagation.doppler).
                imp_phase = imp_phase - 2.0 * np.pi * float(dopp[row]) * tt
                rows_doppler += 1
            trace = None if is_clutter[row] else traces.get(str(row_mode[row]))
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
                    ref_phase_law(tt)
                    - tx_phase(t_e)
                    - 2.0 * np.pi * f_e * tau_e
                    + 2.0 * np.pi * trace.phase_cycles(f_e)
                )
                live = t_e >= 0.0
                samples[lo:hi] += np.where(live, amplitude * imp_amp * np.exp(1j * (phase + imp_phase)), 0.0)
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
                phase = phi_rx[lo:hi] - tx_phase(tt - tau)
                live = tt >= tau
                samples[lo:hi] += np.where(live, amplitude * imp_amp * np.exp(1j * (phase + imp_phase)), 0.0)
                rows_used += 1
                continue
            else:
                phase = 2.0 * np.pi * (
                    f0 * tau + direction * k * tau * tt - direction * k * tau**2 / 2.0
                )
            samples[lo:hi] += amplitude * imp_amp * np.exp(1j * (phase + imp_phase))
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

        # Interference (WP4 Phase 1): carriers listed by the clutter generator,
        # added at the receiver input like the echoes, before the leakage, the
        # noise and the receive chain (tss.l2_instrument.interference).
        carriers = getattr(echoes, "interference", None)
        carrier_info: dict[str, Any] = {"carriers_listed": 0, "carriers_rendered": 0,
                                        "carrier_power_dbm_max": None}
        if carriers:
            from tss.l2_instrument.interference import render_carriers

            if stepped:
                ref_phase_t = phi_rx
            else:
                ref_phase_t = ref_phase_law(t)
            x_c, carrier_info = render_carriers(
                carriers, phi_ref_rad=ref_phase_t, f_ref_hz=f_tx, t_s=t, fs_hz=fs_hz,
                gain_dbi=self.gain_dbi, polarization_loss_db=self.polarization_loss_db,
            )
            samples = samples + x_c

        # Leakage (WP3 Phase 3).  Added at the receiver input, before the noise
        # and the receive chain, and drawn from the record's impairment
        # realization so that it and the echoes share one oscillator.
        leak = None
        leak_uncancelled = None
        g_canc = None
        if self.leakage_enabled and self.leakage is not None:
            from tss.l2_instrument import leakage as _lk

            leak = _lk.leakage_signal(
                sweep, fs_hz, self.leakage, rng, n_samples=n_samples, t0_s=0.0,
                transmit_power_w=self.tx_power_w,
                tx_chain_efficiency_db=self.tx_chain_efficiency_db,
                impairments=imp, realization=realization, primed=self.leakage_primed,
            )
            leak_uncancelled = leak
            g_canc = None
            if self.canceller is not None:
                # WP3 Phase 5: the RF canceller subtracts a weighted sample of
                # the transmitted signal, so the leakage keeps its noise in
                # proportion; it adds its own noise when it is active.
                g_canc = self.canceller.residual_factor(sweep, t, self.leakage, rng)
                leak = leak * g_canc
                c_noise = self.canceller.added_noise(n_samples, fs_hz, rng)
                if c_noise is not None:
                    leak = leak + c_noise
            samples = samples + leak

        # Noise.  The antenna temperature follows the sweep, because the sky is
        # three orders of magnitude colder at the top of the band than at the
        # bottom and a single band-averaged figure would be wrong at both ends.
        t_ant = self.antenna_temperature_k(np.clip(f_tx, 1.0, None) / 1.0e6)
        t_sys = t_ant + self.receiver.noise_temperature_k
        sigma = np.sqrt(K_BOLTZ * t_sys * fs_hz / 2.0)
        samples = samples + rng.normal(0.0, sigma) + 1j * rng.normal(0.0, sigma)

        # Receive chain.  Front end (bandpass-equivalent compression of the
        # composite), the leakage high-pass after the de-chirp (WP3 Phase 3),
        # then gain, clipping against the instantaneous dynamic range referred
        # to the median noise floor, and quantisation.
        hp_state = None
        if leak is not None and self.leakage_primed and self.highpass_enabled:
            hp_state = self.highpass.steady_state(self._front(samples[:1])[0], sweep, fs_hz)
        amplified = self.receiver.apply_gain(self._analog(samples, sweep, fs_hz, state=hp_state))
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
            chain = self.receiver.apply_gain(
                self._analog(samples + extra_noise, sweep, fs_hz, state=hp_state)
            )
            chain, _ = self.receiver.clip(chain, clip_level_w)
            rows_out.append(
                self.receiver.quantise(chain, clip_level_w) if self.quantise else chain
            )
        block = np.vstack(rows_out)

        # Requirement evidence owned by L2 (WP3 Phase 4): the declared values
        # always, and the leakage scores when the leakage is in the record.
        l2_metrics = self._declared_metrics()
        per_step: dict[str, Any] | None = None
        if leak is not None:
            leak_metrics, per_step = self._leakage_metrics(
                leak, sweep, fs_hz, t_sys, clip_level_w, uncancelled=leak_uncancelled,
                canceller_factor=g_canc,
            )
            l2_metrics.update(leak_metrics)

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
            "transmit_impairments": dict(imp.describe()),
            "rows_doppler_applied": rows_doppler,
            # WP4 Phase 1: clutter rows synthesised, and carriers at the input.
            "clutter_rows_used": int(np.count_nonzero(usable & is_clutter)),
            "interference": carrier_info,
            "leakage_modelled": leak is not None,
            "leakage_note": (
                "Transmit-to-receive leakage is in these samples (WP3 Phase 3), "
                f"profile: {getattr(self.leakage, 'label', '') or 'flat'}, "
                f"{'primed (transmitter running between sweeps)' if self.leakage_primed else 'cold start'}.  "
                "INS-11, INS-23 and INS-24 are scored from it (l2_metrics, WP3 Phase 4)."
                if leak is not None
                else "Leakage is off in this configuration, so the receiver never sees it "
                "and INS-11, INS-23 and INS-24 are not evaluated from this run."
            ),
            "leakage_power_dbm_mean": (
                float(w_to_dbm(float(np.mean(np.abs(leak) ** 2)))) if leak is not None else None
            ),
            "front_end": {"enabled": self.front_end_enabled, **dict(self.front_end.describe())},
            "leakage_highpass": {"enabled": self.highpass_enabled, **dict(self.highpass.describe())},
            "leakage_primed": bool(leak is not None and self.leakage_primed),
            "rf_canceller": (
                dict(self.canceller.describe()) if (leak is not None and self.canceller is not None) else None
            ),
            "l2_metrics": l2_metrics,
            "leakage_per_step": per_step,
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

    def _front(self, x: np.ndarray) -> np.ndarray:
        if self.front_end_enabled and not self.front_end.linear:
            return self.front_end.apply(x)
        return np.asarray(x, dtype=complex)

    def _analog(
        self, x: np.ndarray, sweep: SweepSchedule, fs_hz: float, *, state: np.ndarray | None = None
    ) -> np.ndarray:
        """The analog chain ahead of the gain stage: front end, then high-pass."""
        x = self._front(x)
        if self.highpass_enabled:
            x = self.highpass.apply(x, sweep, fs_hz, state=state)
        return x

    # -- requirement evidence (WP3 Phase 4) ------------------------------------

    def _declared_metrics(self) -> dict[str, Any]:
        """INS-14 and INS-26 from the configuration, labelled as declared.

        Neither is measured by a simulated record: the transmit power is what
        the configuration says, and the Allan deviation is the oscillator part's
        (decision of 2026-10-28, "score with evidence").  INS-26 is one number
        scored against a threshold written at 0.05 s and a goal written at
        1.15 s, so the composite is the Allan deviation at the interval of the
        highest level met: at 1.15 s when it meets the goal, otherwise at 0.05 s
        (held just above the goal, so that a part that meets the goal at 0.05 s
        and not at 1.15 s scores threshold, not goal).
        """
        from tss.l2_instrument import oscillator as _osc

        s_t = float(np.asarray(self.oscillator.sigma_y(_osc.INS26_THRESHOLD_INTERVAL_S)))
        s_g = float(np.asarray(self.oscillator.sigma_y(_osc.INS26_GOAL_INTERVAL_S)))
        if s_t > _osc.INS26_ALLAN_THRESHOLD:
            comp = s_t
        elif s_g <= _osc.INS26_ALLAN_GOAL:
            comp = s_g
        else:
            comp = max(s_t, float(np.nextafter(_osc.INS26_ALLAN_GOAL, np.inf)))
        am = getattr(self.impairments, "am_noise_dbc_hz", None)
        return {
            **(
                {"transmit_am_noise_dbc_hz": float(am),
                 "transmit_am_noise_dbc_hz_basis": "declared: the configured AM noise allocation"}
                if am is not None else {}
            ),
            "transmit_power_w": float(self.tx_power_w),
            "transmit_power_w_basis": "declared: the configured transmit power",
            "allan_deviation_composite": comp,
            "allan_deviation_composite_basis": (
                "declared: the configured oscillator part, sigma_y at 0.05 s (threshold) "
                "and 1.15 s (goal)"
            ),
            "allan_deviation_at_0p05s": s_t,
            "allan_deviation_at_1p15s": s_g,
        }

    def _leakage_metrics(
        self,
        leak: np.ndarray,
        sweep: SweepSchedule,
        fs_hz: float,
        t_sys: np.ndarray,
        clip_level_w: float,
        *,
        uncancelled: np.ndarray | None = None,
        canceller_factor: np.ndarray | None = None,
    ) -> tuple[dict[str, Any], dict[str, Any] | None]:
        """INS-11, INS-23 and INS-24 from the leakage of this record.

        The residual is the leakage alone through the analog chain: the front
        end at the leakage's own level (it dominates the composite, so it sets
        the compression) and the high-pass, started settled when the
        transmitter runs between sweeps.  It carries everything the transmitter
        put on the leakage that the reference does not take off again.

        INS-24 (decision G5a) is the margin of that residual below the weakest
        echo INS-10 admits, in each step's echo band; the margin below the
        free-space echo from R_max and the margin behind an ideal canceller of
        ``canceller_bound_db`` are reported beside it.  INS-23 is the least
        isolation of the leakage model over the sweep.  INS-11 is the configured
        dynamic range, verified by the record: it is reduced by any overload of
        the digitiser by the residual and by any step whose high-pass rejection
        falls short of the noise rule (decision G5b), so a range that the
        record does not support is not reported.
        """
        from tss.l2_instrument import leakage as _lk
        from tss.l2_instrument.reference_echo import reference_levels

        def chain(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            x = self._front(x)
            if not self.highpass_enabled:
                return x, x
            st = self.highpass.steady_state(x[0], sweep, fs_hz) if self.leakage_primed else None
            return x, self.highpass.apply(x, sweep, fs_hz, state=st)

        leak_in, res = chain(leak)
        stepped = sweep_law.is_stepped(sweep)
        if stepped:
            f_eval = np.asarray(sweep_law.step_centres_hz(sweep), dtype=float)
        else:
            f_eval = np.linspace(float(sweep.f_start_hz), float(sweep.f_stop_hz), 256)

        out: dict[str, Any] = {}
        per_step: dict[str, Any] | None = None
        if stepped:
            levels = reference_levels(self, sweep)
            resid = _lk.dechirped_residual(res, sweep, fs_hz, levels=levels)
            canc = None
            if self.canceller_bound_db is not None:
                base = leak if uncancelled is None else uncancelled
                _, res_c = chain(base * 10.0 ** (-self.canceller_bound_db / 20.0))
                canc = _lk.dechirped_residual(res_c, sweep, fs_hz, levels=levels)
            sc = _lk.score_ins23_ins24(
                f_eval, self.leakage, resid, levels=levels, canceller_residual=canc
            )
            per_step = {
                "f_hz": np.asarray(levels["f_hz"], dtype=float).tolist(),
                "residual_level_dbw": np.asarray(resid["level_dbw"]).tolist(),
                "residual_pedestal_dbw": np.asarray(resid["pedestal_dbw"]).tolist(),
                "residual_line_dbw": np.asarray(resid["line_dbw"]).tolist(),
                "weakest_echo_dbw": np.asarray(levels["weakest_echo_dbw"]).tolist(),
                "echo_rmax_dbw": np.asarray(levels["echo_rmax_dbw"]).tolist(),
                "ins24_margin_db": np.asarray(sc.pop("ins24_per_step_margin_db")).tolist(),
                "ins24_margin_db_rmax": np.asarray(sc.pop("ins24_per_step_margin_db_rmax")).tolist(),
                "isolation_db": np.asarray(_lk.isolation_profile_db(f_eval, self.leakage)).tolist(),
            }
            if canc is not None:
                per_step["residual_level_dbw_canceller"] = np.asarray(canc["level_dbw"]).tolist()
            # INS-29 (WP3 Phase 5): the noise the leakage brings, against the
            # system noise in the same cell.
            over = np.asarray(resid["pedestal_dbw"]) - np.asarray(levels["noise_cell_dbw"], dtype=float)
            per_step["leakage_noise_over_system_noise_db"] = over.tolist()
            if np.any(np.isfinite(over)):
                out["leakage_noise_over_system_noise_db"] = float(np.nanmax(over))
                out["leakage_noise_over_system_noise_freq_hz"] = float(
                    np.asarray(levels["f_hz"])[int(np.nanargmax(over))]
                )
        else:
            sc = _lk.score_ins23_ins24(f_eval, self.leakage)
            out["ins24_not_scored_reason"] = "INS-24 is scored per step and needs a stepped sweep"

        out["isolation_db_min"] = sc["isolation_db_min"]
        out["isolation_db_min_freq_hz"] = sc["isolation_db_min_freq_hz"]
        for key, name in (
            ("ins24_margin_db", "leakage_residual_margin_db"),
            ("ins24_worst_freq_hz", "leakage_residual_margin_freq_hz"),
            ("ins24_margin_db_rmax", "leakage_residual_margin_db_rmax"),
            ("ins24_steps_failing", "leakage_residual_steps_failing"),
            ("ins24_steps_scored", "leakage_residual_steps_scored"),
            ("ins24_margin_db_canceller", "leakage_residual_margin_db_canceller"),
            ("ins24_steps_failing_canceller", "leakage_residual_steps_failing_canceller"),
        ):
            if sc.get(key) is not None:
                out[name] = sc[key]
        if self.canceller_bound_db is not None and stepped:
            out["leakage_canceller_bound_db"] = self.canceller_bound_db
        if self.canceller is not None:
            out["rf_canceller_nominal_depth_db"] = float(self.canceller.nominal_depth_db())
            if canceller_factor is not None and stepped:
                # INS-31: the depth the canceller achieved at its worst step,
                # the carrier it leaves averaged over each dwell.
                dwell = float(sweep.step_dwell_s)
                n_st = sweep_law.step_count(sweep)
                idx = np.clip((np.arange(canceller_factor.size) / float(fs_hz) / dwell).astype(int), 0, n_st - 1)
                g2 = np.bincount(idx, weights=np.abs(canceller_factor) ** 2, minlength=n_st)
                cnt = np.bincount(idx, minlength=n_st)
                depth = -10.0 * np.log10(np.maximum(g2[cnt > 0] / cnt[cnt > 0], 1e-30))
                # Rounded to a micro-decibel, so that an ideal canceller at a
                # threshold scores at it rather than a rounding error below.
                out["rf_canceller_min_depth_db"] = round(float(np.min(depth)), 6)
                out["rf_canceller_min_depth_db_basis"] = (
                    "measured in the record: the carrier left per dwell; the continuity part of "
                    "INS-31 is not scored here"
                )

        # INS-11: the configured range, verified by the record.
        dr = float(self.receiver.instantaneous_dynamic_range_db)
        nfs_t = K_BOLTZ * np.asarray(t_sys, dtype=float) * float(fs_hz)
        shortfall = -np.inf
        worst_f = float("nan")
        for s, e, _fc in self.highpass.segments(sweep, 0.0, res.size, fs_hz):
            if e - s < 8:
                continue
            tail = slice(s + int(0.75 * (e - s)), e)
            lin = float(np.mean(np.abs(leak_in[tail]) ** 2))
            lout = max(float(np.mean(np.abs(res[tail]) ** 2)), 1e-300)
            nfs = float(np.median(nfs_t[s:e]))
            need = 10.0 * math.log10(lin / nfs) - dr + _INS11_HEADROOM_DB
            got = 10.0 * math.log10(lin / lout)
            if need - got > shortfall:
                shortfall = need - got
                worst_f = float(sweep_law.frequency_hz(sweep, (s + e) / 2.0 / float(fs_hz)))
        clip_in = float(clip_level_w) / float(undb(self.receiver.gain_db))
        over = 10.0 * math.log10(max(float(np.max(np.abs(res) ** 2)), 1e-300) / clip_in)
        out["dynamic_range_db_at_digitiser"] = dr - max(0.0, float(shortfall), over)
        out["dynamic_range_db_at_digitiser_basis"] = (
            "configured, verified by the record: residual leakage against the clip "
            "level and the high-pass rejection against the noise rule at every step"
        )
        out["ins11_configured_db"] = dr
        out["ins11_rejection_margin_db"] = float(-shortfall)
        out["ins11_rejection_margin_freq_hz"] = worst_f
        out["ins11_residual_peak_over_clip_db"] = over
        return out, per_step

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
            "transmit_impairments": dict(self.impairments.describe()),
            "resources": dict(self.resources),
            "leakage_modelled": bool(self.leakage_enabled),
            "leakage_primed": bool(self.leakage_primed),
            "canceller_bound_db": self.canceller_bound_db,
            "rf_canceller": None if self.canceller is None else dict(self.canceller.describe()),
            "leakage": None if self.leakage is None else dict(self.leakage.describe()),
            "front_end": {"enabled": self.front_end_enabled, **dict(self.front_end.describe())},
            "leakage_highpass": {"enabled": self.highpass_enabled, **dict(self.highpass.describe())},
            "matching_network_loss": (
                "carried in the antenna radiation efficiency only (decision of "
                "2026-10-17); tx_chain_efficiency_db excludes it"
            ),
        }
