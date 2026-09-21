"""Layer contracts for the L0 to L5 simulator.

This module is the architectural keystone of the package.  It defines the data
objects that pass between tiers and the protocol each tier implementation must
satisfy.  Nothing in here computes physics.  The intent is that any tier can be
replaced by a different fidelity implementation without disturbing the rest,
which is the property Section 3.4 of the proposal asks for.

The layers are:

    L0  Ionosphere state      truth Ne, collision frequency, geomagnetic field
    L1  Propagation           echo generation from the truth state
    L2  Instrument            waveform, antenna, receive chain, impairments
    L3  Signal processing     de-chirp, transform, ionogram formation
    L4  Retrieval             virtual height to true height inversion
    L5  Assessment            score against L0 truth and against requirements

Every tier implementation registers itself with :mod:`tss.registry` under a
name, and a campaign configuration selects tiers by those names.  A tier must
also return a JSON-serialisable ``describe()`` so the run manifest records
exactly which implementation produced a result.

Contract rules that apply to every tier:

1.  Inputs are never mutated.  A tier returns new objects.
2.  Units are SI unless the attribute name carries a suffix (``_km``, ``_mhz``,
    ``_db``).  See :mod:`tss.units`.
3.  A tier that cannot answer for a given input returns a result with the
    relevant entries masked and a reason recorded in ``flags``, rather than
    raising.  Physical impossibility is a result, not an error.  The band edges
    are the motivating case: below the local plasma frequency at the spacecraft
    the wave is evanescent and there is no measurement to report.
4.  Randomness is drawn only from a ``numpy.random.Generator`` passed in, so a
    run is reproducible from its recorded seed.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Protocol, runtime_checkable

import numpy as np

__all__ = [
    "IonosphereState",
    "PlatformState",
    "EchoSet",
    "ReceivedSignal",
    "Ionogram",
    "DensityProfile",
    "Assessment",
    "L0Ionosphere",
    "L1Propagation",
    "L2Instrument",
    "L3SignalProcessing",
    "L4Retrieval",
    "L5Assessment",
    "Tier",
]


# ---------------------------------------------------------------------------
# Data objects
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class IonosphereState:
    """Truth ionosphere at one epoch and location.  Output of L0.

    Attributes
    ----------
    alt_m
        Monotonically increasing altitude grid, metres above the ellipsoid.
    ne
        Electron density on ``alt_m``, m^-3.
    nu
        Electron-neutral collision frequency on ``alt_m``, s^-1.  Absorption
        needs this, so a backend that cannot supply it should derive one from
        an atmosphere model rather than leave it None.
    te
        Electron temperature on ``alt_m``, K.  Optional.
    b_field_t
        Geomagnetic flux density at each altitude, tesla, as (3, N) in local
        north-east-down components.  Optional, but the magnetoionic splitting
        of L1 and the X-mode cutoff of the band-edge study both need it.
    lat_deg, lon_deg
        Geodetic location the profile applies to.
    epoch
        ISO 8601 UTC timestamp.
    source
        Free text naming the model and its options, e.g.
        ``"IRI-2020 topside=NeQuick R12=105"``.  Carried into the manifest.
    extra
        Scalar diagnostics a backend chooses to expose, such as ``nmf2``,
        ``hmf2_km``, ``solar_zenith_deg``.  Studies may read these but no tier
        may depend on a key being present.
    """

    alt_m: np.ndarray
    ne: np.ndarray
    nu: np.ndarray | None = None
    te: np.ndarray | None = None
    b_field_t: np.ndarray | None = None
    lat_deg: float = 0.0
    lon_deg: float = 0.0
    epoch: str = ""
    source: str = ""
    extra: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PlatformState:
    """Spacecraft position, velocity and attitude at one instant.

    C3 rests on the platform's motion being large and precisely known, so this
    object is a first-class input to L1 and to the novel L3 algorithms rather
    than a detail buried inside an orbit propagator.
    """

    alt_m: float
    lat_deg: float
    lon_deg: float
    velocity_ecef_mps: np.ndarray | None = None
    attitude_quat: np.ndarray | None = None
    epoch: str = ""


@dataclass(frozen=True)
class EchoSet:
    """Echoes returned to the platform.  Output of L1.

    One entry per (sounding frequency, magnetoionic mode, ray).  Arrays are
    parallel and one-dimensional so that an analytic tier and a ray-traced tier
    can produce the same object with different row counts.

    ``mode`` is ``"O"`` or ``"X"``.  ``flags`` carries per-row physical states,
    notably ``"evanescent"`` (below local plasma frequency) and
    ``"penetrating"`` (above foF2), which are results and not failures.
    """

    freq_hz: np.ndarray
    mode: np.ndarray
    group_delay_s: np.ndarray
    amplitude_db: np.ndarray
    absorption_db: np.ndarray
    doppler_hz: np.ndarray | None = None
    aoa_deg: np.ndarray | None = None
    reflection_alt_m: np.ndarray | None = None
    flags: np.ndarray | None = None
    source: str = ""


@dataclass(frozen=True)
class ReceivedSignal:
    """Complex baseband at the receiver.  Output of L2.

    Attributes
    ----------
    samples
        Complex samples, shape ``(n_channels, n_samples)``.  A scalar receive
        configuration has one channel; the vector configuration of C5 has
        three.
    fs_hz
        Sample rate.
    sweep
        The sweep schedule actually transmitted, so L3 de-chirps against what
        was sent rather than against what was requested.  INS-21's
        single-reference architecture means the same object is the de-chirp
        reference, which is what makes the range-correlation credit of
        Technical Volume Chapter 5 legitimate.
    noise_floor_dbm_hz
        The external noise density assumed, for traceability back to the
        Franklin and MacLean reference.
    """

    samples: np.ndarray
    fs_hz: float
    sweep: SweepSchedule
    noise_floor_dbm_hz: np.ndarray | None = None
    flags: Mapping[str, Any] = field(default_factory=dict)
    source: str = ""


@dataclass(frozen=True)
class SweepSchedule:
    """The frequency-versus-time law the instrument transmits.

    ``f_start_hz`` to ``f_stop_hz`` at rate ``k_f_hz_per_s``, repeated with
    period ``sweep_period_s``.  ``direction`` is +1 for up-sweep and -1 for
    down-sweep; the up/down pairing of C3 needs both, so a campaign supplies a
    sequence of these rather than one.
    """

    f_start_hz: float
    f_stop_hz: float
    k_f_hz_per_s: float
    sweep_period_s: float
    direction: int = 1
    nonlinearity_ppm: float = 0.0


@dataclass(frozen=True)
class Ionogram:
    """Virtual range against frequency.  Output of L3.

    ``virtual_range_m`` has shape ``(n_freq, n_modes)`` with NaN where no trace
    was detected.  ``snr_db`` is parallel to it.  ``doppler_mps`` is populated
    only when the L3 implementation resolved it, which is the point of C3.
    """

    freq_hz: np.ndarray
    virtual_range_m: np.ndarray
    snr_db: np.ndarray
    modes: Sequence[str] = ("O", "X")
    doppler_mps: np.ndarray | None = None
    aoa_deg: np.ndarray | None = None
    residual_range_bias_m: np.ndarray | None = None
    flags: Mapping[str, Any] = field(default_factory=dict)
    source: str = ""


@dataclass(frozen=True)
class DensityProfile:
    """Retrieved electron density against true height.  Output of L4."""

    alt_m: np.ndarray
    ne: np.ndarray
    ne_sigma: np.ndarray | None = None
    alt_sigma_m: np.ndarray | None = None
    flags: Mapping[str, Any] = field(default_factory=dict)
    source: str = ""


@dataclass(frozen=True)
class Assessment:
    """Scored result.  Output of L5.

    ``requirements`` maps a requirement ID (``"INS-01"``, ``"SCI-09"``) to a
    record carrying the achieved value, the threshold and goal it was scored
    against, and the verdict.  ``errors`` carries the truth-versus-retrieval
    comparison.  ``configuration`` names the configuration the result belongs
    to, because Section 1.3 of the proposal reports per configuration and a
    result without its configuration is not a result.
    """

    configuration: str
    requirements: Mapping[str, Mapping[str, Any]]
    errors: Mapping[str, float] = field(default_factory=dict)
    flags: Mapping[str, Any] = field(default_factory=dict)
    source: str = ""


# ---------------------------------------------------------------------------
# Tier protocols
# ---------------------------------------------------------------------------


@runtime_checkable
class Tier(Protocol):
    """Common part of every tier: it can describe itself for the manifest."""

    name: str

    def describe(self) -> Mapping[str, Any]:
        """Return a JSON-serialisable record of this tier and its settings."""
        ...


@runtime_checkable
class L0Ionosphere(Tier, Protocol):
    """Generates truth ionosphere states."""

    def state_at(
        self,
        lat_deg: float,
        lon_deg: float,
        epoch: str,
        alt_m: np.ndarray,
    ) -> IonosphereState: ...


@runtime_checkable
class L1Propagation(Tier, Protocol):
    """Turns a truth state plus a platform state into echoes."""

    def echoes(
        self,
        state: IonosphereState,
        platform: PlatformState,
        freq_hz: np.ndarray,
        rng: np.random.Generator,
    ) -> EchoSet: ...


@runtime_checkable
class L2Instrument(Tier, Protocol):
    """Applies the instrument to echoes and returns received samples."""

    def receive(
        self,
        echoes: EchoSet,
        sweep: SweepSchedule,
        rng: np.random.Generator,
    ) -> ReceivedSignal: ...


@runtime_checkable
class L3SignalProcessing(Tier, Protocol):
    """Forms an ionogram from received samples."""

    def form_ionogram(
        self,
        signal: ReceivedSignal,
        platform: PlatformState | None = None,
    ) -> Ionogram: ...


@runtime_checkable
class L4Retrieval(Tier, Protocol):
    """Inverts an ionogram to a density profile."""

    def invert(
        self,
        ionogram: Ionogram,
        platform: PlatformState,
    ) -> DensityProfile: ...


@runtime_checkable
class L5Assessment(Tier, Protocol):
    """Scores a retrieval against truth and against the requirement set."""

    def score(
        self,
        retrieved: DensityProfile,
        truth: IonosphereState,
        ionogram: Ionogram,
        configuration: str,
    ) -> Assessment: ...
