"""Layer: L2
WP: WP3
Claims: C1, C4, C5
Status: partial
Requirements: INS-01, INS-11, INS-14, INS-23, INS-24, INS-26

L2, instrument: everything between the echoes arriving at the spacecraft and
complex baseband samples leaving the receiver.

The package is thin, as ``CONTRIBUTING.md`` requires of a layer package.  It
imports its tier modules so that registration happens, and re-exports the names
another layer or a study legitimately needs.  There is no physics here.

What is delivered and what is not:

======================  ========  ==========================================
Module                  Status    Note
======================  ========  ==========================================
``linkbudget``          done      The migrated closed-form budget with its
                                  own tests.  A known seam: kept whole and
                                  not split along layer lines, see
                                  ``docs/architecture.md``.
``antenna``             done      Electrically short dipole, directivity and
                                  efficiency carried separately.
``noise``               partial   Franklin and MacLean as primary; the Cane
                                  and Dulk galactic background as the WP12
                                  candidate below 0.6 MHz.
``waveform``            done      Sweep law and its exact phase integral,
                                  with an optional nonlinearity.
``receiver``            done      Gain, noise figure, clipping, quantisation.
``oscillator``          done      Allan deviation to phase noise and the
                                  disputed range-correlation credit.
``baseline``            partial   The registered tier, ``l2/baseline``.
                                  Partial because leakage is missing.
``leakage``             stub      Transmit-to-receive coupling and its
                                  de-chirped residual.  WP3, critical path.
======================  ========  ==========================================

The stub is importable so that a campaign naming it can be validated before it
exists.  It raises :class:`NotImplementedError` when called, with a message
naming the work package.

THE SEAM
--------
``linkbudget`` contains terms that belong to L0 and L1 as well as to L2, and it
is kept whole because its 21 passing tests are written against its current
signatures.  The full reasoning, the cost and the scheduled split are in
``docs/architecture.md``, section "Known seams".
"""

from __future__ import annotations

from tss.l2_instrument import (  # noqa: F401
    antenna,
    baseline,
    leakage,
    linkbudget,
    noise,
    oscillator,
    receiver,
    waveform,
)
from tss.l2_instrument.antenna import (  # noqa: F401
    ShortDipole,
    antenna_noise_temperature_k,
    non_sky_temperature_k,
    resonant_length_m,
    short_dipole_efficiency,
    sky_dominance_margin_db,
)
from tss.l2_instrument.baseline import BaselineInstrument  # noqa: F401
from tss.l2_instrument.leakage import (  # noqa: F401
    LeakageModel,
    dechirped_residual,
    isolation_profile_db,
    leakage_signal,
    score_ins23_ins24,
)
from tss.l2_instrument.linkbudget import (  # noqa: F401
    Antenna,
    BudgetResult,
    IonosphereProfile,
    SounderConfig,
    TabulatedProfile,
    absorption_db,
    legacy_snr_ionosonde,
    link_budget,
    reflection_height,
    required_tx_power,
    snr_vs_frequency,
)
from tss.l2_instrument.noise import (  # noqa: F401
    cane_dulk_brightness,
    cane_dulk_temperature_k,
    describe_reference,
    is_extrapolated,
    noise_density_dbm_per_hz,
    sky_temperature,
)
from tss.l2_instrument.oscillator import (  # noqa: F401
    Oscillator,
    phase_jitter_rad,
    range_correlation_credit,
    range_correlation_credit_db,
    range_correlation_suppression,
    sideband_power_ratio,
)
from tss.l2_instrument.receiver import Receiver  # noqa: F401
from tss.l2_instrument.waveform import (  # noqa: F401
    beat_frequency_hz,
    generate_sweep,
    schedule_from_config,
    sweep_bandwidth_hz,
    sweep_duration_s,
    sweep_frequency_hz,
    sweep_phase_rad,
)

__all__ = [
    # tier
    "BaselineInstrument",
    # link budget, kept whole; see the seam note above
    "SounderConfig",
    "BudgetResult",
    "IonosphereProfile",
    "TabulatedProfile",
    "Antenna",
    "link_budget",
    "snr_vs_frequency",
    "required_tx_power",
    "absorption_db",
    "reflection_height",
    "legacy_snr_ionosonde",
    # antenna
    "ShortDipole",
    "short_dipole_efficiency",
    "resonant_length_m",
    "antenna_noise_temperature_k",
    "non_sky_temperature_k",
    "sky_dominance_margin_db",
    # noise
    "sky_temperature",
    "cane_dulk_brightness",
    "cane_dulk_temperature_k",
    "noise_density_dbm_per_hz",
    "is_extrapolated",
    "describe_reference",
    # waveform
    "sweep_frequency_hz",
    "sweep_phase_rad",
    "generate_sweep",
    "sweep_duration_s",
    "sweep_bandwidth_hz",
    "beat_frequency_hz",
    "schedule_from_config",
    # receive chain
    "Receiver",
    # oscillator
    "Oscillator",
    "phase_jitter_rad",
    "sideband_power_ratio",
    "range_correlation_suppression",
    "range_correlation_credit",
    "range_correlation_credit_db",
    # WP3 stub
    "LeakageModel",
    "isolation_profile_db",
    "leakage_signal",
    "dechirped_residual",
    "score_ins23_ins24",
]
