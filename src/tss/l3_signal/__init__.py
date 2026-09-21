"""Layer: L3
WP: WP4, WP5, WP11
Claims: C1, C2, C3, C4, C5
Status: partial
Requirements: INS-01, INS-04, INS-05, INS-07, INS-09, INS-10

L3, signal processing: de-chirp against the transmitted sweep, transform, map
bin index to apparent range, extract the trace and form the ionogram.

The package is thin, as ``CONTRIBUTING.md`` requires of a layer package.  It
imports its tier modules so that registration happens, and re-exports the names
another layer or a study legitimately needs.  There is no signal processing
here.

What is delivered and what is not:

======================  ========  ==========================================
Module                  Status    Note
======================  ========  ==========================================
``base``                done      Sweep law, phase, beat-to-range mapping.
``dechirp``             done      De-chirp and the range-correlation
                                  treatment of waveform purity.
``transform``           done      Windowed transform and the per-cell noise
                                  bandwidth that makes C1 work.
``trace_extract``       partial   Threshold and connect, valid on clean
                                  data.  The robust extractor is WP4.
``ionogram``            partial   Registered as ``l3/baseline``.  No O and X
                                  labelling, no Doppler.
``vector``              stub      Three-channel separation and arrival
                                  angle, WP11, C5, goal level.
``motion_separated``    stub      Registered as ``l3/motion_separated``,
                                  the tier c3_doppler.yaml names.  WP5.
``novel/``              stub      The four claimed algorithms.  See its own
                                  docstring.
======================  ========  ==========================================

A note on module names.  ``docs/architecture.md`` lists this layer's modules as
``detect.py``, ``doppler_decouple.py``, ``leakage_suppress.py`` and ``aoa.py``,
and lists the tiers as ``baseline``, ``updown``, ``leakage_suppress`` and
``aoa``.  The delivered layout is the one above, with the claimed algorithms
grouped under ``novel/`` so that the boundary between conventional processing
and the contribution is visible in the directory listing, and with the tier the
campaign actually names spelled ``motion_separated``.  The architecture document
is stale on this point and, by its own rule, the code is right and the document
needs the correction.

Both stubs are registered or importable so that a campaign naming them can be
validated before they exist.  They raise :class:`NotImplementedError` when
called, with a message naming the work package.
"""

from __future__ import annotations

# Importing the tier modules is what populates the registry.  ``dechirp`` and
# ``transform`` are deliberately absent from this list: each is a module whose
# principal export has the same name, and the function is what a caller wants
# from ``from tss.l3_signal import dechirp``.  Both modules are imported
# transitively by the function imports below and remain reachable as
# ``tss.l3_signal.dechirp`` and ``tss.l3_signal.transform`` through the ordinary
# ``import`` statement.
from tss.l3_signal import (  # noqa: F401
    base,
    ionogram,
    motion_separated,
    novel,
    trace_extract,
    vector,
)
from tss.l3_signal.base import (  # noqa: F401
    SignalProcessingBase,
    beat_from_range_m,
    instantaneous_frequency_hz,
    range_from_beat_hz,
    sweep_bandwidth_hz,
    sweep_duration_s,
    sweep_phase_rad,
)
from tss.l3_signal.dechirp import (  # noqa: F401
    Dechirped,
    dechirp,
    phase_residual_beta_rad,
    range_correlation_suppression,
    range_correlation_suppression_db,
    reference_waveform,
    residual_phase_noise_dbc_hz,
    two_way_delay_s,
)
from tss.l3_signal.ionogram import (  # noqa: F401
    UNRESOLVED_MODE,
    BaselineSignalProcessing,
)
from tss.l3_signal.motion_separated import MotionSeparatedSignalProcessing  # noqa: F401
from tss.l3_signal.trace_extract import (  # noqa: F401
    TraceSample,
    connect_trace,
    extract_trace,
    leading_edge,
    noise_floor_db,
    robust_trace_extract,
)
from tss.l3_signal.transform import (  # noqa: F401
    WINDOWS,
    RangeSpectrum,
    cell_noise_bandwidth_hz,
    range_bin_index,
    transform,
    unambiguous_range_m,
    window_gains,
    window_taper,
)

__all__ = [
    # tiers
    "BaselineSignalProcessing",
    "MotionSeparatedSignalProcessing",
    "SignalProcessingBase",
    # sweep law and mapping
    "sweep_bandwidth_hz",
    "sweep_duration_s",
    "instantaneous_frequency_hz",
    "sweep_phase_rad",
    "beat_from_range_m",
    "range_from_beat_hz",
    # de-chirp
    "Dechirped",
    "dechirp",
    "reference_waveform",
    "two_way_delay_s",
    "range_correlation_suppression",
    "range_correlation_suppression_db",
    "phase_residual_beta_rad",
    "residual_phase_noise_dbc_hz",
    # transform
    "RangeSpectrum",
    "transform",
    "window_taper",
    "window_gains",
    "range_bin_index",
    "unambiguous_range_m",
    "cell_noise_bandwidth_hz",
    "WINDOWS",
    # trace extraction
    "TraceSample",
    "leading_edge",
    "extract_trace",
    "connect_trace",
    "noise_floor_db",
    "robust_trace_extract",
    # constants
    "UNRESOLVED_MODE",
]
