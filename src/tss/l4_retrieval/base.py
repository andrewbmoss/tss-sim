"""Layer: L4
WP: WP6
Claims: C2
Status: implemented
Requirements: -

Shared helpers for the L4 retrieval tiers.

Every L4 tier starts from the same three things: a trace read off an ionogram,
a spacecraft altitude to run the inversion downward from, and a starting
density at that altitude.  Those live here so that the lamination tier and
whatever replaces it cannot disagree about them.

There is no inversion physics in this module.  The lamination algebra is in
:mod:`tss.l4_retrieval.lamination`.

The starting-value problem, stated once
----------------------------------------
Topside inversion runs downward from the spacecraft, and every reflection
height it computes is measured from the previous one.  The first lamina is
therefore measured from the spacecraft, and it needs the plasma frequency at
the spacecraft.  There are three ways to get it and they give materially
different error budgets:

``truth``
    Taken from the L0 state.  Only available in the simulator.  Used
    deliberately, to isolate the inversion error from the starting-value error,
    which is the only way to attribute an error at L5 to the inversion rather
    than to its initial condition.
``in_situ``
    From a Langmuir probe or plasma frequency probe on the spacecraft.  This is
    what a real mission would carry and what Alouette and ISIS carried.  The
    baseline 6U configuration does not currently include one, which is a
    mission design question and not an inversion question, but it is the
    inversion that pays for the omission.
``lowest_trace``
    The lowest sounding frequency at which an echo was detected, used as a
    proxy for the local plasma frequency.  This is the operational fallback and
    it is biased: the lowest detected frequency is above the local plasma
    frequency by however much the sensitivity and the frequency step missed,
    so it understates the starting density and the whole profile inherits the
    error.

Which one was used is recorded in ``DensityProfile.flags["starting_value"]``.
``docs/architecture.md`` requires that, and the reason is that the three give
different error budgets and a retrieved profile that does not say which one it
used is not usable as evidence.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from tss.constants import density_from_plasma_frequency, plasma_frequency
from tss.interfaces import Ionogram, IonosphereState, PlatformState

__all__ = [
    "RetrievalBase",
    "STARTING_VALUE_SOURCES",
    "trace_points",
    "starting_plasma_frequency_hz",
    "density_from_trace",
]

#: How the plasma frequency at the spacecraft was obtained.  See the module
#: docstring; the three are not interchangeable.
STARTING_VALUE_SOURCES: tuple[str, ...] = ("truth", "in_situ", "lowest_trace")


def trace_points(
    ionogram: Ionogram,
    *,
    mode_index: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Detected trace from one mode column, sorted by increasing frequency.

    Returns ``(freq_hz, virtual_range_m, snr_db)`` with the undetected
    frequencies removed.  Undetected points are dropped rather than
    interpolated across, because an interpolated trace point would enter the
    lamination as though it were a measurement.

    ``mode_index`` selects the column.  The baseline L3 tier produces one
    column labelled ``"unresolved"``, because a scalar receiver cannot label
    O and X; see :mod:`tss.l3_signal.ionogram`.  A tier that does label them
    produces two, and which one is being inverted matters, because the X trace
    reflects at a different density for the same frequency and inverting it
    with the O-mode relation is wrong by the gyrofrequency shift.
    """
    freq = np.asarray(ionogram.freq_hz, dtype=float)
    vr = np.asarray(ionogram.virtual_range_m, dtype=float)
    snr = np.asarray(ionogram.snr_db, dtype=float)
    if vr.ndim != 2:
        raise ValueError(f"Ionogram.virtual_range_m must be 2-D; got {vr.shape}")
    if not 0 <= int(mode_index) < vr.shape[1]:
        raise ValueError(
            f"mode_index {mode_index} outside 0..{vr.shape[1] - 1}; "
            f"this ionogram carries modes {tuple(ionogram.modes)}"
        )
    column = vr[:, int(mode_index)]
    snr_column = snr[:, int(mode_index)] if snr.ndim == 2 else np.full_like(column, np.nan)

    keep = np.isfinite(freq) & np.isfinite(column)
    f = freq[keep]
    r = column[keep]
    s = snr_column[keep]
    order = np.argsort(f)
    return f[order], r[order], s[order]


def starting_plasma_frequency_hz(
    ionogram: Ionogram,
    *,
    given_hz: float | None = None,
    truth: IonosphereState | None = None,
    platform: PlatformState | None = None,
    mode_index: int = 0,
) -> tuple[float, str]:
    """Plasma frequency at the spacecraft, and where it came from.

    Resolution order: an explicitly given value, then the truth state if one is
    supplied, then the lowest detected sounding frequency.  The returned string
    is one of :data:`STARTING_VALUE_SOURCES` and belongs in
    ``DensityProfile.flags``.

    The ``lowest_trace`` fallback is biased high, which biases the starting
    density low.  It is the fallback and not the default choice for that
    reason.
    """
    if given_hz is not None:
        return float(given_hz), "in_situ"

    if truth is not None and platform is not None:
        alt = np.asarray(truth.alt_m, dtype=float)
        ne = np.asarray(truth.ne, dtype=float)
        ne_at_sc = float(np.interp(float(platform.alt_m), alt, ne))
        return float(plasma_frequency(ne_at_sc)), "truth"

    f, _, _ = trace_points(ionogram, mode_index=mode_index)
    if f.size == 0:
        raise ValueError(
            "no detected trace points and no starting value supplied, so the "
            "plasma frequency at the spacecraft cannot be established"
        )
    return float(f[0]), "lowest_trace"


def density_from_trace(freq_hz: np.ndarray) -> np.ndarray:
    """Reflection density for each sounding frequency, m^-3.

    ``Ne = (f / 8.98)^2`` in the usual units, written through
    :func:`tss.constants.density_from_plasma_frequency`.  This is the O-mode
    relation, which holds when the reflection condition is ``f = f_N``.  The X
    mode reflects where ``f = f_x``, which is above ``f_N`` by roughly half the
    gyrofrequency, so applying this to an X trace overstates the density.
    """
    return density_from_plasma_frequency(np.asarray(freq_hz, dtype=float))


class RetrievalBase:
    """Common base for L4 tiers.

    Holds the starting-value policy and the description machinery.  A subclass
    registers itself with :func:`tss.registry.register`, which sets ``name``
    and ``layer`` on the class, and implements ``invert`` to satisfy
    :class:`tss.interfaces.L4Retrieval`.
    """

    name: str = "base"
    layer: str = "l4"

    def __init__(
        self,
        *,
        mode_index: int = 0,
        start_plasma_hz: float | None = None,
    ) -> None:
        self.mode_index = int(mode_index)
        self.start_plasma_hz = None if start_plasma_hz is None else float(start_plasma_hz)

    def _settings(self) -> dict[str, Any]:
        return {
            "mode_index": self.mode_index,
            "start_plasma_hz": self.start_plasma_hz,
        }

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable identity and settings, for the run manifest."""
        record: dict[str, Any] = {
            "layer": "l4",
            "name": getattr(self, "name", type(self).__name__),
            "class": type(self).__name__,
            "module": type(self).__module__,
        }
        record.update(self._settings())
        return record
