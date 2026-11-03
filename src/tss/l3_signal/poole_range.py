"""Layer: L3
WP: WP4 (P5 acceptance); WP5 (baseline comparator for C3)
Claims: C1, C3
Status: implemented
Requirements: INS-09, SCI-07

Poole's stationary-phase group range and phase velocity for FM-CW sounding.

Poole (1985, "Advanced sounding 1. The FMCW alternative", Radio Sci. 20(6),
1609-1616; "paper 1") showed that the phase of one discrete Fourier component
of the de-chirped receiver output behaves like the phase of a pulsed,
fixed-frequency echo: it is ``2 pi (2 f R*/c - 2k (R*/c)^2)`` (his eq. 14),
with ``R*`` the phase range.  Comparing that phase between cells:

* between two cells at the same start frequency, ``T`` apart, gives the phase
  (Doppler) velocity ``V*`` (eqs. 60, 72-73);
* between two cells whose start frequencies differ by ``df``, corrected for
  the Doppler change, gives the group range by stationary phase,
  ``R' = c (N + phi_R / 2 pi) / (2 df)`` (eqs. 47-49, 70-71), with the integer
  ``N`` resolved from the beat frequency of the selected component.

That separates group range from velocity with one antenna, which a single
sweep's beat frequency cannot do (``f_D = 2kR'/c + 2 f_c V*/c``, eq. 23).  The
resolution of ``N`` holds while ``f_dopp + f' < f_R1 / 2`` (eq. 59), with
``f_R1 = k / df`` (eq. 51) and ``f'`` the offset of the echo from its bin.

The complex components are taken with the same taper and the same scaling as
:func:`tss.l3_signal.transform.transform`, so a phase read here belongs to the
cell the chain's detector sees.  Poole's three-cell structure uses two
receivers and four antennas for arrival angle and O/X; with one antenna, as
here, the cells give ``R'``, ``V*`` and the amplitude.

WP5 uses :func:`three_cell` as the baseline comparator for C3: it is the
published single-antenna way to separate range from radial velocity.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np

from tss.l3_signal.transform import window_gains, window_taper

__all__ = [
    "C_LIGHT",
    "ambiguity_interval_m",
    "range_frequency_hz",
    "max_unambiguous_doppler_hz",
    "max_unambiguous_velocity_mps",
    "complex_bins",
    "group_range_m",
    "doppler_velocity_mps",
    "three_cell",
    "phase_error_rms_rad",
]

C_LIGHT = 299_792_458.0


def ambiguity_interval_m(delta_f_hz: float) -> float:
    """Group-range interval over which the cell-3 phase wraps once, ``c / (2 df)`` (eq. 48)."""
    return C_LIGHT / (2.0 * float(delta_f_hz))


def range_frequency_hz(k_hz_per_s: float, delta_f_hz: float) -> float:
    """``f_R1 = k / df`` (eq. 51): the beat-frequency width of one ambiguity interval."""
    return float(k_hz_per_s) / float(delta_f_hz)


def max_unambiguous_doppler_hz(k_hz_per_s: float, delta_f_hz: float, f_prime_hz: float = 0.0) -> float:
    """Largest Doppler frequency for which ``N`` is resolved, ``f_R1 / 2 - |f'|`` (eq. 59)."""
    return 0.5 * range_frequency_hz(k_hz_per_s, delta_f_hz) - abs(float(f_prime_hz))


def max_unambiguous_velocity_mps(k_hz_per_s: float, delta_f_hz: float, f_c_hz: float,
                                 f_prime_hz: float = 0.0) -> float:
    """The same bound as a phase velocity, through ``f_dopp = 2 f_c V* / c`` (eq. 45)."""
    return C_LIGHT * max_unambiguous_doppler_hz(k_hz_per_s, delta_f_hz, f_prime_hz) / (2.0 * float(f_c_hz))


def complex_bins(samples: np.ndarray, fs_hz: float, *, window: str = "hann",
                 n_fft: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Beat frequencies and complex amplitudes of one de-chirped cell.

    Same taper and scaling as the chain's transform: a unit tone filling the
    cell reads unit magnitude.  Returns the non-negative half.
    """
    x = np.asarray(samples)
    x = x[0] if x.ndim == 2 else x
    n = int(x.size)
    length = n if n_fft is None else int(n_fft)
    taper = window_taper(window, n)
    cg, _ = window_gains(taper)
    full = np.fft.fft(x * taper, n=length)
    keep = length // 2 + 1
    return np.arange(keep) * float(fs_hz) / length, full[:keep] / (n * cg)


def _wrap(phi: float) -> float:
    return float((phi + math.pi) % (2.0 * math.pi) - math.pi)


def group_range_m(phi_r_rad: float, f_i_hz: float, k_hz_per_s: float, delta_f_hz: float) -> tuple[float, int]:
    """Group range by stationary phase: ``R' = c (N + phi_R/2pi) / (2 df)``, eqs. 70-71.

    ``N`` is the nearest integer to ``f_i df / k - phi_R / 2 pi``.  Returns
    ``(R', N)``.
    """
    phi = _wrap(float(phi_r_rad))
    n = int(round(float(f_i_hz) * float(delta_f_hz) / float(k_hz_per_s) - phi / (2.0 * math.pi)))
    return C_LIGHT * (n + phi / (2.0 * math.pi)) / (2.0 * float(delta_f_hz)), n


def doppler_velocity_mps(phi_d_rad: float, f_i_hz: float, f_r_hz: float, cell_s: float,
                         f_ref_hz: float) -> tuple[float, int]:
    """Phase velocity from the Doppler phase between two cells ``T`` apart, eqs. 72-73.

    ``V* = c (M + phi_D/2pi) / (2 T f_ref)``, ``M`` the nearest integer to
    ``T (f_i - f_R) - phi_D / 2 pi``.  ``f_ref`` is Poole's ``f_c``, the cell's
    mean frequency: the phase of a component of a tapered transform is that of
    the taper's centroid, the middle of the cell, where the instantaneous
    frequency is ``f_c``.  (Referring it to ``f_0`` is 0.4 % wrong at Poole
    and Evans's settings, checked on synthetic cells.)
    """
    phi = _wrap(float(phi_d_rad))
    m = int(round(float(cell_s) * (float(f_i_hz) - float(f_r_hz)) - phi / (2.0 * math.pi)))
    return C_LIGHT * (m + phi / (2.0 * math.pi)) / (2.0 * float(cell_s) * float(f_ref_hz)), m


def three_cell(cells: Sequence[np.ndarray], fs_hz: float, *, k_hz_per_s: float, f0_hz: float,
               delta_f_hz: float, window: str = "hann", n_fft: int | None = None,
               bin_index: int | None = None, f_ref: str = "fc") -> dict[str, Any]:
    """Poole's single-antenna three-cell analysis.

    ``cells`` are the de-chirped samples of cell 1 and cell 2 (start frequency
    ``f0``, ``T`` apart) and cell 3 (start frequency ``f0 + df``, ``T`` after
    cell 2).  The component is the strongest of cell 1 unless ``bin_index`` is
    given (Poole and Evans pick it in cell 1).  ``phi_D = phi_2 - phi_1``
    (eq. 34) and ``phi_R = (phi_3 - phi_2) - phi_D`` (eq. 42).  Poole's
    algebra drops a term of order ``V* T`` in ``R'`` (the target moves while
    the cells run): evaluated at each cell's centroid, ``phi_R`` carries
    ``5 df V* T / c`` beyond eq. 42, so ``R'`` reads ``2.5 V* T`` long (-0.25 km
    at -200 m/s, +0.31 km at +250 m/s, at T = 0.5 s), and eq. 59's velocity
    bound moves by about 2 % (``fmcw_acceptance`` A5).
    """
    if len(cells) != 3:
        raise ValueError("three cells are needed: two at f0, one at f0 + df")
    spectra = [complex_bins(c, fs_hz, window=window, n_fft=n_fft) for c in cells]
    freqs = spectra[0][0]
    i = int(np.argmax(np.abs(spectra[0][1]))) if bin_index is None else int(bin_index)
    ph = [float(np.angle(s[1][i])) for s in spectra]
    n = int(np.asarray(cells[0]).shape[-1])
    cell_s = n / float(fs_hz)
    phi_d = _wrap(ph[1] - ph[0])
    phi_r = _wrap(ph[2] - ph[1] - phi_d)
    f_i = float(freqs[i])
    r_g, n_amb = group_range_m(phi_r, f_i, k_hz_per_s, delta_f_hz)
    f_r = 2.0 * k_hz_per_s * r_g / C_LIGHT
    f_c = float(f0_hz) + 0.5 * k_hz_per_s * cell_s
    ref = float(f0_hz) if f_ref == "f0" else f_c
    v, m = doppler_velocity_mps(phi_d, f_i, f_r, cell_s, ref)
    return {"group_range_m": r_g, "phase_velocity_mps": v, "n": n_amb, "m": m, "phi_r_rad": phi_r,
            "phi_d_rad": phi_d, "f_i_hz": f_i, "bin_index": i, "amplitude": float(np.abs(spectra[0][1][i])),
            "f_c_hz": f_c, "cell_s": cell_s}


def phase_error_rms_rad(snr_db: float, n: int = 20001) -> float:
    """Poole's phase error of a Fourier component with a noise phasor at random phase (eqs. 80-82).

    ``X = N/S`` (voltage), ``theta = atan(X sin e / (1 + X cos e))``, and the rms
    over ``e`` uniform in ``(-pi, pi)``.  20 dB gives 4.05 degrees.
    """
    x = 10.0 ** (-float(snr_db) / 20.0)
    e = np.linspace(-math.pi, math.pi, int(n))
    th = np.arctan2(x * np.sin(e), 1.0 + x * np.cos(e))
    return float(np.sqrt(np.mean(th**2)))
