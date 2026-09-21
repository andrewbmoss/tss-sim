"""Physical constants and the two relations the requirement set is built on.

CODATA 2018 values.  These are duplicated nowhere else in the package; a module
that needs a constant imports it from here so that there is one place to check.

The two relations of Section 3.2 of the proposal are implemented here rather
than in a layer, because both L0 and L5 need them and neither owns them.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "C_LIGHT",
    "K_BOLTZ",
    "Q_E",
    "M_E",
    "EPS0",
    "MU0",
    "R_EARTH_M",
    "PLASMA_CONST",
    "T0_K",
    "plasma_frequency",
    "density_from_plasma_frequency",
    "density_precision_from_freq_step",
    "x_mode_cutoff",
    "gyrofrequency",
    "beat_frequency",
    "range_cell_noise_bandwidth",
    "doppler_range_bias",
]

C_LIGHT = 299_792_458.0          # m/s, exact
K_BOLTZ = 1.380_649e-23          # J/K, exact
Q_E = 1.602_176_634e-19          # C, exact
M_E = 9.109_383_7015e-31         # kg
EPS0 = 8.854_187_8128e-12        # F/m
MU0 = 1.256_637_062_12e-6        # H/m
R_EARTH_M = 6_371_000.0          # mean radius, m
T0_K = 290.0                     # reference temperature for noise figure, K

#: ``f_N [Hz] = PLASMA_CONST * sqrt(Ne [m^-3])``, the 8.98 of the proposal.
PLASMA_CONST = 1.0 / (2.0 * math.pi) * math.sqrt(Q_E**2 / (EPS0 * M_E))


def plasma_frequency(ne):
    """Electron plasma frequency in Hz from density in m^-3."""
    return PLASMA_CONST * np.sqrt(np.asarray(ne, dtype=float))


def density_from_plasma_frequency(f_hz):
    """Electron density in m^-3 from plasma frequency in Hz."""
    return (np.asarray(f_hz, dtype=float) / PLASMA_CONST) ** 2


def density_precision_from_freq_step(delta_f, f):
    """Fractional density precision delivered by a fractional frequency step.

    Because ``f_N`` goes as the square root of density, a fractional frequency
    step delivers twice the fractional density step:

        ``dNe/Ne = 2 * df/f``

    This is the relation that makes density precision and height resolution one
    requirement in two variables (INS-02, SCI-02).
    """
    return 2.0 * np.asarray(delta_f, dtype=float) / np.asarray(f, dtype=float)


def gyrofrequency(b_tesla):
    """Electron gyrofrequency in Hz from flux density in tesla."""
    return Q_E * np.asarray(b_tesla, dtype=float) / (2.0 * math.pi * M_E)


def x_mode_cutoff(f_p_hz, f_ce_hz):
    """Extraordinary-mode cutoff frequency.

    ``f_x = f_ce/2 + sqrt((f_ce/2)^2 + f_p^2)``.  The X trace therefore always
    lies above the O trace, which is what makes polarization separation a
    workable way to label the two (INS-08, C5).
    """
    f_p_hz = np.asarray(f_p_hz, dtype=float)
    half_ce = np.asarray(f_ce_hz, dtype=float) / 2.0
    return half_ce + np.sqrt(half_ce**2 + f_p_hz**2)


def beat_frequency(range_m, k_f_hz_per_s):
    """De-chirped beat frequency for a two-way range.

    ``f_b = 2 * K_f * R / c``.  This is the relation that turns range into a
    Fourier bin index and is the reason an FMCW sounder can work at low power.
    """
    return 2.0 * np.asarray(k_f_hz_per_s, dtype=float) * np.asarray(range_m, dtype=float) / C_LIGHT


def range_cell_noise_bandwidth(range_resolution_m, k_f_hz_per_s):
    """Noise bandwidth admitted by one range cell, Hz.

    ``B = 2 * K_f * dR / c``.  For a 1 km cell at K_f = 5.6 MHz/s this is about
    37 Hz, against the tens of kilohertz a pulse sounder admits.  The whole
    sub-watt argument of C1 comes from this number.
    """
    return 2.0 * np.asarray(k_f_hz_per_s, dtype=float) * np.asarray(range_resolution_m, dtype=float) / C_LIGHT


def doppler_range_bias(velocity_mps, f_hz, k_f_hz_per_s):
    """Apparent-range bias produced by a Doppler shift in a linear sweep.

    ``dd = v * f / K_f``.  It is independent of range and grows linearly with
    sounding frequency, so it is worst at the top of the band.  Removing it
    from a single transmit and receive pair is claim C3, the principal
    contribution proposed.
    """
    return (
        np.asarray(velocity_mps, dtype=float)
        * np.asarray(f_hz, dtype=float)
        / np.asarray(k_f_hz_per_s, dtype=float)
    )
