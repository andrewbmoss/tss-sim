"""Unit conventions and conversions.

The convention is SI everywhere, with a suffix on the attribute name whenever a
quantity is carried in anything else.  A bare ``alt`` would be ambiguous, so the
codebase uses ``alt_m`` and ``alt_km`` and never ``alt``.  Suffixes in use:

    ``_m``  metres        ``_km``   kilometres
    ``_s``  seconds       ``_hz``   hertz        ``_mhz``  megahertz
    ``_k``  kelvin        ``_w``    watts        ``_dbm``  dBm
    ``_db`` decibels      ``_deg``  degrees      ``_rad``  radians
    ``_mps`` metres per second

The helpers here exist so that conversions are written once and can be tested,
rather than scattered as magic factors.  ``db`` and ``undb`` are power ratios,
which is the convention the link budget uses throughout.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "db",
    "undb",
    "db_amplitude",
    "km_to_m",
    "m_to_km",
    "mhz_to_hz",
    "hz_to_mhz",
    "dbm_to_w",
    "w_to_dbm",
    "deg_to_rad",
    "rad_to_deg",
]


def db(x):
    """Power ratio to decibels.  Non-positive input returns ``-inf``."""
    x = np.asarray(x, dtype=float)
    out = np.full(np.shape(x), -np.inf, dtype=float)
    np.log10(x, out=out, where=x > 0)
    return 10.0 * out if out.ndim else float(10.0 * out)


def undb(x_db):
    """Decibels to power ratio."""
    return 10.0 ** (np.asarray(x_db, dtype=float) / 10.0)


def db_amplitude(x):
    """Amplitude ratio to decibels, which is ``20*log10``."""
    x = np.asarray(x, dtype=float)
    out = np.full(np.shape(x), -np.inf, dtype=float)
    np.log10(np.abs(x), out=out, where=np.abs(x) > 0)
    return 20.0 * out if out.ndim else float(20.0 * out)


def km_to_m(x):
    return np.asarray(x, dtype=float) * 1.0e3


def m_to_km(x):
    return np.asarray(x, dtype=float) / 1.0e3


def mhz_to_hz(x):
    return np.asarray(x, dtype=float) * 1.0e6


def hz_to_mhz(x):
    return np.asarray(x, dtype=float) / 1.0e6


def dbm_to_w(x_dbm):
    return 10.0 ** (np.asarray(x_dbm, dtype=float) / 10.0) * 1.0e-3


def w_to_dbm(x_w):
    return db(np.asarray(x_w, dtype=float) / 1.0e-3)


def deg_to_rad(x):
    return np.deg2rad(np.asarray(x, dtype=float))


def rad_to_deg(x):
    return np.rad2deg(np.asarray(x, dtype=float))
