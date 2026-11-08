"""Layer: support
WP: WP3, WP4
Claims: C1
Status: implemented
Requirements: -

Window tapers and their gains, shared by L2 (the reference echo of WP3) and
L3 (the range transform).

They sit below the layers so that neither imports the other for them.  Until
2026-11-08 they lived in :mod:`tss.l3_signal.transform`, which still
re-exports them, and L2 imported them from there, an import from L2 into L3
(code review F4).
"""

from __future__ import annotations

import numpy as np

__all__ = ["WINDOWS", "window_taper", "window_gains"]

#: Window tapers implemented here.  They are written out in numpy rather than
#: taken from ``scipy.signal`` because this module may not import scipy at
#: module scope and a lazy import for four cosine coefficients is not worth the
#: indirection.  The periodic (``n / N``) convention is used throughout, which
#: is the correct one for spectral analysis; the symmetric (``n / (N-1)``)
#: convention that :func:`numpy.hanning` uses belongs to filter design and
#: leaks a little more.
WINDOWS: tuple[str, ...] = ("rect", "hann", "hamming", "blackman", "blackmanharris")

_COSINE_COEFFS: dict[str, tuple[float, ...]] = {
    "rect": (1.0,),
    "hann": (0.5, -0.5),
    "hamming": (0.54, -0.46),
    "blackman": (0.42, -0.5, 0.08),
    "blackmanharris": (0.35875, -0.48829, 0.14128, -0.01168),
}


def window_taper(name: str, n: int) -> np.ndarray:
    """Return the named taper of length ``n``, periodic convention.

    Raises
    ------
    ValueError
        For an unknown name.  A silent fall back to a rectangular window would
        change the sidelobe floor by tens of decibels without saying so, which
        is exactly the sort of quiet substitution that makes a dynamic-range
        result untrustworthy.
    """
    key = str(name).lower()
    if key not in _COSINE_COEFFS:
        raise ValueError(f"unknown window {name!r}; expected one of {WINDOWS}")
    if int(n) < 1:
        raise ValueError(f"window length must be at least 1; got {n}")
    n = int(n)
    coeffs = _COSINE_COEFFS[key]
    k = np.arange(n, dtype=float)
    taper = np.zeros(n, dtype=float)
    for order, a in enumerate(coeffs):
        taper += a * np.cos(2.0 * np.pi * order * k / n)
    return taper


def window_gains(taper: np.ndarray) -> tuple[float, float]:
    """Coherent gain and equivalent noise bandwidth in bins for a taper.

    Coherent gain is ``sum(w) / n``: what the taper does to the amplitude of a
    tone sitting exactly on a bin.  Equivalent noise bandwidth is
    ``n sum(w^2) / (sum w)^2``: how many bins' worth of white noise one bin
    collects.  For a rectangular window both are 1.0 and 1.0; for Hann they are
    0.5 and 1.5.

    They are returned together because using one without the other is how a
    signal-to-noise figure acquires a factor of the window gain.
    """
    w = np.asarray(taper, dtype=float)
    n = w.size
    total = float(np.sum(w))
    if total == 0.0:
        raise ValueError("taper sums to zero, so its gains are undefined")
    coherent = total / n
    enbw = n * float(np.sum(w**2)) / (total**2)
    return coherent, enbw
