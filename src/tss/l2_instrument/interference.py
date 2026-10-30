"""Layer: L2
WP: WP4
Claims: C1, C4
Status: implemented
Requirements: INS-10

Narrowband carriers at the receiver input (WP4 Phase 1, broadcast interference).

A carrier at ``f_c`` with power ``P`` at the receiver input, after the de-chirp
against the reference (the receiver's convention: baseband = reference times
the conjugate of the received signal, so an echo of delay ``tau`` beats at
``+K tau`` on an up-sweep), is

    x(t) = sqrt(P) exp( j [ phi_ref(t) - 2 pi f_c t - theta ] )

whose instantaneous frequency is ``f_ref(t) - f_c``.  It sweeps through the
de-chirped band as the reference passes the carrier, which is why a carrier
shows on an ionogram as a line across every range at its frequency.  The
digitiser's anti-alias filter passes ``|f_ref - f_c| < fs/2`` only; the carrier
is rendered there, with a raised-cosine edge over the outer ``taper`` fraction
of the band, and is zero elsewhere.

The power at the receiver input is the carrier's isotropic power plus the
antenna gain at ``f_c`` (which carries the radiation efficiency) plus the
instrument's polarization loss.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any

import numpy as np

__all__ = ["carrier_baseband", "render_carriers"]


def carrier_baseband(phi_ref_rad: np.ndarray, f_ref_hz: np.ndarray, t_s: np.ndarray, f_c_hz: float,
                     power_w: float, phase_rad: float, fs_hz: float, taper: float = 0.1) -> np.ndarray:
    """One carrier, de-chirped, on the samples given (module docstring)."""
    fb = np.asarray(f_ref_hz, dtype=float) - float(f_c_hz)
    half = 0.5 * float(fs_hz)
    a = np.abs(fb) / half
    edge = 1.0 - float(taper)
    w = np.where(a <= edge, 1.0,
                 np.where(a < 1.0, 0.5 * (1.0 + np.cos(np.pi * (a - edge) / max(taper, 1e-12))), 0.0))
    out = np.zeros(np.shape(t_s), dtype=complex)
    live = w > 0.0
    if np.any(live):
        ph = np.asarray(phi_ref_rad, dtype=float)[live] - 2.0 * np.pi * float(f_c_hz) * np.asarray(t_s)[live] \
            - float(phase_rad)
        out[live] = np.sqrt(float(power_w)) * w[live] * np.exp(1j * ph)
    return out


def render_carriers(carriers: Mapping[str, Any] | None, *, phi_ref_rad: np.ndarray, f_ref_hz: np.ndarray,
                    t_s: np.ndarray, fs_hz: float, gain_dbi: Callable[[float], float],
                    polarization_loss_db: float = 0.0) -> tuple[np.ndarray, dict[str, Any]]:
    """Sum of every carrier in ``carriers`` (an ``EchoSet.interference`` mapping)."""
    out = np.zeros(np.shape(t_s), dtype=complex)
    info: dict[str, Any] = {"carriers_listed": 0, "carriers_rendered": 0, "carrier_power_dbm_max": None}
    if not carriers:
        return out, info
    f_c = np.asarray(carriers.get("freq_hz", []), dtype=float)
    p_iso = np.asarray(carriers.get("power_dbm_iso", []), dtype=float)
    ph = np.asarray(carriers.get("phase_rad", np.zeros(f_c.size)), dtype=float)
    info["carriers_listed"] = int(f_c.size)
    f_lo, f_hi = float(np.min(f_ref_hz)), float(np.max(f_ref_hz))
    p_max = None
    for f, p, th in zip(f_c, p_iso, ph, strict=True):
        if not (np.isfinite(f) and np.isfinite(p)) or f < f_lo - fs_hz or f > f_hi + fs_hz:
            continue
        p_dbm = float(p) + float(gain_dbi(f / 1.0e6)) + float(polarization_loss_db)
        x = carrier_baseband(phi_ref_rad, f_ref_hz, t_s, f, 10.0 ** ((p_dbm - 30.0) / 10.0), th, fs_hz)
        if np.any(x != 0.0):
            out += x
            info["carriers_rendered"] += 1
            p_max = p_dbm if p_max is None else max(p_max, p_dbm)
    info["carrier_power_dbm_max"] = p_max
    return out, info
