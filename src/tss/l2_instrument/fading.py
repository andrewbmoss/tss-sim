"""Layer: L2
WP: WP5
Claims: C3
Status: partial
Requirements: SCI-15, INS-13

Echo fading: the complex gain an echo is multiplied by, sweep to sweep and
frequency to frequency.

STATUS
------
Standalone (WP5 Phase 1, 2026-11-24).  Nothing in the L2 tier or the
route-one estimator calls it yet; Phase 2 wires it in.  ``Status: partial``
until then.

WHAT THE LITERATURE GIVES, AND WHAT IT DOES NOT
-----------------------------------------------
No measurement of echo decorrelation from a topside sounder in orbit was found
(search of 2026-11-24, docs/decorrelation.md).  Every number comes from a
stationary receiver on the ground, mostly looking at the bottomside: ITU-R
F.1487 Doppler spreads (field 1/e times of 4.5, 0.9 and 0.45 s for quiet,
moderate and disturbed mid-latitude paths, down to 15 to 45 ms at low and high
latitudes when disturbed), Digisonde practice (coherent integration of 0.16 to
2.56 s, coherence for tens of seconds "not rare") and Ganguly et al. (2001)
(tens of ms to several s).  Those numbers bound only the *intrinsic* change of
the irregularities.  Gate G1 (2026-11-24) takes them as the starting point,
held at low confidence and high variance: the intrinsic time of each echo is
drawn about the carried value with a log-normal spread,
:data:`LOW_CONFIDENCE_LOG_SIGMA`, and the carried values are
:data:`CARRY_TAU_I_S`.

THE MODEL
---------
The gain is a Rician sum of two unit-power parts,

    g = sqrt(K / (K+1)) s + sqrt(1 / (K+1)) d,

so its mean power is one for every ``K``.

``s``, the specular part, is the facet's own reflection.  With no intrinsic
time it is the constant 1, and ``K -> inf`` with no intrinsic time is no
fading at all.  With an intrinsic time ``tau_i`` it is a complex Gaussian
process of that correlation time, so over times short against ``tau_i`` the
echo is Rician and over long times the specular part wanders too.  ``K -> inf``
with ``tau_i = 0.5 s`` is the Rayleigh model the route-one budget has used
since October.

``d``, the diffuse part, is the randomly phased scatter around the facet.  Its
decorrelation is tied to the facet (G1 decision 2).  A randomly phased patch of
along-track extent ``L`` at apparent range ``P'`` has a Doppler spectrum
uniform over ``B = 2 v_s L / (lambda P')`` in the convention of route one, so
its field correlation is ``sinc(pi B tau)`` and falls to 1/e at

    tau_d = x_e lambda P' / (2 pi v_s L),     x_e = 2.199,

the root of ``sin x / x = 1/e`` (:func:`halo_decorrelation_s`).  The same
geometry makes the power envelope of a coherent facet of size ``D`` have the
half-power width ``W = x_h lambda P' / (pi D v_s)`` with ``x_h = 1.392``,
so

    tau_d / W = (x_e / 2 x_h) D / L = 0.79 D / L,

whatever the frequency and range (:func:`diffuse_to_envelope_ratio`).  A
diffuse halo wider than about 0.8 D decorrelates faster than the facet's
envelope.  This is the van Cittert-Zernike duality: motion past a randomly
phased patch and motion past a coherent facet of the same size take the same
time.  It is why the 0.1 s case that the route-one budget carried until
2026-11-24 is an extended reflector (a halo of about 14 km at 5 MHz and
500 km) and not a compact facet.

The ``"exp"`` kind gives both parts the correlation ``exp(-|tau| / tau_c)``,
realised exactly on any increasing time grid as a first-order autoregression.
The ``"gauss"`` kind gives ``exp(-tau^2 / tau_c^2)``, realised on a uniform
grid by filtering white noise with a Gaussian kernel of standard deviation
``tau_c / 2``.  Either way ``tau_c`` is the field 1/e time.  A part whose time
is shorter than the grid step, or is ``None`` for the diffuse part, is white.

ACROSS FREQUENCY
----------------
Each part may carry a coherence bandwidth ``B_c``.  Its correlation between
two frequencies is then ``exp(-|df| / B_c)``, realised as a first-order
autoregression across the sorted frequencies of independent time processes of
the same kind, which keeps the time correlation and makes the two separable.
With no ``B_c`` one draw serves every frequency.  The halo gives the diffuse
part a bandwidth from the spread of its excess delays: for a uniform patch the
two-way excess delay ``x^2 / (c P')`` has standard deviation
``L^2 / (6 sqrt 5 c P')``, and ``B_c = 1 / (2 pi sigma_tau)``
(:func:`halo_coherence_bandwidth_hz`).  That is an order-of-magnitude figure;
nothing has measured it.

RANDOM STREAMS
--------------
Each echo gets its own generator (:func:`echo_streams`, G1 decision of
2026-11-24 "per-echo streams").  Inside an echo the intrinsic-time draw, the
specular part and the diffuse part each take a child spawned from that
generator, so switching the diffuse part on or off leaves the specular draw as
it was (code review F3).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.optimize import brentq
from scipy.signal import lfilter

from tss.constants import C_LIGHT

__all__ = [
    "CARRY_TAU_I_S",
    "CARRY_K",
    "CARRY_HALO_RATIOS",
    "LOW_CONFIDENCE_LOG_SIGMA",
    "DEFAULT_V_S_MPS",
    "FadingModel",
    "FadingDraw",
    "echo_streams",
    "halo_decorrelation_s",
    "halo_extent_m",
    "diffuse_to_envelope_ratio",
    "halo_coherence_bandwidth_hz",
]

#: Intrinsic (specular) decorrelation times carried, s.  G1, 2026-11-24.
CARRY_TAU_I_S = (0.1, 0.5, 1.0)

#: Rician K factors carried: no diffuse part, then 10, 3, 1 and Rayleigh.
CARRY_K = (math.inf, 10.0, 3.0, 1.0, 0.0)

#: Diffuse halo extent as a multiple of the facet size.
CARRY_HALO_RATIOS = (1.0, 2.0, 4.0)

#: Standard deviation of ``ln tau_i`` across echoes: a factor of 3 at one
#: sigma.  The literature spans 15 ms to 4.5 s, all of it from the ground, so
#: G1 holds it at low confidence and high variance.
LOW_CONFIDENCE_LOG_SIGMA = math.log(3.0)

#: Orbital speed used by route one, m/s.
DEFAULT_V_S_MPS = 7600.0

_KINDS = ("exp", "gauss")

#: Root of sin(x)/x = 1/e: the 1/e point of a uniform patch's field correlation.
_X_E = float(brentq(lambda x: math.sin(x) / x - math.exp(-1.0), 1.5, 3.0, xtol=1e-14))

#: Root of sin(x)/x = 1/sqrt(2): the half-power point of a facet's envelope.
_X_H = float(brentq(lambda x: math.sin(x) / x - math.sqrt(0.5), 1.0, 2.0, xtol=1e-14))

#: Kernel half-length in kernel standard deviations for the Gaussian kind.
_KERNEL_SIGMAS = 4.0


# ---------------------------------------------------------------------------
# The halo relations
# ---------------------------------------------------------------------------


def halo_decorrelation_s(
    f_hz: float, p_group_m: float, halo_m: float, v_s_mps: float = DEFAULT_V_S_MPS
) -> float:
    """Field 1/e time of a randomly phased patch of extent ``halo_m``, s."""
    lam = C_LIGHT / float(f_hz)
    return _X_E * lam * float(p_group_m) / (2.0 * math.pi * float(v_s_mps) * float(halo_m))


def halo_extent_m(
    tau_s: float, f_hz: float, p_group_m: float, v_s_mps: float = DEFAULT_V_S_MPS
) -> float:
    """The patch extent whose decorrelation time is ``tau_s``, m."""
    lam = C_LIGHT / float(f_hz)
    return _X_E * lam * float(p_group_m) / (2.0 * math.pi * float(v_s_mps) * float(tau_s))


def diffuse_to_envelope_ratio(facet_m: float, halo_m: float) -> float:
    """``tau_d / W`` for a halo of ``halo_m`` around a facet of ``facet_m``."""
    return (_X_E / (2.0 * _X_H)) * float(facet_m) / float(halo_m)


def halo_coherence_bandwidth_hz(p_group_m: float, halo_m: float) -> float:
    """Coherence bandwidth of a uniform patch, ``1 / (2 pi sigma_tau)``, Hz."""
    sigma_tau = float(halo_m) ** 2 / (6.0 * math.sqrt(5.0) * C_LIGHT * float(p_group_m))
    return 1.0 / (2.0 * math.pi * sigma_tau)


# ---------------------------------------------------------------------------
# Streams and processes
# ---------------------------------------------------------------------------


def echo_streams(rng: np.random.Generator, n: int) -> list[np.random.Generator]:
    """``n`` independent child generators of ``rng``, one per echo.

    A generator that cannot spawn is returned ``n`` times, as in the L2 tier.
    """
    try:
        return list(rng.spawn(int(n)))
    except (AttributeError, TypeError, ValueError):
        return [rng] * int(n)


def _white(rng: np.random.Generator, n: int) -> np.ndarray:
    return (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / math.sqrt(2.0)


def _uniform_step(t: np.ndarray) -> float | None:
    if t.size < 2:
        return None
    d = np.diff(t)
    if np.allclose(d, d[0], rtol=1e-9, atol=0.0):
        return float(d[0])
    return None


def _process(t: np.ndarray, tau: float | None, kind: str, rng: np.random.Generator) -> np.ndarray:
    """Unit-power complex Gaussian process on ``t`` with field 1/e time ``tau``."""
    n = t.size
    if n == 0:
        return np.zeros(0, dtype=complex)
    d = np.diff(t)
    if tau is None or (d.size and float(np.min(d)) > float(tau)):
        return _white(rng, n)
    if kind == "exp":
        w = _white(rng, n)
        step = _uniform_step(t)
        if step is not None:
            r = math.exp(-step / float(tau))
            innov = math.sqrt(1.0 - r * r)
            x = w.copy()
            x[0] = w[0] / innov
            return lfilter([innov], [1.0, -r], x)
        g = np.empty(n, dtype=complex)
        g[0] = w[0]
        for k in range(1, n):
            r = math.exp(-d[k - 1] / float(tau))
            g[k] = r * g[k - 1] + math.sqrt(1.0 - r * r) * w[k]
        return g
    # Gaussian: needs a uniform grid.
    step = _uniform_step(t)
    if step is None:
        if n == 1:
            return _white(rng, 1)
        raise ValueError("the gauss kind needs a uniform time grid")
    s = 0.5 * float(tau) / step  # kernel standard deviation in samples
    m = int(math.ceil(_KERNEL_SIGMAS * s))
    k = np.arange(-m, m + 1, dtype=float)
    h = np.exp(-0.5 * (k / s) ** 2)
    h /= math.sqrt(float(np.sum(h * h)))
    w = _white(rng, n + 2 * m)
    return np.convolve(w, h, mode="valid")


def _over_frequency(
    t: np.ndarray,
    tau: float | None,
    kind: str,
    f_hz: np.ndarray | None,
    bw_hz: float | None,
    rng: np.random.Generator,
) -> np.ndarray:
    """One process, or one per frequency correlated ``exp(-|df| / B_c)``."""
    if f_hz is None:
        return _process(t, tau, kind, rng)
    f = np.asarray(f_hz, dtype=float).ravel()
    if bw_hz is None:
        one = _process(t, tau, kind, rng)
        return np.tile(one, (f.size, 1))
    order = np.argsort(f, kind="stable")
    fs = f[order]
    out_sorted = np.empty((f.size, t.size), dtype=complex)
    prev = None
    for j in range(f.size):
        w = _process(t, tau, kind, rng)
        if prev is None:
            cur = w
        else:
            r = math.exp(-abs(fs[j] - fs[j - 1]) / float(bw_hz))
            cur = r * prev + math.sqrt(1.0 - r * r) * w
        out_sorted[j] = cur
        prev = cur
    out = np.empty_like(out_sorted)
    out[order] = out_sorted
    return out


def _rho(lag: np.ndarray, tau: float | None, kind: str, *, white_if_none: bool) -> np.ndarray:
    lag = np.abs(np.asarray(lag, dtype=float))
    if tau is None:
        return (lag == 0.0).astype(float) if white_if_none else np.ones_like(lag)
    if kind == "exp":
        return np.exp(-lag / float(tau))
    return np.exp(-((lag / float(tau)) ** 2))


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FadingDraw:
    """One echo's realisation: the gain, its two parts and the times used."""

    gain: np.ndarray
    specular: np.ndarray
    diffuse: np.ndarray
    tau_i_s: float | None
    tau_d_s: float | None


@dataclass(frozen=True)
class FadingModel:
    """Rician two-time-scale fading of one echo.  See the module docstring.

    ``tau_i_s``: intrinsic time of the specular part, s; ``None`` is steady.
    ``k_factor``: specular to diffuse power; ``inf`` has no diffuse part and
    0 is Rayleigh.  ``tau_d_s``: diffuse time, s; ``None`` is white.
    ``kind``: ``"exp"`` or ``"gauss"``.  ``coherence_bw_hz`` and
    ``specular_coherence_bw_hz``: 1/e bandwidths across frequency of the
    diffuse and specular parts; ``None`` is one draw for all frequencies.
    ``tau_log_sigma``: spread of ``ln tau_i`` across echoes.
    """

    tau_i_s: float | None = None
    k_factor: float = math.inf
    tau_d_s: float | None = None
    kind: str = "exp"
    coherence_bw_hz: float | None = None
    specular_coherence_bw_hz: float | None = None
    tau_log_sigma: float = LOW_CONFIDENCE_LOG_SIGMA
    halo_m: float | None = field(default=None, compare=False)

    def __post_init__(self) -> None:
        for name in ("tau_i_s", "tau_d_s", "coherence_bw_hz", "specular_coherence_bw_hz"):
            v = getattr(self, name)
            if v is not None and not (float(v) > 0.0):
                raise ValueError(f"{name} must be positive or None, got {v!r}")
        if not (float(self.k_factor) >= 0.0):
            raise ValueError(f"k_factor must be >= 0, got {self.k_factor!r}")
        if self.kind not in _KINDS:
            raise ValueError(f"kind must be one of {_KINDS}, got {self.kind!r}")
        if not (float(self.tau_log_sigma) >= 0.0):
            raise ValueError(f"tau_log_sigma must be >= 0, got {self.tau_log_sigma!r}")

    @classmethod
    def from_halo(
        cls,
        *,
        f_hz: float,
        p_group_m: float,
        facet_m: float,
        halo_ratio: float,
        tau_i_s: float | None = 0.5,
        k_factor: float = 3.0,
        kind: str = "exp",
        v_s_mps: float = DEFAULT_V_S_MPS,
        tau_log_sigma: float = LOW_CONFIDENCE_LOG_SIGMA,
        specular_coherence_bw_hz: float | None = None,
    ) -> FadingModel:
        """A model whose diffuse part is a halo of ``halo_ratio`` facets."""
        halo = float(halo_ratio) * float(facet_m)
        return cls(
            tau_i_s=tau_i_s,
            k_factor=k_factor,
            tau_d_s=halo_decorrelation_s(f_hz, p_group_m, halo, v_s_mps),
            kind=kind,
            coherence_bw_hz=halo_coherence_bandwidth_hz(p_group_m, halo),
            specular_coherence_bw_hz=specular_coherence_bw_hz,
            tau_log_sigma=tau_log_sigma,
            halo_m=halo,
        )

    # -- weights and closed forms ------------------------------------------

    def _weights(self) -> tuple[float, float]:
        k = float(self.k_factor)
        if math.isinf(k):
            return 1.0, 0.0
        return math.sqrt(k / (k + 1.0)), math.sqrt(1.0 / (k + 1.0))

    def acf(self, lag_s: Any, df_hz: Any = 0.0) -> np.ndarray | float:
        """Normalised field correlation at a time lag and frequency offset.

        At the carried ``tau_i``, without the per-echo spread.
        """
        a, b = self._weights()
        lag = np.asarray(lag_s, dtype=float)
        df = np.abs(np.asarray(df_hz, dtype=float))
        rs = _rho(lag, self.tau_i_s, self.kind, white_if_none=False)
        rd = _rho(lag, self.tau_d_s, self.kind, white_if_none=True)
        fs = np.ones_like(df) if self.specular_coherence_bw_hz is None else np.exp(
            -df / float(self.specular_coherence_bw_hz)
        )
        fdiff = np.ones_like(df) if self.coherence_bw_hz is None else np.exp(
            -df / float(self.coherence_bw_hz)
        )
        out = a * a * rs * fs + b * b * rd * fdiff
        return float(out) if np.ndim(out) == 0 else out

    # -- realisation -----------------------------------------------------------

    def draw_tau_i_s(self, rng: np.random.Generator) -> float | None:
        """This echo's intrinsic time: log-normal about ``tau_i_s``."""
        if self.tau_i_s is None:
            return None
        if self.tau_log_sigma == 0.0:
            return float(self.tau_i_s)
        return float(self.tau_i_s) * math.exp(float(self.tau_log_sigma) * float(rng.standard_normal()))

    def realise_parts(
        self, t_s: Any, rng: np.random.Generator, *, f_hz: Any = None
    ) -> FadingDraw:
        """Realise one echo on ``t_s`` (and ``f_hz``, if given).

        Arrays are ``(n_t,)`` without frequencies and ``(n_f, n_t)`` with.
        """
        t = np.asarray(t_s, dtype=float).ravel()
        if t.size > 1 and np.any(np.diff(t) <= 0.0):
            raise ValueError("t_s must be strictly increasing")
        f = None if f_hz is None else np.asarray(f_hz, dtype=float).ravel()
        shape = (t.size,) if f is None else (f.size, t.size)
        s_tau, s_spec, s_diff = echo_streams(rng, 3)
        tau_i = self.draw_tau_i_s(s_tau)
        a, b = self._weights()

        if a == 0.0 or self.tau_i_s is None:
            spec = np.ones(shape, dtype=complex)
        else:
            spec = _over_frequency(t, tau_i, self.kind, f, self.specular_coherence_bw_hz, s_spec)
        if b == 0.0:
            diff = np.zeros(shape, dtype=complex)
        else:
            diff = _over_frequency(t, self.tau_d_s, self.kind, f, self.coherence_bw_hz, s_diff)
        gain = a * spec + b * diff
        return FadingDraw(gain=gain, specular=spec, diffuse=diff, tau_i_s=tau_i, tau_d_s=self.tau_d_s)

    def realise(self, t_s: Any, rng: np.random.Generator, *, f_hz: Any = None) -> np.ndarray:
        """The complex gain of one echo.  See :meth:`realise_parts`."""
        return self.realise_parts(t_s, rng, f_hz=f_hz).gain

    def describe(self) -> dict[str, Any]:
        """Parameters for a manifest or a data product."""
        return {
            "tau_i_s": self.tau_i_s,
            "k_factor": self.k_factor,
            "tau_d_s": self.tau_d_s,
            "kind": self.kind,
            "coherence_bw_hz": self.coherence_bw_hz,
            "specular_coherence_bw_hz": self.specular_coherence_bw_hz,
            "tau_log_sigma": self.tau_log_sigma,
            "halo_m": self.halo_m,
        }
