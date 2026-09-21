"""Layer: L2
WP: WP3
Claims: C1
Status: implemented
Requirements: INS-01, INS-10, INS-14, SCI-06

Link budget for a spaceborne HF FMCW topside ionospheric sounder.

This module replaces the legacy ``snr_ionosonde.m`` used in the 2014
TopSideSounder work.  Four things are done differently, and each changes the
answer:

1.  SKY NOISE.  At HF the noise floor is the galactic background, not the
    receiver.  The default reference is the in-situ measurement made by the
    Alouette sounders and reported by Franklin and MacLean (1969): 2e7 K at
    0.6 MHz falling to 4e5 K at 12 MHz, established to +/- 2 dB.  The legacy
    model used a 250 K antenna temperature plus a 15 dB "environmental" term,
    and the ECE 6930 paper used a 300 K receiver noise temperature; both are
    20-35 dB optimistic.

2.  ANTENNA EFFICIENCY IS NOT A SIMPLE LOSS.  An electrically small antenna
    couples to the sky with the same efficiency with which it radiates, so the
    antenna noise temperature is  T_a = eta*T_sky + (1-eta)*T_phys.  When
    eta*T_sky dominates the receiver contribution -- which it does across the
    whole HF band for any efficiency above roughly -30 dB -- the RECEIVE antenna
    efficiency cancels out of the signal-to-noise ratio entirely.  Only the
    transmit efficiency costs anything.  Treating receive efficiency as a
    straight loss, as the legacy model does, double-counts it.

3.  ABSORPTION.  The legacy approach anchors a single measured value (about 5 dB
    two-way at 5 MHz) and scales it as f^-2.  That extrapolation diverges below
    about 2 MHz and produces physically meaningless losses of tens of decibels.
    Here the absorption is integrated properly as 8.686 * (omega/c) * Im(n)
    along the ray path, with the complex refractive index n^2 = 1 - X/(1 - iZ).
    This is finite everywhere, including at the reflection point, and it is the
    term that was responsible for the apparent low-frequency shortfall in
    earlier analyses.

4.  NOISE BANDWIDTH.  In an FMCW sounder the noise that competes with an echo is
    integrated over one DFT bin, not over the swept bandwidth.  The bin width
    follows directly from the range resolution:  B = 2 * K_f * dR / c.  This is
    the entire SWaP argument for the FMCW topology and it must appear explicitly.

Units are SI throughout unless a name ends in _db, _dbw, _dbm, _dbi or _mhz.

THIS MODULE IS A KNOWN SEAM AND IS DELIBERATELY NOT SPLIT ALONG LAYER LINES
---------------------------------------------------------------------------
It sits under ``l2_instrument`` and it contains, mixed together, terms that
belong to three layers: transmitter power, antenna efficiency, receive chain
noise figure and external noise, which are L2; the free-space spreading term and
the absorption integral, which are L1; and the plasma frequency and reflection
condition it uses to decide whether there is an echo at all, which are L0.

It is kept whole because it arrived from the author's earlier work with 37
passing unit tests written against these exact call signatures, and moving code
that those tests cover, at the same time as moving the tests, is the one kind of
refactor whose transcription errors would not be caught.  WP3 leakage work is
the critical path and it builds on this module, so the codebase needs a
trustworthy L2 now more than it needs a tidy one.

The split is scheduled for the start of WP4, in this order: freeze these tests,
extract the L1 terms with the tests duplicated on both sides until they agree,
extract the L0 terms, then delete the re-exports and the duplicated tests.  The
full statement is in ``docs/architecture.md``, section "Known seams",
subsection "The migrated link budget is not split along layer lines".

Migrated verbatim from ``LinkBudget/tss_link_budget.py``.  Only three things
changed: the constants and the ``db``/``undb`` helpers are now imported from
:mod:`tss.constants` and :mod:`tss.units` instead of being defined here, this
docstring, and the seam note above.  No number and no signature moved.

Andrew B. Moss -- PhD dissertation, Utah State University.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from tss.constants import (
    C_LIGHT,
    EPS0,
    K_BOLTZ,
    M_E,
    PLASMA_CONST,
    Q_E,
    T0_K,
    density_from_plasma_frequency,
    plasma_frequency,
)
from tss.units import db, undb

__all__ = [
    "C_LIGHT",
    "K_BOLTZ",
    "Q_E",
    "M_E",
    "EPS0",
    "T0_REF",
    "PLASMA_CONST",
    "db",
    "undb",
    "plasma_frequency",
    "density_from_plasma_frequency",
    "sky_temperature",
    "IonosphereProfile",
    "TabulatedProfile",
    "collision_frequency",
    "refractive_index",
    "reflection_height",
    "absorption_db",
    "Antenna",
    "SounderConfig",
    "BudgetResult",
    "link_budget",
    "snr_vs_frequency",
    "required_tx_power",
    "legacy_snr_ionosonde",
]

#: Reference temperature for noise figure, K.  The same 290 K as
#: :data:`tss.constants.T0_K`; the alias is kept because the migrated tests and
#: the author's scripts spell it ``T0_REF``.
T0_REF = T0_K


# ===========================================================================
# 1.  Sky noise
# ===========================================================================
def sky_temperature(f_mhz, model: str = "franklin"):
    """
    Brightness temperature of the HF sky, in kelvin.

    Parameters
    ----------
    f_mhz : float or array
        Frequency in MHz.
    model : str
        ``"franklin"``   Franklin and MacLean (1969), Proc. IEEE 57(6).  In-situ
                         measurement from Alouette I and II: 2e7 K at 0.6 MHz to
                         4e5 K at 12 MHz, +/- 2 dB.  Log-log interpolated and
                         extrapolated.  This is the default and the recommended
                         reference: it was measured in the exact environment and
                         with the class of antenna in question.
        ``"itu"``        ITU-R P.372 galactic term, Fa = 52 - 23*log10(f_MHz) dB
                         above kT0.  A ground-based model; included for comparison.
        ``"legacy300"``  Flat 300 K -- the assumption in the ECE 6930 paper.
        ``"legacy_sni"`` Flat 250 K raised by 15 dB, i.e. ~7900 K -- the effective
                         assumption in snr_ionosonde.m.

    Notes
    -----
    Solar radio bursts add more than 30 dB sporadically, and bands of ionospheric
    noise below the upper hybrid frequency can add 80 dB (Franklin and MacLean).
    Neither is included here; both are availability questions rather than link
    budget terms.

    Below 0.6 MHz the Franklin curve is an extrapolation of a measurement, not a
    measurement.  See :mod:`tss.l2_instrument.noise` and WP12.
    """
    f = np.asarray(f_mhz, dtype=float)
    if model == "franklin":
        f1, t1 = 0.6, 2.0e7
        f2, t2 = 12.0, 4.0e5
        slope = math.log10(t2 / t1) / math.log10(f2 / f1)     # ~ -1.306
        return t1 * (f / f1) ** slope
    if model == "itu":
        fa_db = 52.0 - 23.0 * np.log10(f)
        return T0_REF * undb(fa_db)
    if model == "legacy300":
        return np.full_like(f, 300.0)
    if model == "legacy_sni":
        return np.full_like(f, 250.0 * undb(15.0))
    raise ValueError(f"unknown sky noise model: {model!r}")


# ===========================================================================
# 2.  Ionosphere: density profile, collisions, absorption
# ===========================================================================
@dataclass
class IonosphereProfile:
    """
    An alpha-Chapman topside used when no external profile is supplied.

    ``plasma_scale_height`` is the PLASMA scale height Hp -- the quantity Meehan
    (2017) measures as 150 +/- 30 km at night and that Stankov et al. (2003)
    show varying from 78 to 244 km depending only on the profiler chosen.  The
    alpha-Chapman topside decays asymptotically with twice the neutral scale
    height, so internally H_chapman = Hp / 2.  Parameterizing by Hp keeps this
    module consistent with the science requirements, which are written in terms
    of Hp.

    Replace with IRI-2020 or NeQuick output using :class:`TabulatedProfile`;
    the absorption integral needs only Ne(h) and Te(h).
    """
    nmf2: float = 5.0e11        # m^-3, F2 peak density (foF2 = 6.35 MHz)
    hmf2: float = 300e3         # m, F2 peak height
    plasma_scale_height: float = 100e3   # m, Hp
    te_peak: float = 1500.0     # K, electron temperature at hmF2
    te_top: float = 3000.0      # K, electron temperature at 1000 km

    @property
    def chapman_scale_height(self):
        return self.plasma_scale_height / 2.0

    def ne(self, h):
        """alpha-Chapman topside above hmF2; below the peak it is not used."""
        h = np.asarray(h, dtype=float)
        z = (h - self.hmf2) / self.chapman_scale_height
        return self.nmf2 * np.exp(0.5 * (1.0 - z - np.exp(-z)))

    def fof2(self):
        """Critical frequency of the F2 layer, Hz."""
        return float(plasma_frequency(self.nmf2))

    def local_plasma_frequency(self, h):
        """Plasma frequency at altitude ``h`` -- the lower bound on sounding."""
        return float(plasma_frequency(self.ne(h)))

    def te(self, h):
        h = np.asarray(h, dtype=float)
        frac = np.clip((h - self.hmf2) / (1000e3 - self.hmf2), 0.0, 1.5)
        return self.te_peak + (self.te_top - self.te_peak) * frac


@dataclass
class TabulatedProfile:
    """Wrap externally computed profiles (IRI, NeQuick, a simulation truth case)."""
    h: np.ndarray          # m, monotonically increasing
    ne_arr: np.ndarray     # m^-3
    te_arr: np.ndarray     # K

    def ne(self, h):
        return np.interp(h, self.h, self.ne_arr)

    def te(self, h):
        return np.interp(h, self.h, self.te_arr)


def collision_frequency(ne, te):
    """
    Electron collision frequency, s^-1.

    Electron-ion Coulomb collisions dominate in the topside F region; the
    Schunk and Nagy form  nu_ei = 54.5 * n_i[cm^-3] / Te^1.5  is used, with
    quasi-neutrality n_i = n_e.  Electron-neutral collisions are negligible above
    about 300 km and are omitted; add them if the profile is extended downward.
    """
    ne_cm3 = np.asarray(ne, dtype=float) * 1e-6
    te = np.asarray(te, dtype=float)
    return 54.5 * ne_cm3 / np.power(te, 1.5)


def refractive_index(f_hz, ne, nu):
    """
    Complex refractive index for a collisional, unmagnetized plasma.

        n^2 = 1 - X / (1 - i Z),    X = (fN/f)^2,   Z = nu / omega

    The magnetic field is neglected.  For a topside sounder this understates the
    difference between the O and X modes but is adequate for a link budget; a
    full Appleton-Hartree treatment belongs in the propagation model (WP2), not
    here.
    """
    omega = 2.0 * math.pi * np.asarray(f_hz, dtype=float)
    x = (plasma_frequency(ne) / np.asarray(f_hz, dtype=float)) ** 2
    z = np.asarray(nu, dtype=float) / omega
    n2 = 1.0 - x / (1.0 - 1j * z)
    return np.sqrt(n2.astype(complex))


def reflection_height(f_hz, profile, h_sat, h_floor=None, n_steps=2000,
                      tol=1e-3):
    """
    Altitude at which a wave of frequency ``f_hz`` reflects, i.e. where the local
    plasma frequency equals the sounding frequency, searching downward from the
    spacecraft.

    Returns
    -------
    float
        Reflection altitude in metres, located by bisection to ``tol`` metres.
    ``"evanescent"``
        The plasma frequency at the spacecraft already exceeds ``f_hz``; the wave
        cannot propagate away from the antenna.  This is the physical floor on
        the sounding band -- the open item that the science requirements flag as
        needing an IRI/NeQuick study to set the minimum sounding frequency.
    None
        The wave penetrates the F2 peak (f > foF2); there is no topside echo.

    The bisection matters: the absorption integrand rises steeply just below the
    reflection point, so a reflection height quantized to a coarse search grid
    puts visible jitter on the absorption-versus-frequency curve.
    """
    h_floor = profile.hmf2 if h_floor is None else h_floor
    g = lambda h: float(plasma_frequency(profile.ne(h))) - f_hz   # noqa: E731

    if g(h_sat) >= 0.0:
        return "evanescent"
    if g(h_floor) < 0.0:
        return None                       # penetrates: f > foF2

    # g is monotonically increasing as altitude decreases toward the peak.
    lo, hi = h_floor, h_sat               # g(lo) >= 0, g(hi) < 0
    for _ in range(int(np.ceil(np.log2(max(h_sat - h_floor, 1.0) / tol))) + 4):
        mid = 0.5 * (lo + hi)
        if g(mid) >= 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def absorption_db(f_hz, profile, h_sat, h_refl=None, n_steps=3000,
                  legacy_model=False, legacy_anchor_db=5.0, legacy_anchor_mhz=5.0):
    """
    Two-way ionospheric absorption in dB for a vertical sounding from ``h_sat``.

    The absorption coefficient is  kappa = (omega/c) * |Im(n)|  in Np/m, integrated
    over the physical path from the spacecraft to the reflection height and
    doubled.  Because the complex refractive index is used, the integrand stays
    finite at the reflection point, where X -> 1 and Im(n) -> sqrt(Z/2).

    Set ``legacy_model=True`` to reproduce the f^-2 anchor scaling used in the
    ECE 6930 paper.  It is retained only so that the difference can be shown; it
    diverges below about 2 MHz and should not be used for results.
    """
    f_hz = float(f_hz)
    if legacy_model:
        return legacy_anchor_db * (legacy_anchor_mhz / (f_hz / 1e6)) ** 2

    if h_refl is None:
        h_refl = reflection_height(f_hz, profile, h_sat)
    if h_refl == "evanescent":
        return float("nan")                 # no propagation; absorption undefined
    if h_refl is None:                      # penetrating ray: integrate to the peak
        h_refl = profile.hmf2

    # Non-uniform grid clustered toward the reflection point.  Im(n) varies as
    # sqrt of the distance from reflection there, so a uniform grid converges
    # slowly and non-monotonically; a cubic stretch fixes it at the same cost.
    s = np.linspace(0.0, 1.0, n_steps) ** 3          # dense at s -> 0
    h = h_refl + (h_sat - h_refl) * s
    ne = profile.ne(h)
    te = profile.te(h)
    nu = collision_frequency(ne, te)
    n_cplx = refractive_index(f_hz, ne, nu)
    omega = 2.0 * math.pi * f_hz
    kappa = (omega / C_LIGHT) * np.abs(np.imag(n_cplx))       # Np/m
    one_way_np = float(np.trapezoid(kappa, h))               # h increasing
    return 2.0 * abs(one_way_np) * 8.685889638


# ===========================================================================
# 3.  Antenna
# ===========================================================================
@dataclass
class Antenna:
    """
    An electrically small antenna described by directivity and radiation
    efficiency rather than by a single lumped gain, because the two enter the
    noise calculation differently.

    ``efficiency_db_at`` maps frequency in MHz to radiation efficiency in dB
    (negative).  The default reproduces the gain assumption used in the ECE 6930
    paper and the TISC proposal -- about -10 dBi at 1 MHz rising to -4 dBi at
    20 MHz -- for a small loop of 1.76 dBi directivity.
    """
    directivity_dbi: float = 1.76                # small loop / short dipole
    eff_ref_mhz: Sequence[float] = (1.0, 20.0)
    eff_ref_db: Sequence[float] = (-11.76, -5.76)
    physical_temperature: float = 290.0          # K, for the lossy part

    def efficiency_db(self, f_mhz):
        return np.interp(np.log10(np.asarray(f_mhz, dtype=float)),
                         np.log10(np.asarray(self.eff_ref_mhz, dtype=float)),
                         np.asarray(self.eff_ref_db, dtype=float))

    def efficiency(self, f_mhz):
        return undb(self.efficiency_db(f_mhz))

    def gain_dbi(self, f_mhz):
        return self.efficiency_db(f_mhz) + self.directivity_dbi

    def noise_temperature(self, f_mhz, t_sky):
        """T_a = eta*T_sky + (1-eta)*T_phys, referred to the antenna terminals."""
        eta = self.efficiency(f_mhz)
        return eta * t_sky + (1.0 - eta) * self.physical_temperature


# ===========================================================================
# 4.  The budget
# ===========================================================================
@dataclass
class SounderConfig:
    """Instrument and geometry.  Defaults follow the TISC point design corrected
    to the requirements of TSS_Science_Requirements_RevA."""
    # --- waveform / processing -------------------------------------------------
    sweep_rate: float = 5.6e6        # Hz/s, K_f
    range_resolution: float = 1000.0  # m, delta R  (SCI-06 threshold is 10 km)
    coherent_sweeps: int = 1         # M, coherently integrated sweeps
    # --- geometry --------------------------------------------------------------
    h_sat: float = 600e3             # m, spacecraft altitude
    true_range: float = 300e3        # m, spacecraft to reflection point
    # --- power and losses ------------------------------------------------------
    tx_power: float = 1.0            # W
    polarization_loss_db: float = -6.0   # -3 dB O/X split, -3 dB elliptical->linear
    extra_loss_db: float = 0.0       # feed, mismatch, anything else (negative)
    rx_noise_figure_db: float = 3.0  # receiver chain
    feed_loss_db: float = 0.0        # between antenna and LNA (negative)
    # --- antennas --------------------------------------------------------------
    tx_antenna: Antenna = field(default_factory=Antenna)
    rx_antenna: Antenna = field(default_factory=Antenna)
    # --- environment -----------------------------------------------------------
    sky_model: str = "franklin"
    ionosphere: IonosphereProfile = field(default_factory=IonosphereProfile)
    absorption_legacy: bool = False

    @property
    def noise_bandwidth(self):
        """Noise bandwidth of one range cell, Hz:  B = 2 K_f dR / c."""
        return 2.0 * self.sweep_rate * self.range_resolution / C_LIGHT

    @property
    def dwell_time(self):
        """Coherent dwell of one transform, s."""
        return 1.0 / self.noise_bandwidth

    def max_beat_frequency(self, r_max):
        return 2.0 * self.sweep_rate * r_max / C_LIGHT

    def min_sample_rate(self, r_max):
        return 2.0 * self.max_beat_frequency(r_max)


@dataclass
class BudgetResult:
    """Every term, in dB, so the result can be printed as a budget table."""
    f_mhz: float
    tx_power_dbw: float
    tx_gain_dbi: float
    eirp_dbw: float
    path_loss_db: float
    absorption_db: float
    polarization_loss_db: float
    extra_loss_db: float
    rx_gain_dbi: float
    rx_power_dbw: float
    sky_temperature_k: float
    antenna_temperature_k: float
    system_temperature_k: float
    noise_bandwidth_hz: float
    noise_power_dbw: float
    snr_single_db: float
    integration_gain_db: float
    snr_db: float
    reflection_height_km: float | None
    rx_efficiency_db: float
    status: str = "ok"          # "ok" | "evanescent" | "penetrates"

    def table(self):
        rows = [
            ("Frequency", f"{self.f_mhz:.3f}", "MHz"),
            ("Transmit power", f"{self.tx_power_dbw:+.2f}", "dBW"),
            ("Transmit antenna gain", f"{self.tx_gain_dbi:+.2f}", "dBi"),
            ("EIRP", f"{self.eirp_dbw:+.2f}", "dBW"),
            ("Path loss (2R free space)", f"{self.path_loss_db:+.2f}", "dB"),
            ("Ionospheric absorption (2-way)", f"{self.absorption_db:+.2f}", "dB"),
            ("Polarization loss", f"{self.polarization_loss_db:+.2f}", "dB"),
            ("Other losses", f"{self.extra_loss_db:+.2f}", "dB"),
            ("Receive antenna gain", f"{self.rx_gain_dbi:+.2f}", "dBi"),
            ("Received signal power", f"{self.rx_power_dbw:+.2f}", "dBW"),
            ("", "", ""),
            ("Sky temperature", f"{self.sky_temperature_k:.3e}", "K"),
            ("Receive efficiency", f"{self.rx_efficiency_db:+.2f}", "dB"),
            ("Antenna temperature", f"{self.antenna_temperature_k:.3e}", "K"),
            ("System temperature", f"{self.system_temperature_k:.3e}", "K"),
            ("Noise bandwidth (one range cell)", f"{self.noise_bandwidth_hz:.2f}", "Hz"),
            ("Noise power", f"{self.noise_power_dbw:+.2f}", "dBW"),
            ("", "", ""),
            ("SNR, single sweep", f"{self.snr_single_db:+.2f}", "dB"),
            ("Coherent integration gain", f"{self.integration_gain_db:+.2f}", "dB"),
            ("SNR", f"{self.snr_db:+.2f}", "dB"),
            ("Status", self.status, ""),
        ]
        w = max(len(r[0]) for r in rows)
        out = []
        for name, val, unit in rows:
            out.append("" if not name else f"{name:<{w}}  {val:>12}  {unit}")
        return "\n".join(out)


def link_budget(f_hz, cfg: SounderConfig) -> BudgetResult:
    """Compute the full link budget at one sounding frequency."""
    f_hz = float(f_hz)
    f_mhz = f_hz / 1e6
    lam = C_LIGHT / f_hz

    gt = float(cfg.tx_antenna.gain_dbi(f_mhz))
    gr = float(cfg.rx_antenna.gain_dbi(f_mhz))
    pt_dbw = db(cfg.tx_power)
    eirp = pt_dbw + gt

    # Specular reflection from a smooth overdense layer: the image source sits at
    # twice the true range and the layer does not focus, so the geometric term is
    # free-space spreading over a path of 2R.
    path = 20.0 * math.log10(lam / (4.0 * math.pi * 2.0 * cfg.true_range))

    h_refl = reflection_height(f_hz, cfg.ionosphere, cfg.h_sat)
    if h_refl == "evanescent":
        status, h_refl_val = "evanescent", None
    elif h_refl is None:
        status, h_refl_val = "penetrates", None
    else:
        status, h_refl_val = "ok", h_refl
    a_db = -absorption_db(f_hz, cfg.ionosphere, cfg.h_sat, h_refl,
                          legacy_model=cfg.absorption_legacy)

    pr = (eirp + path + a_db + cfg.polarization_loss_db + cfg.extra_loss_db
          + gr + cfg.feed_loss_db)

    t_sky = float(sky_temperature(f_mhz, cfg.sky_model))
    t_ant = float(cfg.rx_antenna.noise_temperature(f_mhz, t_sky))
    # Feed loss attenuates the antenna noise and adds its own thermal noise.
    l_feed = undb(cfg.feed_loss_db)                      # <= 1
    t_feed = (1.0 / l_feed - 1.0) * cfg.rx_antenna.physical_temperature if l_feed < 1 else 0.0
    t_rx = (undb(cfg.rx_noise_figure_db) - 1.0) * T0_REF
    t_sys = t_ant + t_feed + t_rx

    b = cfg.noise_bandwidth
    n_dbw = db(K_BOLTZ * t_sys * b)

    if status != "ok":
        # Evanescent: the wave never leaves the antenna.  Penetrating: there is
        # no topside reflection, only a ground echo, whose geometry is a
        # different calculation and is out of scope here.
        pr = float("nan")
    snr1 = pr - n_dbw
    gi = 10.0 * math.log10(max(cfg.coherent_sweeps, 1))
    return BudgetResult(
        f_mhz=f_mhz, tx_power_dbw=pt_dbw, tx_gain_dbi=gt, eirp_dbw=eirp,
        path_loss_db=path, absorption_db=a_db,
        polarization_loss_db=cfg.polarization_loss_db, extra_loss_db=cfg.extra_loss_db,
        rx_gain_dbi=gr, rx_power_dbw=pr, sky_temperature_k=t_sky,
        antenna_temperature_k=t_ant, system_temperature_k=t_sys,
        noise_bandwidth_hz=b, noise_power_dbw=n_dbw,
        snr_single_db=snr1, integration_gain_db=gi, snr_db=snr1 + gi,
        reflection_height_km=None if h_refl_val is None else h_refl_val / 1e3,
        rx_efficiency_db=float(cfg.rx_antenna.efficiency_db(f_mhz)),
        status=status,
    )


def snr_vs_frequency(f_mhz_array, cfg: SounderConfig):
    """SNR in dB at each frequency in ``f_mhz_array``."""
    return np.array([link_budget(f * 1e6, cfg).snr_db for f in np.atleast_1d(f_mhz_array)])


def required_tx_power(f_hz, cfg: SounderConfig, target_snr_db=10.0):
    """
    Transmit power in watts needed to reach ``target_snr_db``.  SNR is linear in
    transmit power in dB, so one evaluation suffices.
    """
    ref = link_budget(f_hz, cfg)
    return cfg.tx_power * undb(target_snr_db - ref.snr_db)


# ===========================================================================
# 5.  Legacy model, retained as a regression case
# ===========================================================================
def legacy_snr_ionosonde(kf, r, r_max, f, pt):
    """
    Faithful Python port of the 2014 ``snr_ionosonde.m`` (J. Pratt).

    Kept so that the difference between the old and new models can be shown
    term by term.  Do not use for results: it takes the noise bandwidth as the
    full swept bandwidth rather than one range cell, treats receive antenna
    efficiency as a straight loss while assuming a 250 K antenna temperature,
    and adds a 15 dB "environmental noise" term in place of a sky model.

    Returns (snr_db, rx_power_db, noise_power_db) exactly as the original did.
    """
    tref, k = 290.0, 1.38e-23
    bw = 2.0 * kf * r_max / C_LIGHT
    lt, gt = -3.0, -6.0
    eirp = 10.0 * np.log10(pt) + lt + gt
    lam = C_LIGHT / np.asarray(f, dtype=float)
    ls = 20.0 * np.log10(lam / (4.0 * np.pi * 2.0 * np.asarray(r, dtype=float)))
    la, lf, lp = -0.3, -5.0, -6.51
    gr, ta, lr = gt, 250.0, -3.0
    tr = (undb(-lr) - 1.0) * tref
    gl, fl = 20.0, 3.0
    tl = (undb(fl) - 1.0) * tref
    len_env = 15.0

    n_a = 10.0 * np.log10(k * ta * bw) + len_env
    n_f = 10.0 * np.log10(undb(n_a + lr) + undb(lr) * k * tr * bw)
    n_r = 10.0 * np.log10(undb(n_f + gl) + undb(gl) * k * tl * bw)

    pwr_a = eirp + gr + ls + la + lf + lp
    pwr_r = pwr_a + lr + gl
    return pwr_r - n_r, pwr_r, n_r
