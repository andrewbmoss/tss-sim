"""Layer: L2
WP: WP3
Claims: C1
Status: implemented
Requirements: INS-01, INS-14

The electrically short dipole: physical length to radiation efficiency, and the
antenna noise temperature that goes with it.

The band edge is a requirement on the instrument.  The antenna that meets it is
a mission decision.  This module turns a physical length into the efficiency
penalty that decision costs, so a requirement can state the ideal size and the
price of anything smaller.

THE MODEL
---------
A centre-fed dipole of tip-to-tip length L and conductor radius a, short
compared with a wavelength, fed through a series matching inductor.

    radiation resistance     R_r = 20 pi^2 (L/lambda)^2      (triangular current)
    conductor loss           R_c = L R_s / (6 pi a),  R_s = sqrt(pi f mu0 / sigma)
    terminal capacitance     C_a = pi eps0 L / (ln(L/2a) - 1)
    matching coil loss       R_m = |X_a| / Q,  X_a = -1/(omega C_a)
    radiation efficiency     eta = R_r / (R_r + R_c + R_m)

Two things this deliberately does not model: the mutual impedance of a deployed
structure against a small conducting bus, and the plasma sheath, which at these
frequencies changes the terminal impedance materially.  Both make the real
number worse than this one, so treat the result as an optimistic bound.

DIRECTIVITY AND EFFICIENCY ARE CARRIED SEPARATELY
-------------------------------------------------
They are never lumped into a single gain figure, because they enter the noise
calculation differently and lumping them is exactly the error the legacy model
made.

An electrically small antenna couples to the external sky with the same
efficiency with which it radiates.  So the antenna noise temperature referred to
its terminals is

    T_a = eta * T_sky + (1 - eta) * T_phys

while the gain is

    G = D + eta            (in decibels)

Signal power at the terminals goes as ``eta``.  Noise power at the terminals
also goes as ``eta``, as long as the first term dominates the second and the
receiver's own contribution.  The ratio of the two is therefore independent of
``eta``, and the RECEIVE efficiency cancels out of the signal-to-noise ratio
entirely.  Only the TRANSMIT efficiency costs anything.  Dropping the radiation
efficiency of both antennas by 3 dB costs 3 dB of SNR, not 6 dB.  A model that
treats receive efficiency as a straight loss double-counts it.

This is the single most consequential correction in the WP3 link budget and
``tests/unit/test_l2_linkbudget.py::test_receive_efficiency_cancels_while_the_sky_dominates``
is the test that holds it.

WHERE THE CANCELLATION STOPS
----------------------------
It holds only while the sky term dominates the rest of the system temperature.
The rest is the lossy part of the antenna plus the receiver:

    (1 - eta) * T_phys + (undb(NF) - 1) * T0

At the baseline noise figure of 3 dB and a physical temperature of 290 K, and
with ``eta`` small enough that ``(1 - eta) -> 1``, that is

    290 + (10**0.3 - 1) * 290 = 290 + 288.6 = 578.6 K

which is the ~580 K quoted in the requirement discussion.  Sky dominance means
``eta * T_sky >> 580 K``.  Below roughly -30 dB of efficiency at the top of the
band, that stops being true and the receiver sets the floor;
``test_receiver_noise_takes_over_at_very_low_efficiency`` holds that end.

THE 5 m CANDIDATE
-----------------
The author's reported result is that the 5 m candidate clears the sky-dominance
threshold by only 0.5 dB on the measured sky, that is at the bottom of the band
where the measurement exists.  That is a thin margin and it should be read as
one: half a decibel is inside the +/- 2 dB the Franklin and MacLean reference
itself carries, so on that number sky dominance at the band bottom is asserted
rather than demonstrated.

Recomputing it here with :func:`short_dipole_efficiency` at its default
conductor radius of 5 mm and matching-coil Q of 150 gives a larger margin, about
10 dB at 0.6 MHz and about 6 dB at the 0.33 MHz band edge.  The two are
reconcilable only with a more pessimistic antenna than this module's defaults,
roughly a factor of three lower coil Q or about 5 dB of additional
matching-network loss.  Which of those the reported figure assumed is not
recorded in the source this module was migrated from, and it is flagged here
rather than quietly resolved in favour of the more comfortable number.  Use
:func:`sky_dominance_margin_db` to recompute the margin for whatever antenna a
configuration actually specifies.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from tss.constants import C_LIGHT, EPS0, MU0, T0_K
from tss.units import db, undb

__all__ = [
    "SIGMA_CU",
    "SHORT_DIPOLE_DIRECTIVITY_DBI",
    "ShortDipole",
    "short_dipole_efficiency",
    "resonant_length_m",
    "antenna_noise_temperature_k",
    "non_sky_temperature_k",
    "sky_dominance_margin_db",
]

#: Conductivity of copper, S/m.
SIGMA_CU = 5.8e7

#: Directivity of a short dipole, dBi.  ``10*log10(1.5)``.  Lossless, and it is
#: the same 1.76 for a small loop, which is why the link budget's ``Antenna``
#: default carries the same figure.
SHORT_DIPOLE_DIRECTIVITY_DBI = 1.76


def short_dipole_efficiency(
    f_hz: np.ndarray | float,
    length_m: float,
    radius_m: float = 5e-3,
    sigma: float = SIGMA_CU,
    coil_q: float = 150.0,
    include_match: bool = True,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Radiation efficiency of a short dipole, as a ratio, with its terms.

    Migrated from ``ISISValidation/small_antenna.py`` with the arithmetic
    unchanged.  The second return value carries ``r_rad``, ``r_cond``,
    ``r_match`` and ``lam`` so that a trade study can see which loss dominates
    rather than only the total, which matters because the conductor loss and the
    matching-coil loss respond to completely different design choices.

    The ``max(ln(L/2a) - 1, 0.1)`` floor in the terminal capacitance is kept
    from the original.  It stops the capacitance running away for a fat, short
    conductor where the thin-wire expression has stopped being valid; it is a
    guard, not a model, and a configuration that reaches it is outside the
    model's range.
    """
    f = np.asarray(f_hz, dtype=float)
    lam = C_LIGHT / f
    length = float(length_m)
    radius = float(radius_m)

    r_rad = 20.0 * np.pi**2 * (length / lam) ** 2
    r_s = np.sqrt(np.pi * f * MU0 / sigma)
    r_c = length * r_s / (6.0 * np.pi * radius)

    if include_match:
        c_a = np.pi * EPS0 * length / max(np.log(length / (2.0 * radius)) - 1.0, 0.1)
        x_a = 1.0 / (2.0 * np.pi * f * c_a)   # magnitude of the capacitive reactance
        r_m = x_a / coil_q
    else:
        r_m = np.zeros_like(f)

    eta = r_rad / (r_rad + r_c + r_m)
    return eta, {"r_rad": r_rad, "r_cond": r_c, "r_match": r_m, "lam": lam}


def resonant_length_m(f_hz: np.ndarray | float) -> np.ndarray:
    """Tip-to-tip length of a resonant half-wave dipole, m.

    Reported alongside the chosen length so that a requirement can state how
    electrically short the antenna is.  At 0.33 MHz this is 454 m, against a 5 m
    candidate, which is the scale of the problem in one number.
    """
    return C_LIGHT / np.asarray(f_hz, dtype=float) / 2.0


def antenna_noise_temperature_k(
    eta: np.ndarray | float,
    t_sky_k: np.ndarray | float,
    t_phys_k: float = 290.0,
) -> np.ndarray:
    """``T_a = eta*T_sky + (1-eta)*T_phys``, referred to the antenna terminals.

    This is the relation the whole receive-efficiency cancellation rests on.  It
    is the same expression as ``Antenna.noise_temperature`` in
    :mod:`tss.l2_instrument.linkbudget`; it is written again here because this
    module is where the argument for it lives, and the two are checked against
    each other in ``tests/unit/test_l2_linkbudget.py``.
    """
    e = np.asarray(eta, dtype=float)
    return e * np.asarray(t_sky_k, dtype=float) + (1.0 - e) * float(t_phys_k)


def non_sky_temperature_k(
    eta: np.ndarray | float,
    noise_figure_db: float = 3.0,
    t_phys_k: float = 290.0,
) -> np.ndarray:
    """Everything in the system temperature that is not the sky, K.

    ``(1 - eta) * T_phys + (undb(NF) - 1) * T0``.  At the baseline 3 dB noise
    figure and small ``eta`` this is the ~580 K the sky term has to beat for the
    receive-efficiency cancellation to hold.
    """
    e = np.asarray(eta, dtype=float)
    return (1.0 - e) * float(t_phys_k) + (undb(noise_figure_db) - 1.0) * T0_K


def sky_dominance_margin_db(
    eta: np.ndarray | float,
    t_sky_k: np.ndarray | float,
    noise_figure_db: float = 3.0,
    t_phys_k: float = 290.0,
) -> np.ndarray:
    """How far the coupled sky noise sits above everything else, in dB.

    ``10 log10( eta*T_sky / non_sky_temperature )``.  Positive means the sky
    dominates and the receive-efficiency cancellation applies.  A margin inside
    the +/- 2 dB the Franklin and MacLean reference carries means the
    cancellation is being asserted rather than demonstrated at that frequency,
    and a result that depends on it should say so.
    """
    e = np.asarray(eta, dtype=float)
    return db(e * np.asarray(t_sky_k, dtype=float) / non_sky_temperature_k(e, noise_figure_db, t_phys_k))


class ShortDipole:
    """An electrically short dipole, described by length rather than by gain.

    Parameters
    ----------
    length_m
        Tip-to-tip length.
    radius_m
        Conductor radius.
    directivity_dbi
        Carried separately from the efficiency and never folded into it.  The
        default is the 1.76 dBi of a lossless short dipole.
    coil_q
        Quality factor of the series matching inductor.  This is usually the
        dominant loss at the bottom of the band, so it is a design parameter and
        not a detail.
    sigma
        Conductor conductivity.
    physical_temperature_k
        Physical temperature of the lossy part, for the antenna noise
        temperature.
    include_match
        False removes the matching inductor entirely, which is not a buildable
        antenna and is useful for isolating how much of the loss the match
        costs.

    The class exposes the same three methods the link budget's ``Antenna``
    exposes, ``efficiency_db``, ``gain_dbi`` and ``noise_temperature``, so it can
    be passed wherever that dataclass is expected.  The difference is that this
    one derives the efficiency from a physical length instead of interpolating a
    pair of reference points, which is what makes the length a design variable a
    trade study can sweep.
    """

    def __init__(
        self,
        length_m: float = 5.0,
        radius_m: float = 5.0e-3,
        *,
        directivity_dbi: float = SHORT_DIPOLE_DIRECTIVITY_DBI,
        coil_q: float = 150.0,
        sigma: float = SIGMA_CU,
        physical_temperature_k: float = 290.0,
        include_match: bool = True,
    ) -> None:
        if length_m <= 0.0:
            raise ValueError("length_m must be positive")
        if radius_m <= 0.0:
            raise ValueError("radius_m must be positive")
        if coil_q <= 0.0:
            raise ValueError("coil_q must be positive")
        self.length_m = float(length_m)
        self.radius_m = float(radius_m)
        self.directivity_dbi = float(directivity_dbi)
        self.coil_q = float(coil_q)
        self.sigma = float(sigma)
        self.physical_temperature = float(physical_temperature_k)
        self.include_match = bool(include_match)

    # -- efficiency ---------------------------------------------------------

    def efficiency(self, f_mhz: np.ndarray | float) -> np.ndarray:
        """Radiation efficiency as a ratio, from the physical model."""
        eta, _ = short_dipole_efficiency(
            np.asarray(f_mhz, dtype=float) * 1.0e6,
            self.length_m,
            self.radius_m,
            self.sigma,
            self.coil_q,
            self.include_match,
        )
        return eta

    def efficiency_db(self, f_mhz: np.ndarray | float) -> np.ndarray:
        """Radiation efficiency in dB.  Negative."""
        return db(self.efficiency(f_mhz))

    def terms(self, f_mhz: np.ndarray | float) -> dict[str, Any]:
        """The resistance terms, so a trade study can see which loss dominates."""
        _, detail = short_dipole_efficiency(
            np.asarray(f_mhz, dtype=float) * 1.0e6,
            self.length_m,
            self.radius_m,
            self.sigma,
            self.coil_q,
            self.include_match,
        )
        return detail

    # -- gain and noise -----------------------------------------------------

    def gain_dbi(self, f_mhz: np.ndarray | float) -> np.ndarray:
        """``G = D + eta`` in decibels.  The two are added here and nowhere else."""
        return self.efficiency_db(f_mhz) + self.directivity_dbi

    def noise_temperature(
        self, f_mhz: np.ndarray | float, t_sky: np.ndarray | float
    ) -> np.ndarray:
        """``T_a = eta*T_sky + (1-eta)*T_phys``."""
        return antenna_noise_temperature_k(
            self.efficiency(f_mhz), t_sky, self.physical_temperature
        )

    def sky_dominance_margin_db(
        self,
        f_mhz: np.ndarray | float,
        t_sky: np.ndarray | float,
        noise_figure_db: float = 3.0,
    ) -> np.ndarray:
        """Margin by which the coupled sky beats everything else, dB."""
        return sky_dominance_margin_db(
            self.efficiency(f_mhz), t_sky, noise_figure_db, self.physical_temperature
        )

    def electrical_length(self, f_mhz: np.ndarray | float) -> np.ndarray:
        """``L / lambda``.  Below about 0.1 the short-dipole forms apply."""
        lam = C_LIGHT / (np.asarray(f_mhz, dtype=float) * 1.0e6)
        return self.length_m / lam

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable record for the run manifest."""
        return {
            "kind": "electrically_short_dipole",
            "length_m": self.length_m,
            "radius_m": self.radius_m,
            "directivity_dbi": self.directivity_dbi,
            "coil_q": self.coil_q,
            "sigma_s_per_m": self.sigma,
            "physical_temperature_k": self.physical_temperature,
            "include_match": self.include_match,
            "note": (
                "Optimistic bound: mutual impedance against the bus and the plasma "
                "sheath are not modelled and both make the real number worse."
            ),
        }

    @classmethod
    def from_config(cls, cfg: Mapping[str, Any]) -> ShortDipole:
        """Build from the ``antenna:`` block of an instrument configuration.

        The baseline configuration states ``radiation_efficiency_db`` directly as
        well as a length.  That stated figure is a requirement threshold and this
        class computes an efficiency from physics, so the two can disagree.  When
        they do, the stated figure wins and is carried as
        ``stated_efficiency_db``, because a run against the baseline
        configuration has to be a run against the threshold the requirement set
        actually writes down.  :func:`efficiency_db` is unaffected; the caller
        decides which to use and :mod:`tss.l2_instrument.baseline` records which
        it used.
        """
        antenna = cls(
            length_m=float(cfg.get("length_m", 5.0)),
            radius_m=float(cfg.get("radius_m", 5.0e-3)),
            directivity_dbi=float(cfg.get("directivity_dbi", SHORT_DIPOLE_DIRECTIVITY_DBI)),
            coil_q=float(cfg.get("coil_q", 150.0)),
            physical_temperature_k=float(cfg.get("physical_temperature_k", 290.0)),
        )
        stated = cfg.get("radiation_efficiency_db")
        antenna.stated_efficiency_db = None if stated is None else float(stated)
        return antenna

    #: Set by :meth:`from_config` when the configuration states an efficiency.
    stated_efficiency_db: float | None = None
