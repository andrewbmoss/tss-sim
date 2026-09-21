"""Layer: L2
WP: WP3, WP12
Claims: C1
Status: partial
Requirements: INS-01

The external noise reference, and the open question at the bottom of the band.

At HF the noise floor is the sky, not the receiver.  Everything about whether a
sub-watt transmitter closes the link depends on which sky reference is used, and
the references disagree with each other by tens of decibels below a megahertz.
This module holds the reference, states its provenance, and marks exactly where
it stops being a measurement.

THE PRIMARY REFERENCE
---------------------
Franklin and MacLean (1969), Proc. IEEE 57(6).  An in-situ measurement made by
the Alouette sounders: 2e7 K at 0.6 MHz falling to 4e5 K at 12 MHz, established
to +/- 2 dB.  It is the right reference for this instrument for a reason that
has nothing to do with it being convenient.  It was measured in the exact
environment, above the ionosphere rather than below it, and with the class of
antenna in question, an electrically short dipole on a small spacecraft.  A
ground-based model measures a sky seen through the ionosphere, which is a
different thing at these frequencies.

The implementation is :func:`tss.l2_instrument.linkbudget.sky_temperature`, and
it is re-exported here rather than reimplemented.  The link budget is the module
with the 21 tests on it and there is no version of this that is improved by
having two implementations of the same power law.

WHERE IT STOPS BEING A MEASUREMENT
----------------------------------
Below 0.6 MHz the Franklin curve is an extrapolation.  The band edge the
instrument is specified to, 0.33 MHz, is below that.  So the noise floor at the
bottom of the band, and therefore the link margin at the bottom of the band, and
therefore the part of C1 that concerns the lowest sounding frequencies, rests on
a power law continued past its data.  This is the open part of C1 and it is
stated as open.  :data:`MEASURED_FLOOR_MHZ` is the frequency below which
:func:`is_extrapolated` returns True, and any L5 verdict on INS-01 coverage
below it carries a flag saying so.

WP12 exists to replace the extrapolation with something better.

THE CANDIDATE EXTENSION
-----------------------
Cane (1979) and Dulk et al. (2001) give the galactic background as a brightness

    I(nu) = I_g nu^-0.52 [1 - exp(-tau)] / tau  +  I_eg nu^-0.80 exp(-tau)

    tau = 5.0 nu^-2.1

with nu in MHz, I_g = 2.48e-20 and I_eg = 1.06e-20 W m^-2 Hz^-1 sr^-1.  The
first term is the galactic disc emission attenuated by free-free absorption in
the interstellar medium, the second the extragalactic component seen through it,
and ``tau`` is the free-free optical depth, which is what turns the spectrum
over near a few megahertz and makes an extrapolated power law the wrong shape
below it.

That turnover is the whole reason this is a candidate.  A power law fitted
between 0.6 and 12 MHz and continued downward keeps rising; the Cane and Dulk
form peaks around 3 MHz and then falls, because the galaxy becomes optically
thick to its own emission.  The two therefore diverge in the direction that
matters, and the extrapolation is pessimistic rather than optimistic if the
Cane and Dulk shape is right down there.

It is offered as a candidate and not adopted, for two reasons.  It is a model of
the galactic background alone, and the Alouette measurement includes whatever
else is present above the ionosphere, which at these frequencies includes
auroral kilometric radiation and the ionospheric noise bands below the upper
hybrid frequency that Franklin and MacLean report can add 80 dB sporadically.
And it has never been checked against an in-situ measurement in the 0.3 to 0.6
MHz decade from a topside platform, which is precisely the comparison WP12 is
for.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from tss.constants import C_LIGHT, K_BOLTZ
from tss.l2_instrument.linkbudget import sky_temperature

__all__ = [
    "sky_temperature",
    "MEASURED_FLOOR_MHZ",
    "MEASURED_CEILING_MHZ",
    "FRANKLIN_UNCERTAINTY_DB",
    "CANE_I_G",
    "CANE_I_EG",
    "is_extrapolated",
    "cane_dulk_brightness",
    "cane_dulk_optical_depth",
    "cane_dulk_temperature_k",
    "brightness_to_temperature_k",
    "noise_density_dbm_per_hz",
    "describe_reference",
]

#: Lowest frequency at which the Franklin and MacLean reference is a
#: measurement.  Below this it is an extrapolation and the flag says so.
MEASURED_FLOOR_MHZ = 0.6

#: Highest frequency of the measurement.  Above it the curve is again an
#: extrapolation, but a far less consequential one: the band top is 20 MHz
#: against a measurement that reaches 12, and the sky is not the limiting noise
#: source up there in any case.
MEASURED_CEILING_MHZ = 12.0

#: Stated uncertainty on the Franklin and MacLean values, dB.
FRANKLIN_UNCERTAINTY_DB = 2.0

#: Galactic disc term, W m^-2 Hz^-1 sr^-1 at 1 MHz.
CANE_I_G = 2.48e-20

#: Extragalactic term, W m^-2 Hz^-1 sr^-1 at 1 MHz.
CANE_I_EG = 1.06e-20


def is_extrapolated(f_mhz: np.ndarray | float) -> np.ndarray:
    """True where the primary reference is extrapolated rather than measured.

    A result computed at a frequency where this is True carries the
    ``noise_reference_extrapolated`` flag into L5, per the "Sub-megahertz noise
    is extrapolated" seam in ``docs/architecture.md``.
    """
    f = np.asarray(f_mhz, dtype=float)
    return (f < MEASURED_FLOOR_MHZ) | (f > MEASURED_CEILING_MHZ)


def cane_dulk_optical_depth(f_mhz: np.ndarray | float) -> np.ndarray:
    """Free-free optical depth of the interstellar medium, ``tau = 5.0 nu^-2.1``.

    Dimensionless, with ``nu`` in MHz.  It passes through unity near 2.2 MHz,
    which is where the spectrum turns over.
    """
    return 5.0 * np.asarray(f_mhz, dtype=float) ** -2.1


def cane_dulk_brightness(f_mhz: np.ndarray | float) -> np.ndarray:
    """Galactic background brightness, W m^-2 Hz^-1 sr^-1.

        ``I(nu) = I_g nu^-0.52 [1 - exp(-tau)]/tau + I_eg nu^-0.80 exp(-tau)``

    The ``[1 - exp(-tau)]/tau`` factor is the radiative transfer solution for a
    uniformly emitting, uniformly absorbing slab seen from inside it: it tends
    to one for small ``tau``, so the galactic term is the bare power law at high
    frequency, and to ``1/tau`` for large ``tau``, so the term saturates and then
    falls at low frequency.  Writing it that way rather than as a fitted broken
    power law is the point of the form, because the turnover then comes out of
    the physics rather than being placed by hand.
    """
    nu = np.asarray(f_mhz, dtype=float)
    tau = cane_dulk_optical_depth(nu)
    # For very small tau the ratio is 0/0 in floating point.  The limit is 1.
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        slab = np.where(tau > 1.0e-8, (1.0 - np.exp(-tau)) / np.where(tau > 0.0, tau, 1.0), 1.0)
    return CANE_I_G * nu**-0.52 * slab + CANE_I_EG * nu**-0.80 * np.exp(-tau)


def brightness_to_temperature_k(
    brightness: np.ndarray | float, f_mhz: np.ndarray | float
) -> np.ndarray:
    """Rayleigh-Jeans brightness temperature, K, from ``W m^-2 Hz^-1 sr^-1``.

    ``T = I c^2 / (2 k nu^2)``.  Rayleigh-Jeans is exact to well within the
    quoted uncertainties everywhere in this band; ``h nu / k T`` is of order
    1e-11 at 1 MHz and 1e5 K.
    """
    nu_hz = np.asarray(f_mhz, dtype=float) * 1.0e6
    return np.asarray(brightness, dtype=float) * C_LIGHT**2 / (2.0 * K_BOLTZ * nu_hz**2)


def cane_dulk_temperature_k(f_mhz: np.ndarray | float) -> np.ndarray:
    """Brightness temperature of the Cane and Dulk galactic background, K.

    This is the sub-megahertz extension CANDIDATE for WP12.  It is not the
    package default and no scored run uses it, because the primary reference is
    an in-situ measurement and this is a model.  Compare the two with
    :func:`describe_reference` before using either below 0.6 MHz.
    """
    return brightness_to_temperature_k(cane_dulk_brightness(f_mhz), f_mhz)


def noise_density_dbm_per_hz(t_sys_k: np.ndarray | float) -> np.ndarray:
    """External noise density in dBm/Hz from a system temperature.

    ``10 log10(k T * 1000)``.  This is what goes into
    ``ReceivedSignal.noise_floor_dbm_hz`` so that a received signal carries the
    noise reference it was generated against, and a result can be traced back to
    Franklin and MacLean rather than to an unnamed number.
    """
    t = np.asarray(t_sys_k, dtype=float)
    return 10.0 * np.log10(K_BOLTZ * t * 1.0e3)


def describe_reference(f_mhz: np.ndarray | float) -> Mapping[str, Any]:
    """Both references at the given frequencies, with the extrapolation flag.

    Returned as a JSON-serialisable mapping so it can go straight into a run
    manifest.  ``ratio_db`` is the primary reference over the candidate, so a
    positive number means the extrapolated Franklin curve is the more
    pessimistic of the two.
    """
    f = np.atleast_1d(np.asarray(f_mhz, dtype=float))
    primary = np.asarray(sky_temperature(f, "franklin"), dtype=float)
    candidate = np.asarray(cane_dulk_temperature_k(f), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio_db = 10.0 * np.log10(primary / candidate)
    return {
        "f_mhz": f.tolist(),
        "primary": "Franklin and MacLean (1969), Alouette in-situ",
        "primary_k": primary.tolist(),
        "primary_uncertainty_db": FRANKLIN_UNCERTAINTY_DB,
        "extrapolated": is_extrapolated(f).tolist(),
        "measured_range_mhz": [MEASURED_FLOOR_MHZ, MEASURED_CEILING_MHZ],
        "candidate": "Cane (1979) / Dulk et al. (2001) galactic background",
        "candidate_k": candidate.tolist(),
        "ratio_db": ratio_db.tolist(),
        "open_item": (
            "C1 is resolved only where a measured reference exists.  Below "
            f"{MEASURED_FLOOR_MHZ} MHz the primary reference is extrapolated and the "
            "candidate has never been checked against an in-situ topside "
            "measurement in that decade.  WP12."
        ),
    }
