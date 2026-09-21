"""Layer: L0
WP: WP1
Claims: C1, C2
Status: implemented
Requirements: INS-01, SCI-02

The alpha-Chapman topside, parameterised by the PLASMA scale height.

This is the default L0 tier named in ``config/default.yaml``, and it is the only
one that is guaranteed to work on a machine with no optional dependency
installed.  It uses numpy and nothing else.  Everything it needs is closed form,
so it is also the tier a limiting-case check can be written against: the profile
has an analytic peak, an analytic critical frequency and an analytic asymptotic
decay, and a result that disagrees with those is a bug in the layer being tested
rather than a question about the ionosphere.

The profile
-----------
Above the F2 peak the density is the alpha-Chapman form

    Ne(h) = NmF2 exp( 0.5 [ 1 - z - exp(-z) ] ),   z = (h - hmF2) / H_chapman

and far above the peak the bracket tends to ``-z``, so the density decays as
``exp(-z/2)``, that is with a scale height of ``2 H_chapman``.

WHICH SCALE HEIGHT
------------------
The constructor takes the PLASMA scale height ``Hp``, not the neutral scale
height, and internally sets ``H_chapman = Hp / 2``.  This is deliberate and it
is the single easiest thing in this module to get wrong, so it is stated here
rather than left in a comment.

``Hp`` is what the topside literature measures and what the science requirements
are written in terms of.  Meehan (2017) measures it as 150 +/- 30 km at night.
Stankov et al. (2003) show it running from 78 to 244 km for the same data
depending only on which profiler was used to extract it, which is itself part of
the case for a sounder that measures the profile rather than fitting one.  An
alpha-Chapman layer decays asymptotically with twice its Chapman scale height,
so a module parameterised by ``H_chapman`` and a requirement written in ``Hp``
differ by a factor of two, which is a factor of two in topside density at a
scale height above the peak and considerably more further up.  Taking ``Hp`` as
the argument removes the chance of that mismatch at the cost of one division.

Below the peak
--------------
A topside sounder does not see below the F2 peak, so what this tier does below
``hmF2`` is a convention and not a model.  The same alpha-Chapman expression is
evaluated there, which gives a bottomside that is symmetric in ``z`` about the
peak in the sense of the Chapman function and is roughly the right shape for a
single layer with no E or F1 ledge.  It is fine as a lower boundary for an
absorption integral that stops at the reflection height and it is not a
bottomside model.  Any result that depends on the shape below the peak needs
IRI, and the tier says so in its ``extra`` under ``bottomside_is_convention``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from tss.constants import plasma_frequency
from tss.interfaces import IonosphereState
from tss.l0_ionosphere.base import (
    IonosphereBase,
    collision_frequency,
    default_te_k,
    validate_alt_grid,
)
from tss.registry import register

__all__ = ["ChapmanIonosphere", "chapman_density", "chapman_scale_height_m"]


def chapman_scale_height_m(plasma_scale_height_m: float) -> float:
    """``H_chapman = Hp / 2``.  See the module docstring for why."""
    return float(plasma_scale_height_m) / 2.0


def chapman_density(
    alt_m: np.ndarray | float,
    nmf2: float,
    hmf2_m: float,
    plasma_scale_height_m: float,
) -> np.ndarray:
    """alpha-Chapman electron density in m^-3.

    Parameters
    ----------
    alt_m
        Altitudes, metres.
    nmf2
        Peak density, m^-3.
    hmf2_m
        Peak height, metres.
    plasma_scale_height_m
        The PLASMA scale height ``Hp``, metres.  Not the neutral scale height
        and not the Chapman scale height.

    The exponent is clipped before it is exponentiated.  Far below the peak the
    ``exp(-z)`` term overflows for ``z`` below about -700, and the density there
    is a convention rather than a result, so the clip changes nothing that
    matters and stops a warning that would otherwise appear on every call with a
    grid that reaches down to 90 km.
    """
    h = np.asarray(alt_m, dtype=float)
    hc = chapman_scale_height_m(plasma_scale_height_m)
    if hc <= 0.0:
        raise ValueError("plasma_scale_height_m must be positive")
    z = (h - float(hmf2_m)) / hc
    arg = 0.5 * (1.0 - z - np.exp(np.clip(-z, -700.0, 700.0)))
    return float(nmf2) * np.exp(np.clip(arg, -700.0, 700.0))


@register("l0", "chapman")
class ChapmanIonosphere(IonosphereBase):
    """Closed-form alpha-Chapman topside.  The default tier.

    Parameters
    ----------
    nmf2
        Peak density, m^-3.  The default of 5e11 gives foF2 = 6.35 MHz, which is
        the value the author's link budget module uses as its own default.
    hmf2_m
        Peak height, metres.
    plasma_scale_height_m
        The plasma scale height ``Hp``, metres.  See the module docstring.
    te_peak_k, te_top_k
        Electron temperature at the peak and at 1000 km, K.  Used to derive the
        collision frequency, which L1 absorption needs.
    geomagnetic_model
        Passed to :mod:`tss.l0_ionosphere.geomag`.  ``"auto"`` uses ppigrf when
        it is importable and the centred dipole when it is not, so this tier
        keeps working with nothing installed.

    The tier ignores ``lat_deg``, ``lon_deg`` and ``epoch`` for the density, and
    says so in ``extra["varies_with_location"]``.  A tier whose profile does not
    depend on where it is asked has to advertise that, otherwise a campaign that
    sweeps latitude will produce a flat line and the flatness will be read as a
    result.  The geomagnetic field does vary with location, so the O and X
    separation this tier produces is real even though the density is not.
    """

    def __init__(
        self,
        *,
        nmf2: float = 5.0e11,
        hmf2_m: float = 300.0e3,
        plasma_scale_height_m: float = 100.0e3,
        te_peak_k: float = 1500.0,
        te_top_k: float = 3000.0,
        geomagnetic_model: str = "auto",
    ) -> None:
        super().__init__(geomagnetic_model=geomagnetic_model)
        if nmf2 <= 0.0:
            raise ValueError("nmf2 must be positive")
        if plasma_scale_height_m <= 0.0:
            raise ValueError("plasma_scale_height_m must be positive")
        self.nmf2 = float(nmf2)
        self.hmf2_m = float(hmf2_m)
        self.plasma_scale_height_m = float(plasma_scale_height_m)
        self.te_peak_k = float(te_peak_k)
        self.te_top_k = float(te_top_k)
        self.label = (
            f"alpha-Chapman(NmF2={self.nmf2:.3e} m^-3, hmF2={self.hmf2_m / 1e3:.0f} km, "
            f"Hp={self.plasma_scale_height_m / 1e3:.0f} km)"
        )
        self.provenance = (
            "Closed-form alpha-Chapman topside parameterised by the plasma scale "
            "height.  Analytic, not a model of any particular ionosphere; use it "
            "for limiting cases and for the fast tier of a trade study."
        )

    # -- derived quantities --------------------------------------------------

    @property
    def chapman_scale_height_m(self) -> float:
        """``Hp / 2``, the scale height that appears in the Chapman exponent."""
        return chapman_scale_height_m(self.plasma_scale_height_m)

    def fof2_hz(self) -> float:
        """Critical frequency of the F2 layer, Hz."""
        return float(plasma_frequency(self.nmf2))

    def local_plasma_frequency_hz(self, alt_m: float) -> float:
        """Plasma frequency at one altitude, the lower bound on sounding there.

        This is the quantity the band-edge study of WP0 is a distribution of, and
        it is exposed here so a study can ask the tier for it directly rather
        than sampling a profile and reading it back off.
        """
        return float(
            plasma_frequency(
                chapman_density(
                    float(alt_m), self.nmf2, self.hmf2_m, self.plasma_scale_height_m
                )
            )
        )

    # -- contract ------------------------------------------------------------

    def state_at(
        self,
        lat_deg: float,
        lon_deg: float,
        epoch: str,
        alt_m: np.ndarray,
    ) -> IonosphereState:
        """Truth state on ``alt_m``.  The input grid is copied, never modified."""
        alt = validate_alt_grid(alt_m)
        ne = chapman_density(alt, self.nmf2, self.hmf2_m, self.plasma_scale_height_m)
        te = default_te_k(alt, self.hmf2_m, self.te_peak_k, self.te_top_k)
        nu = collision_frequency(ne, te)
        b = self.magnetic_field(lat_deg, lon_deg, alt, epoch)

        extra: dict[str, Any] = {
            "nmf2": self.nmf2,
            "hmf2_km": self.hmf2_m / 1.0e3,
            "fof2_hz": self.fof2_hz(),
            "plasma_scale_height_km": self.plasma_scale_height_m / 1.0e3,
            "chapman_scale_height_km": self.chapman_scale_height_m / 1.0e3,
            "varies_with_location": False,
            "bottomside_is_convention": True,
            "te_is_stand_in": True,
        }

        return IonosphereState(
            alt_m=alt,
            ne=np.asarray(ne, dtype=float),
            nu=np.asarray(nu, dtype=float),
            te=np.asarray(te, dtype=float),
            b_field_t=b,
            lat_deg=float(lat_deg),
            lon_deg=float(lon_deg),
            epoch=str(epoch),
            source=self.label,
            extra=extra,
        )

    def _settings(self) -> dict[str, Any]:
        return {
            "nmf2": self.nmf2,
            "hmf2_m": self.hmf2_m,
            "plasma_scale_height_m": self.plasma_scale_height_m,
            "chapman_scale_height_m": self.chapman_scale_height_m,
            "te_peak_k": self.te_peak_k,
            "te_top_k": self.te_top_k,
            "fof2_hz": self.fof2_hz(),
        }
