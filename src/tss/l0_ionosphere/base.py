"""Layer: L0
WP: WP1
Claims: C1, C2
Status: implemented
Requirements: INS-01

Shared machinery for the L0 ionosphere tiers.

Every L0 tier answers the same question, ``state_at(lat, lon, epoch, alt_m)``,
and every one of them has to do the same four chores on the way to answering it:
validate the altitude grid, put an electron temperature on it when the model does
not supply one, derive a collision frequency because L1 absorption cannot run
without one, and attach a geomagnetic field.  Those live here so that the
analytic stand-in, the Chapman tier and the three IRI routes agree about them
rather than each carrying a copy that can drift.

There is no ionospheric model in this module.  The Chapman parameterisation is
in :mod:`tss.l0_ionosphere.chapman` and the migrated backends are in
:mod:`tss.l0_ionosphere.backends`.

Two conversions matter enough to name.

**Legacy profile objects to IonosphereState.**  The author's link budget module
carries ``IonosphereProfile`` and ``TabulatedProfile``, which are callables of
altitude rather than sampled arrays.  :func:`state_from_profile` samples one of
those onto a grid and returns the interface object.  It is duck-typed on ``ne``
and ``te`` deliberately, so that anything with those two methods can be lifted
into the layered package without being rewritten first.  That is the migration
path for the 37-test link budget module and for the hand-scaled ISIS-1 profiles.

**Collision frequency.**  L0 has to supply ``nu`` because the absorption
integral needs it and because a tier that leaves it ``None`` silently degrades
every absorption number downstream.  The implementation is imported from
:mod:`tss.l1_propagation.absorption` rather than copied.  That import runs from
L0 into L1, which is the wrong direction, and it is the lesser evil: the
alternative is two implementations of the Schunk and Nagy electron-ion Coulomb
form that can disagree.  It is the same trade recorded in the "known seams"
section of ``docs/architecture.md`` for the link budget module, and it comes out
in the same wash when WP4 splits the physics along layer lines.
"""

from __future__ import annotations

import datetime as dt
import math
from collections.abc import Mapping
from typing import Any, Protocol

import numpy as np

from tss.constants import PLASMA_CONST, plasma_frequency
from tss.interfaces import IonosphereState
from tss.l0_ionosphere import geomag

# Seam: the collision-frequency form lives with the absorption integral that
# consumes it.  See the module docstring and docs/architecture.md.
from tss.l1_propagation.absorption import collision_frequency  # noqa: F401

__all__ = [
    "IonosphereBase",
    "ProfileLike",
    "collision_frequency",
    "default_te_k",
    "parse_epoch",
    "state_from_profile",
    "validate_alt_grid",
    "R12_TO_F107_NOTE",
    "r12_to_f107",
    "IRI_F107_COEFFS",
    "iri_f107_from_r12",
    "iri_r12_from_f107",
]

R12_TO_F107_NOTE = (
    "Davies, Ionospheric Radio, 1990, p442.  This is the wrapper's own "
    "conversion and is reported so that a result records the solar index it was "
    "actually run at."
)


def r12_to_f107(r12: float | np.ndarray) -> np.ndarray | float:
    """The wrapper's own conversion (Davies, Ionospheric Radio, 1990, p442).

    Carried across verbatim from ``MinFreqStudy/iri_adapter.py``.  It is used to
    report the solar index a run was made at, and by the MATLAB route, which
    passes a user R12 rather than an epoch.
    """
    return 63.7 + 0.728 * r12 + 0.00089 * r12**2



#: IRI's own relation between the 12-month sunspot number and F10.7,
#: ``F10.7 = 63.75 + R12 (0.728 + 0.00089 R12)``.  The PHaRLAP ``iri2020`` mex
#: uses it to turn a user R12 into the F10.7 IRI runs on, and PyIRI uses it the
#: other way.  It differs from :func:`r12_to_f107` only in the constant (63.75
#: against Davies' 63.7, 0.05 SFU); the ``iri`` tier converts with this one so
#: that a stated F10.7 reaches the Fortran code unchanged.
IRI_F107_COEFFS: tuple[float, float, float] = (63.75, 0.728, 0.00089)


def iri_f107_from_r12(r12: float | np.ndarray) -> np.ndarray | float:
    """F10.7, SFU, for a 12-month sunspot number, by IRI's own relation."""
    a, b, c = IRI_F107_COEFFS
    return a + b * r12 + c * r12**2


def iri_r12_from_f107(f107: float | np.ndarray) -> np.ndarray | float:
    """The R12 that :func:`iri_f107_from_r12` maps to ``f107``: its positive root.

    Below 63.75 SFU the root is negative, which no caller should accept; it is
    returned rather than clipped so that the caller's range check sees it.
    """
    a, b, c = IRI_F107_COEFFS
    return (-b + np.sqrt(b * b - 4.0 * c * (a - np.asarray(f107, dtype=float)))) / (2.0 * c)

class ProfileLike(Protocol):
    """Anything with the legacy profile interface: ``ne(h)`` and ``te(h)``."""

    def ne(self, h: np.ndarray | float) -> np.ndarray: ...

    def te(self, h: np.ndarray | float) -> np.ndarray: ...


def validate_alt_grid(alt_m: np.ndarray | float) -> np.ndarray:
    """Return the altitude grid as a fresh, strictly increasing float array.

    A copy is made rather than a view, because the returned array goes into a
    frozen :class:`~tss.interfaces.IonosphereState` and contract rule 1 says a
    tier never hands back a reference to something a caller can still modify.
    """
    alt = np.array(np.atleast_1d(np.asarray(alt_m, dtype=float)), dtype=float, copy=True)
    if alt.ndim != 1:
        raise ValueError(f"alt_m must be one-dimensional, got shape {alt.shape}")
    if alt.size < 2:
        raise ValueError("alt_m needs at least two samples")
    if not np.all(np.diff(alt) > 0.0):
        raise ValueError("alt_m must be strictly increasing")
    if not np.all(np.isfinite(alt)):
        raise ValueError("alt_m must be finite")
    return alt


def parse_epoch(epoch: str | dt.datetime | None) -> dt.datetime:
    """ISO 8601 UTC text, or a datetime, to a naive UTC :class:`datetime`.

    The backends migrated from ``iri_adapter.py`` take datetimes, because IRI
    takes a calendar date and a UT hour.  The interface takes ISO text, because
    a manifest that stores a timestamp as text is readable years later and one
    that stores a pickled datetime is not.  This is the one place the two meet.
    """
    if isinstance(epoch, dt.datetime):
        return epoch.replace(tzinfo=None) if epoch.tzinfo is not None else epoch
    text = str(epoch or "").strip()
    if not text:
        raise ValueError("epoch is required; pass an ISO 8601 UTC timestamp")
    if text.endswith("Z"):
        text = text[:-1]
    try:
        parsed = dt.datetime.fromisoformat(text)
    except ValueError:
        for fmt in ("%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M", "%Y-%m-%d"):
            try:
                parsed = dt.datetime.strptime(text, fmt)
                break
            except ValueError:
                continue
        else:
            raise ValueError(f"could not read {epoch!r} as an ISO 8601 timestamp") from None
    return parsed.replace(tzinfo=None) if parsed.tzinfo is not None else parsed


def default_te_k(
    alt_m: np.ndarray,
    hmf2_m: float = 300.0e3,
    te_peak_k: float = 1500.0,
    te_top_k: float = 3000.0,
) -> np.ndarray:
    """Electron temperature on a grid when the model does not supply one.

    This is the linear ramp used by ``IonosphereProfile.te`` in the author's link
    budget module, carried across unchanged so that the absorption numbers this
    package produces stay comparable with the tested implementation: 1500 K at
    the F2 peak rising to 3000 K at 1000 km, clipped to the range 0 to 1.5 in
    the normalised variable so it neither runs backwards below the peak nor
    grows without bound above 1000 km.

    It is a stand-in, not a model.  Absorption goes as ``Te**-1.5`` through the
    Coulomb collision frequency, so a factor of two error in Te is about 4.5 dB
    on a term that is itself a few decibels.  A tier with a real Te, which means
    IRI, should use it and does.
    """
    alt = np.asarray(alt_m, dtype=float)
    span = max(1000.0e3 - float(hmf2_m), 1.0)
    frac = np.clip((alt - float(hmf2_m)) / span, 0.0, 1.5)
    return float(te_peak_k) + (float(te_top_k) - float(te_peak_k)) * frac


def state_from_profile(
    profile: ProfileLike,
    alt_m: np.ndarray,
    *,
    lat_deg: float = 0.0,
    lon_deg: float = 0.0,
    epoch: str = "",
    source: str = "",
    b_field_t: np.ndarray | None = None,
    extra: Mapping[str, Any] | None = None,
) -> IonosphereState:
    """Sample a legacy ``IonosphereProfile`` or ``TabulatedProfile`` onto a grid.

    The profile is asked for density and temperature and nothing else, which is
    why the protocol is two methods wide.  The collision frequency is derived
    here rather than asked for, because neither legacy class carries one.

    The profile object is not modified and nothing about it is retained in the
    returned state beyond the sampled arrays, so the state is safe to hold after
    the profile has been thrown away.
    """
    alt = validate_alt_grid(alt_m)
    ne = np.asarray(profile.ne(alt), dtype=float).reshape(-1)
    if ne.shape != alt.shape:
        raise ValueError(
            f"profile.ne returned shape {ne.shape} for a grid of shape {alt.shape}"
        )
    te = np.asarray(profile.te(alt), dtype=float).reshape(-1)
    nu = np.asarray(collision_frequency(ne, te), dtype=float)

    record: dict[str, Any] = dict(extra or {})
    record.setdefault("nmf2", float(np.max(ne)))
    record.setdefault("hmf2_km", float(alt[int(np.argmax(ne))] / 1.0e3))
    record.setdefault("fof2_hz", float(plasma_frequency(np.max(ne))))

    return IonosphereState(
        alt_m=alt,
        ne=ne,
        nu=nu,
        te=te,
        b_field_t=None if b_field_t is None else np.array(b_field_t, dtype=float, copy=True),
        lat_deg=float(lat_deg),
        lon_deg=float(lon_deg),
        epoch=str(epoch),
        source=str(source),
        extra=record,
    )


class IonosphereBase:
    """Common base for L0 tiers.

    Holds the geomagnetic-field choice and the description machinery.  A
    subclass registers itself with :func:`tss.registry.register`, which sets
    ``name`` and ``layer`` on the class, and implements ``state_at`` to satisfy
    :class:`tss.interfaces.L0Ionosphere`.

    Subclasses add their own settings to :meth:`describe` by overriding
    :meth:`_settings`.  Anything that can change a result belongs in there,
    because the run manifest records that mapping verbatim and a setting that is
    not in it is a setting that cannot be traced back from a result.
    """

    #: Set by the registry decorator.  Declared here so that a tier built
    #: directly, as the tests do, still has the attribute.
    name: str = "base"
    layer: str = "l0"

    #: Free text naming the model and its options, filled in by the subclass and
    #: copied into ``IonosphereState.source``.
    label: str = ""

    #: Longer provenance note, carried into ``describe()``.
    provenance: str = ""

    def __init__(self, *, geomagnetic_model: str = "auto") -> None:
        if geomagnetic_model not in ("auto", "igrf", "dipole", "none"):
            raise ValueError(
                f"unknown geomagnetic_model {geomagnetic_model!r}; expected "
                "auto, igrf, dipole or none"
            )
        self.geomagnetic_model = geomagnetic_model
        self._geomag_provenance: dict[str, Any] = {"model": "not evaluated"}

    # -- helpers a subclass uses -------------------------------------------

    def magnetic_field(
        self, lat_deg: float, lon_deg: float, alt_m: np.ndarray, epoch: str
    ) -> np.ndarray | None:
        """Field for the ``b_field_t`` slot, or ``None`` when it is switched off.

        The provenance of whichever route answered is stored for
        :meth:`describe`, so a result records not only that it had a field but
        where the field came from.
        """
        if self.geomagnetic_model == "none":
            self._geomag_provenance = {"model": "none", "reason": "switched off by configuration"}
            return None
        b, provenance = geomag.field_ned_t(
            lat_deg, lon_deg, alt_m, epoch, prefer=self.geomagnetic_model
        )
        self._geomag_provenance = dict(provenance)
        return b

    @staticmethod
    def collision_frequency(ne: np.ndarray, te: np.ndarray) -> np.ndarray:
        """Electron collision frequency, from the L1 absorption module."""
        return collision_frequency(ne, te)

    @staticmethod
    def peak_diagnostics(alt_m: np.ndarray, ne: np.ndarray) -> dict[str, float]:
        """``nmf2``, ``hmf2_km`` and ``fof2_hz`` read off a sampled profile.

        Reading the peak off a sampled grid is less accurate than taking it from
        a model that knows its own peak, which is why the IRI routes override
        these with the values IRI reports.  On a 5 km grid the height is out by
        up to half a step and the density by the curvature over that step, which
        is small for a Chapman layer and is not small for a structured profile.
        """
        alt = np.asarray(alt_m, dtype=float)
        dens = np.asarray(ne, dtype=float)
        finite = np.isfinite(dens)
        if not np.any(finite):
            return {"nmf2": float("nan"), "hmf2_km": float("nan"), "fof2_hz": float("nan")}
        idx = int(np.nanargmax(np.where(finite, dens, -np.inf)))
        peak = float(dens[idx])
        return {
            "nmf2": peak,
            "hmf2_km": float(alt[idx] / 1.0e3),
            "fof2_hz": float(PLASMA_CONST * math.sqrt(max(peak, 0.0))),
        }

    # -- manifest -----------------------------------------------------------

    def _settings(self) -> dict[str, Any]:
        """Settings this tier wants recorded.  Overridden by subclasses."""
        return {}

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable identity and settings, for the run manifest."""
        record: dict[str, Any] = {
            "layer": "l0",
            "name": getattr(self, "name", type(self).__name__),
            "class": type(self).__name__,
            "module": type(self).__module__,
            "label": self.label,
            "provenance": self.provenance,
            "geomagnetic_model": self.geomagnetic_model,
            "geomagnetic_provenance": dict(self._geomag_provenance),
        }
        record.update(self._settings())
        return record

    # -- contract -----------------------------------------------------------

    def state_at(
        self,
        lat_deg: float,
        lon_deg: float,
        epoch: str,
        alt_m: np.ndarray,
    ) -> IonosphereState:
        raise NotImplementedError(
            f"{type(self).__name__} does not implement state_at"
        )
