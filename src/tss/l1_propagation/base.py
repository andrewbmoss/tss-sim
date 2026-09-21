"""Layer: L1
WP: WP2
Claims: C1, C2, C3
Status: implemented
Requirements: -

Shared helpers for the L1 propagation tiers.

Every L1 tier needs the same handful of operations on an
:class:`~tss.interfaces.IonosphereState`: sample the density at an arbitrary
altitude, get a collision frequency whether or not the L0 backend supplied one,
find the F2 peak, and report itself for the run manifest.  Those live here so
that the analytic tier and the ray-traced tier agree about them, rather than
each carrying its own copy that can drift.

There is no physics in this module.  The absorption integral is in
:mod:`tss.l1_propagation.absorption`, the magnetoionic algebra is in
:mod:`tss.l1_propagation.magnetoionic`, and the echo assembly is in the tier
modules.

One convention that is worth stating once.  A tier reports an altitude as a
float in metres, and reports the two physical states that are not altitudes as
strings: ``"evanescent"`` when the sounding frequency is below the local plasma
frequency at the spacecraft, and ``"penetrating"`` when it is above the F2
critical frequency.  Those are results, not errors, and contract rule 3 in
:mod:`tss.interfaces` says they must not raise.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Literal

import numpy as np

from tss.constants import PLASMA_CONST, plasma_frequency
from tss.interfaces import IonosphereState

__all__ = [
    "PropagationBase",
    "ReflectionState",
    "EVANESCENT",
    "PENETRATING",
    "density_at",
    "collision_frequency_at",
    "plasma_frequency_at",
    "peak_altitude_m",
    "critical_frequency_hz",
]

#: Sentinel returned instead of an altitude when the wave cannot leave the
#: antenna because the local plasma frequency already exceeds the sounding
#: frequency.
EVANESCENT = "evanescent"

#: Sentinel returned instead of an altitude when the wave passes through the F2
#: peak and there is no topside echo.
PENETRATING = "penetrating"

ReflectionState = Literal["evanescent", "penetrating"]


def _as_state_arrays(state: IonosphereState) -> tuple[np.ndarray, np.ndarray]:
    """Return the altitude and density grids as read-only float arrays."""
    alt = np.asarray(state.alt_m, dtype=float)
    ne = np.asarray(state.ne, dtype=float)
    if alt.ndim != 1 or ne.shape != alt.shape:
        raise ValueError(
            "IonosphereState.alt_m and .ne must be parallel one-dimensional arrays; "
            f"got {alt.shape} and {ne.shape}"
        )
    if alt.size < 2:
        raise ValueError("IonosphereState needs at least two altitude samples")
    if not np.all(np.diff(alt) > 0.0):
        raise ValueError("IonosphereState.alt_m must be strictly increasing")
    return alt, ne


def density_at(
    state: IonosphereState,
    alt_m: np.ndarray | float,
    *,
    interp: str = "linear",
) -> np.ndarray:
    """Electron density in m^-3 at arbitrary altitudes.

    ``interp`` selects how the tabulated profile is read between its nodes.

    ``"linear"``
        Straight linear interpolation of Ne against altitude.  This reproduces
        :class:`TabulatedProfile` in the author's link budget module exactly,
        which is why it is the default: the absorption numbers this package
        produces are then directly comparable with the tested implementation
        they were taken from.
    ``"log"``
        Linear interpolation of log Ne against altitude, which is closer to the
        exponential shape of a real topside and is the better choice on a coarse
        grid.  It is not the default only because switching it on would move the
        absorption numbers away from the reference implementation.

    Outside the tabulated range the value is held at the nearest endpoint, which
    is what :func:`numpy.interp` does.  That is a deliberate flat extrapolation
    rather than a physical model, so a caller that integrates above the top of
    the grid is integrating through a constant density and should know it.
    """
    alt, ne = _as_state_arrays(state)
    query = np.asarray(alt_m, dtype=float)
    if interp == "linear":
        return np.interp(query, alt, ne)
    if interp == "log":
        floor = np.maximum(ne, np.finfo(float).tiny)
        return np.exp(np.interp(query, alt, np.log(floor)))
    raise ValueError(f"unknown interp {interp!r}, expected 'linear' or 'log'")


def plasma_frequency_at(
    state: IonosphereState,
    alt_m: np.ndarray | float,
    *,
    interp: str = "linear",
) -> np.ndarray:
    """Local electron plasma frequency in Hz at arbitrary altitudes."""
    return plasma_frequency(density_at(state, alt_m, interp=interp))


def collision_frequency_at(
    state: IonosphereState,
    alt_m: np.ndarray | float,
    *,
    interp: str = "linear",
) -> np.ndarray:
    """Electron collision frequency in s^-1 at arbitrary altitudes.

    ``IonosphereState.nu`` is used when the L0 backend supplied it.  Otherwise
    the Coulomb form of :func:`tss.l1_propagation.absorption.collision_frequency`
    is evaluated from ``ne`` and ``te``.  A state that carries neither ``nu``
    nor ``te`` cannot support an absorption calculation, and this raises rather
    than inventing a number, because a silently assumed collision frequency
    would put an unsourced factor straight into every absorption result.
    """
    from tss.l1_propagation.absorption import collision_frequency

    alt, _ = _as_state_arrays(state)
    query = np.asarray(alt_m, dtype=float)

    if state.nu is not None:
        nu_grid = np.asarray(state.nu, dtype=float)
        if nu_grid.shape != alt.shape:
            raise ValueError("IonosphereState.nu must be parallel to alt_m")
        return np.interp(query, alt, nu_grid)

    if state.te is not None:
        te_grid = np.asarray(state.te, dtype=float)
        if te_grid.shape != alt.shape:
            raise ValueError("IonosphereState.te must be parallel to alt_m")
        te = np.interp(query, alt, te_grid)
        ne = density_at(state, query, interp=interp)
        return collision_frequency(ne, te)

    raise ValueError(
        "IonosphereState carries neither nu nor te, so the collision frequency "
        "cannot be established.  Absorption needs one of them; see the note on "
        "IonosphereState.nu in tss.interfaces."
    )


def peak_altitude_m(state: IonosphereState) -> float:
    """Altitude of the density maximum on the tabulated grid, metres.

    This is hmF2 for a state whose grid reaches the peak.  For a state that is
    topside only, it is the bottom of the grid, and the reflection search then
    treats the bottom of the grid as the floor.  That is the honest answer: a
    profile that does not contain the peak cannot say where the peak is.
    """
    alt, ne = _as_state_arrays(state)
    return float(alt[int(np.argmax(ne))])


def critical_frequency_hz(state: IonosphereState) -> float:
    """foF2 from the tabulated maximum density, Hz.

    Read off the grid rather than from ``state.extra["nmf2"]``, because
    contract rule in :mod:`tss.interfaces` says no tier may depend on an
    ``extra`` key being present.  A grid that undersamples the peak will
    understate foF2 by the amount the grid misses, which is a resolution
    question for L0 and not something L1 can repair.
    """
    _, ne = _as_state_arrays(state)
    return float(PLASMA_CONST * np.sqrt(float(np.max(ne))))


class PropagationBase:
    """Common base for L1 tiers.

    Holds the interpolation choice and the description machinery.  A subclass
    registers itself with :func:`tss.registry.register`, which sets ``name`` and
    ``layer`` on the class, and implements ``echoes`` to satisfy
    :class:`tss.interfaces.L1Propagation`.

    Subclasses add their own settings to :meth:`describe` by overriding
    :meth:`_settings`.  Everything a setting changes about a result belongs in
    there, because the run manifest records that mapping verbatim and a setting
    that is not in it is a setting that cannot be traced from a result.
    """

    #: Set by the registry decorator.  Declared here so that a tier that is
    #: instantiated without going through the registry still has the attribute.
    name: str = "base"
    layer: str = "l1"

    def __init__(self, *, interp: str = "linear") -> None:
        if interp not in ("linear", "log"):
            raise ValueError(f"unknown interp {interp!r}, expected 'linear' or 'log'")
        self.interp = interp

    # -- profile access -----------------------------------------------------

    def density_at(self, state: IonosphereState, alt_m: np.ndarray | float) -> np.ndarray:
        return density_at(state, alt_m, interp=self.interp)

    def plasma_frequency_at(
        self, state: IonosphereState, alt_m: np.ndarray | float
    ) -> np.ndarray:
        return plasma_frequency_at(state, alt_m, interp=self.interp)

    def collision_frequency_at(
        self, state: IonosphereState, alt_m: np.ndarray | float
    ) -> np.ndarray:
        return collision_frequency_at(state, alt_m, interp=self.interp)

    @staticmethod
    def peak_altitude_m(state: IonosphereState) -> float:
        return peak_altitude_m(state)

    @staticmethod
    def critical_frequency_hz(state: IonosphereState) -> float:
        return critical_frequency_hz(state)

    # -- manifest -----------------------------------------------------------

    def _settings(self) -> dict[str, Any]:
        """Settings this tier wants recorded.  Overridden by subclasses."""
        return {"interp": self.interp}

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable identity and settings, for the run manifest."""
        record: dict[str, Any] = {
            "layer": "l1",
            "name": getattr(self, "name", type(self).__name__),
            "class": type(self).__name__,
            "module": type(self).__module__,
        }
        record.update(self._settings())
        return record
