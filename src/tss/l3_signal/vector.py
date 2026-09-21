"""Layer: L3
WP: WP11
Claims: C5
Status: stub
Requirements: INS-08, INS-18

Three-channel vector receive: polarization-based O and X separation, and
arrival-angle estimation.  Not implemented.  Goal level, WP11.

Why it is wanted
----------------
Two requirements sit on this module and neither can be met by a scalar
receiver.

INS-08 asks that the ordinary and extraordinary traces be detected and
labelled.  They are physically distinguishable: the X trace reflects at the
cutoff ``f_x = f_ce/2 + sqrt((f_ce/2)^2 + f_p^2)`` of
:func:`tss.constants.x_mode_cutoff`, always above the O cutoff, and the two
modes arrive with opposite senses of circular polarization about the
geomagnetic field.  A single linear channel sees their sum and has no basis on
which to split it.  Three orthogonal channels measure the incident field
vector, and the projection onto the two characteristic polarizations separates
them.  With the modes labelled, the 150 kHz cap in the INS-02 step law can be
dropped, which is the A2 observing mode named in the requirement note, so the
benefit is not only that INS-08 is met but that the frequency stepping gets
cheaper.

INS-18 asks for 2 degrees of channel phase match, which exists to support the
2 degree arrival-angle goal of SCI-20.  Arrival angle is what would let an
off-vertical return from a depletion wall be recognised as off-vertical rather
than being inverted as though it came from below.  That connects this module to
:mod:`tss.l1_propagation.offvertical` and to :mod:`tss.l4_retrieval.depletion`,
and it is why all three are goal level together: none of them is useful without
the other two.

Intended approach
-----------------
For separation, form the Stokes parameters or the coherency matrix of the
three-channel de-chirped signal per range cell, and project onto the O and X
characteristic polarizations predicted from the local geomagnetic field
direction and the sounding frequency.  The projection is not onto pure circular
states except in the quasi-longitudinal limit; away from that limit the
characteristic polarizations are elliptical and depend on the angle between the
ray and the field, so the projection basis is a function of geometry and not a
constant.

For arrival angle, an interferometric or vector-sensor estimate over the three
channels, with the achievable accuracy set by the channel phase match of
INS-18 and by the signal to noise ratio per range cell.

Open questions
--------------
1.  Whether three short, electrically small antennas on a 6U bus are
    sufficiently orthogonal in practice.  Mutual coupling between elements that
    are a small fraction of a wavelength apart, on a conducting bus that is
    itself a small fraction of a wavelength, is not a small correction.  This
    is an L2 question that decides whether the L3 algorithm has anything to
    work with, and it belongs to :mod:`tss.l2_instrument` and the WP9 bench
    work.
2.  What the separation is worth near the cusp.  That is where O and X are
    closest in range and where labelling them matters most for the inversion,
    and it is also where both returns are spread and weak.
3.  Whether the arrival-angle accuracy the geometry allows is good enough to be
    useful at all, given that the expected off-vertical angles from a depletion
    wall are themselves only a few degrees.  If the estimate is 2 degrees and
    the signal is 3 degrees, the measurement exists but does not discriminate.

Raises
------
NotImplementedError
    On every call.  WP11.
"""

from __future__ import annotations

import numpy as np

from tss.interfaces import Ionogram, ReceivedSignal

__all__ = [
    "characteristic_polarizations",
    "separate_modes",
    "arrival_angle_deg",
    "label_ionogram_modes",
]

_MESSAGE = (
    "WP11: vector receive processing is not implemented.  Three-channel O and X "
    "separation and arrival-angle estimation are goal level and depend on the "
    "three-channel L2 receive chain, which does not exist either.  See "
    "docs/architecture.md section 'L3 Signal processing' and docs/work-packages.md "
    "row WP11.  Requirements INS-08 and INS-18, claim C5."
)


def characteristic_polarizations(
    f_hz: float | np.ndarray,
    b_field_t: np.ndarray,
    u_los: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Unit Jones vectors of the O and X characteristic polarizations.

    Not implemented.  Intended to return the two complex two-vectors, in a
    frame transverse to ``u_los``, that the ordinary and extraordinary waves
    take at sounding frequency ``f_hz`` for the geomagnetic field
    ``b_field_t``.  In the quasi-longitudinal limit these are the two circular
    states; in general they are elliptical and the ellipticity is what makes
    the separation basis geometry dependent rather than fixed.

    Raises
    ------
    NotImplementedError
        Always.  WP11.
    """
    del f_hz, b_field_t, u_los
    raise NotImplementedError(_MESSAGE)


def separate_modes(
    signal: ReceivedSignal,
    b_field_t: np.ndarray,
    u_los: np.ndarray,
) -> tuple[ReceivedSignal, ReceivedSignal]:
    """Split a three-channel record into its O and X parts.  Not implemented.

    Intended to return two single-channel :class:`~tss.interfaces.ReceivedSignal`
    objects that the ordinary baseline chain can then process independently,
    so that mode separation is a pre-processing step rather than a change to
    the ionogram former.  Keeping it that way is what lets the same trace
    extractor be used on both and lets the two be compared.

    Raises
    ------
    NotImplementedError
        Always.  WP11.
    """
    del signal, b_field_t, u_los
    raise NotImplementedError(_MESSAGE)


def arrival_angle_deg(
    signal: ReceivedSignal,
    range_m: np.ndarray,
) -> np.ndarray:
    """Arrival angle per range cell from the three channels.  Not implemented.

    Intended to populate :attr:`tss.interfaces.Ionogram.aoa_deg`.  The accuracy
    is bounded by the channel phase match of INS-18 and by the per-cell signal
    to noise ratio, and the estimator has to report its own uncertainty or the
    downstream depletion test of WP10 has nothing to threshold on.

    Raises
    ------
    NotImplementedError
        Always.  WP11.
    """
    del signal, range_m
    raise NotImplementedError(_MESSAGE)


def label_ionogram_modes(
    ionogram: Ionogram,
    signal: ReceivedSignal,
    b_field_t: np.ndarray,
    u_los: np.ndarray,
) -> Ionogram:
    """Return an ionogram with O and X columns labelled.  Not implemented.

    Intended to replace the single ``"unresolved"`` column that
    :class:`tss.l3_signal.ionogram.BaselineSignalProcessing` produces with the
    two labelled columns the :class:`~tss.interfaces.Ionogram` default assumes,
    and to set ``flags["mode_discrimination"]`` True so that L5 scores INS-08
    as met on evidence rather than on assertion.

    Raises
    ------
    NotImplementedError
        Always.  WP11.
    """
    del ionogram, signal, b_field_t, u_los
    raise NotImplementedError(_MESSAGE)
