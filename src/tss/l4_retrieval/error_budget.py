"""Layer: L4
WP: WP6
Claims: C2, C3
Status: stub
Requirements: INS-02, INS-05, INS-09

Propagated uncertainty through the lamination inversion.  Not implemented.  WP6.

What it has to produce
-----------------------
:class:`~tss.interfaces.DensityProfile` carries ``ne_sigma`` and
``alt_sigma_m`` and :class:`tss.l4_retrieval.lamination.LaminationRetrieval`
leaves both ``None``.  Leaving them ``None`` is the tier saying it has not
computed them, which is correct and is also the reason the dissertation cannot
yet quote a retrieval accuracy.  Section 1.3 of the proposal reports per
configuration and an accuracy without an uncertainty is not a report.

What goes in
------------
Four independent sources, and they do not enter the same way.

1.  **Virtual range error.**  One range cell of INS-05, 10 km at threshold and
    5 km at goal, plus the leading-edge bias of
    :mod:`tss.l3_signal.trace_extract`, plus the residual Doppler bias of
    INS-09 where C3 has been applied and the full first-order bias where it has
    not.  This enters the right-hand side of the lower-triangular system.
2.  **Frequency error.**  The INS-02 step law sets how finely the frequency
    axis is sampled, and ``dNe/Ne = 2 df/f`` from
    :func:`tss.constants.density_precision_from_freq_step` turns a fractional
    frequency step into a fractional density step.  This enters the matrix, not
    the right-hand side, because the frequencies are the lamina boundaries.
3.  **Starting value error.**  ``f_0`` displaces every height below it.  The
    three ways of obtaining it, set out in :mod:`tss.l4_retrieval.base`, have
    different magnitudes, and the ``lowest_trace`` fallback is biased rather
    than merely uncertain.
4.  **Model error from the lamina law.**  Taking the plasma frequency linear in
    height inside each lamina is an approximation whose error depends on how
    curved the real profile is across a lamina.  It is not random and it does
    not average down with more laminae in any simple way.

Why it is not a per-point error bar
------------------------------------
The inversion is forward substitution down a lower-triangular system, so
``t_k`` depends on every ``t_j`` for ``j < k``.  An error in an early lamina
displaces everything beneath it and never washes out.  The propagated
covariance is therefore dense and strongly correlated down the profile, and
quoting the diagonal of it as an error bar would badly understate the
uncertainty in any quantity that integrates down the profile, which includes
total electron content and the height of any feature.

The intended route is to propagate the full covariance through the triangular
solve, ``C_t = L^-1 C_R L^-T`` with the matrix contributions handled by a
first-order expansion in the frequency errors, and to report both the diagonal
and the correlation length down the profile.  Reporting only the diagonal would
be the usual mistake.

The non-random part
-------------------
The leading-edge bias, the starting-value bias when the ``lowest_trace``
fallback was used, the no-field group index of
:mod:`tss.l4_retrieval.lamination`, and any uncorrected Doppler bias are all
biases and not variances.  They have to be carried separately and reported
separately, because adding a bias into a variance in quadrature understates it
when it is small and overstates the random part when it is large.  The error
budget of WP6 is therefore two tables and not one.

Raises
------
NotImplementedError
    On every call.  WP6.
"""

from __future__ import annotations

import numpy as np

from tss.interfaces import DensityProfile, Ionogram, PlatformState

__all__ = [
    "range_covariance_m2",
    "propagate_covariance",
    "bias_budget_m",
    "with_uncertainty",
]

_MESSAGE = (
    "WP6: the retrieval error budget is not implemented.  The lamination inversion "
    "leaves DensityProfile.ne_sigma and .alt_sigma_m as None, which is the tier "
    "saying it has not computed them.  The propagation is a dense correlated "
    "covariance down the profile, not a per-point error bar, and the biases "
    "(leading edge, starting value, no-field group index, uncorrected Doppler) are "
    "carried separately from the variances.  See the module docstring of "
    "tss.l4_retrieval.error_budget and docs/work-packages.md row WP6."
)


def range_covariance_m2(
    ionogram: Ionogram,
    *,
    range_resolution_m: float | None = None,
    include_doppler_bias: bool = True,
) -> np.ndarray:
    """Covariance of the measured virtual ranges.  Not implemented.

    Intended to combine the per-cell range uncertainty from INS-05 and the
    signal to noise ratio, the leading-edge bias, and the residual Doppler bias
    carried on ``Ionogram.residual_range_bias_m``.  When that field is ``None``
    the L3 tier did not estimate the Doppler, and the correct treatment is the
    full first-order bias of :func:`tss.constants.doppler_range_bias` rather
    than zero.

    Raises
    ------
    NotImplementedError
        Always.  WP6.
    """
    del ionogram, range_resolution_m, include_doppler_bias
    raise NotImplementedError(_MESSAGE)


def propagate_covariance(
    range_covariance: np.ndarray,
    freq_hz: np.ndarray,
    f_start_hz: float,
    *,
    f_start_sigma_hz: float = 0.0,
) -> np.ndarray:
    """Covariance of the lamina thicknesses.  Not implemented.

    Intended to evaluate ``C_t = L^-1 C_R L^-T`` with ``L`` the lamination
    matrix of :func:`tss.l4_retrieval.lamination.lamination_matrix`, plus the
    first-order contribution of the frequency errors, which enter the matrix
    rather than the right-hand side.  The full matrix is returned, not its
    diagonal; see the module docstring for why the diagonal alone would
    understate anything that integrates down the profile.

    Raises
    ------
    NotImplementedError
        Always.  WP6.
    """
    del range_covariance, freq_hz, f_start_hz, f_start_sigma_hz
    raise NotImplementedError(_MESSAGE)


def bias_budget_m(
    ionogram: Ionogram,
    profile: DensityProfile,
) -> dict[str, float]:
    """The non-random part of the height error, term by term.  Not implemented.

    Intended to return one entry per identified bias, in metres of height
    displacement, so that the budget can be read as a table rather than as a
    single number.  Kept separate from the variances on purpose.

    Raises
    ------
    NotImplementedError
        Always.  WP6.
    """
    del ionogram, profile
    raise NotImplementedError(_MESSAGE)


def with_uncertainty(
    profile: DensityProfile,
    ionogram: Ionogram,
    platform: PlatformState,
) -> DensityProfile:
    """Return the profile with ``ne_sigma`` and ``alt_sigma_m`` filled.  Not implemented.

    Intended to be composed after any L4 tier rather than built into one, so
    that the error budget is the same whichever inversion produced the profile.
    The input is not modified.

    Raises
    ------
    NotImplementedError
        Always.  WP6.
    """
    del profile, ionogram, platform
    raise NotImplementedError(_MESSAGE)
