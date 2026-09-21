"""Separating Doppler from range using spacecraft motion.  THE principal contribution.

Layer:        study
WP:           WP5, WP6
Claims:       C3
Status:       stub
Requirements: INS-09, INS-15, INS-26, SCI-23

This is the study that resolves the dissertation's principal contribution, and it
is unstarted.  Its absence is the open risk in the plan, so the design is written
down here before the code exists.

The problem.  In a linear-sweep FMCW radar a Doppler shift displaces the beat
frequency, and the resulting apparent-range bias dd = v*f/K_f is independent of
range and grows linearly with sounding frequency.  At the top of the band it is
several range cells, which corrupts a retrieved profile rather than merely
degrading it.  The established remedy needs four antennas in two orthogonal pairs
and a dual-channel receiver, which a nanosatellite cannot carry.

The proposal.  The spacecraft's own motion does that work through processing
instead of aperture, because the velocity is large and known in advance to a
precision far exceeding what the measurement requires.  Three methods, in the
order they should be attempted:

1.  Up-sweep and down-sweep pairing.  The bias changes sign with sweep
    direction, so the half-sum of a matched pair recovers unbiased range and the
    half-difference recovers Doppler.  The cancellation is exact only to first
    order; the residual second-order term is what INS-09 bounds at one range
    cell of INS-05.
2.  Motion-compensated coherent integration over the known ephemeris.
3.  Along-track aperture synthesis, as an exploratory extension.

What makes this more than routine application.  SCI-23 demands a coherent
interval of one to two seconds, which is long enough that the geometry changes
materially within a single velocity estimate.  The platform therefore cannot be
treated as quasi-stationary, which is the assumption every ground-based
treatment makes.

What a first result looks like.  A noise-free simulation with a known injected
Doppler, up/down pairing applied, and the recovered range compared against truth.
Success is a residual bias under one range cell across the band.  That result
needs no ionosphere model beyond a smooth profile and no noise model at all, so
it is cheap, and it is the single most valuable thing to produce next because it
is the only claim with no evidence of tractability behind it.

A bounded negative result would also be a contribution.  Establishing that a
single-aperture FMCW sounder cannot separate plasma drift from platform motion
below some stated velocity resolution would settle a question the literature has
not asked.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from tss.studies import StudyResult, register_study


@register_study(
    "doppler_separation",
    description="Recover unbiased range and Doppler from a single transmit and receive pair",
    work_packages=["WP5", "WP6"],
    claims=["C3"],
    status="stub",
    requirements=["INS-09", "INS-15", "INS-26", "SCI-23"],
    blocked_on="WP5",
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    raise NotImplementedError(
        "WP5 is unstarted.  This study needs tss.l3_signal.novel.sweep_pairing and "
        "tss.l3_signal.novel.motion_compensation, both of which are stubs.  See "
        "this module's docstring for the design and for what a first result looks "
        "like, and docs/work-packages.md for where WP5 sits on the critical path."
    )
