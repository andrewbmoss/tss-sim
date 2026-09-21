"""The sub-megahertz sky noise reference.

Layer:        study
WP:           WP12
Claims:       C1
Status:       partial
Requirements: INS-01, INS-10

Franklin and MacLean measured the galactic background in orbit from 0.6 to
12 MHz.  The lower edge of the sounding band is 0.33 MHz, so across the bottom
0.27 MHz of the band the noise reference this work relies on is an
extrapolation rather than a measurement.  That is the open part of C1.

The risk runs in the unintuitive direction.  A sounder is usually helped by a
quieter sky, but here the receive-efficiency cancellation depends on the antenna
temperature staying sky-dominated, and that condition needs eta*T_sky to stay
well above about 580 K.  A sky that is lower than extrapolated therefore breaks
the cancellation and costs the small-antenna argument its margin.  So the useful
output of this study is a lower bound on the sky, not a best estimate.

This study compares the candidate extensions against each other and reports the
spread as the uncertainty.  It does not pick a winner, because picking one is
what WP12 is for and the evidence to do it is not in yet.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study


@register_study(
    "sky_noise",
    description="Bound the sky noise reference below 0.6 MHz for WP12",
    work_packages=["WP12"],
    claims=["C1"],
    status="partial",
    requirements=["INS-01", "INS-10"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    from tss.l2_instrument import linkbudget as lb
    from tss.l2_instrument import noise

    f_mhz = np.linspace(0.2, 2.0, 91)

    cane_dulk = np.asarray(noise.cane_dulk_temperature_k(f_mhz), dtype=float)
    measured = np.asarray(lb.sky_temperature(f_mhz, "franklin"), dtype=float)
    extrapolated = np.asarray(noise.is_extrapolated(f_mhz), dtype=bool)

    # A plain power-law continuation of the measured curve, anchored on the two
    # lowest measured points.  This is the least defensible of the candidates
    # and is included so the spread is honest about how wide it is.
    anchor = f_mhz >= 0.6
    if anchor.sum() >= 2:
        f_a = f_mhz[anchor][:2]
        t_a = measured[anchor][:2]
        slope = np.log(t_a[1] / t_a[0]) / np.log(f_a[1] / f_a[0])
        power_law = t_a[0] * (f_mhz / f_a[0]) ** slope
    else:  # pragma: no cover - the grid always contains two anchor points
        slope = float("nan")
        power_law = np.full_like(f_mhz, np.nan)

    candidates = np.column_stack([f_mhz, measured, cane_dulk, power_law,
                                  extrapolated.astype(float)])
    path = out_dir / "sky_noise_candidates.csv"
    np.savetxt(
        path, candidates, delimiter=",",
        header="f_mhz,franklin_maclean_k,cane_dulk_k,power_law_k,is_extrapolated",
        comments="",
    )

    below = f_mhz < 0.6
    stack = np.vstack([cane_dulk[below], power_law[below]])
    lower = np.nanmin(stack, axis=0)
    upper = np.nanmax(stack, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        spread_db = 10.0 * np.log10(upper / lower)

    #: The condition the receive-efficiency cancellation needs.  Below this the
    #: antenna temperature stops being sky-dominated and the cancellation fails.
    sky_dominance_floor_k = 580.0
    worst = float(np.nanmin(lower)) if lower.size else float("nan")

    metrics = {
        "power_law_slope": float(slope),
        "spread_db_max_below_0p6mhz": float(np.nanmax(spread_db)) if spread_db.size else None,
        "lowest_candidate_sky_k_below_0p6mhz": worst,
        "sky_dominance_floor_k": sky_dominance_floor_k,
        "lowest_candidate_clears_floor": bool(worst > sky_dominance_floor_k),
        "candidates": ["cane_dulk", "power_law"],
        "settled": False,
    }

    summary = f"""
Three curves were compared from 0.2 to 2.0 MHz: the Franklin and MacLean in-orbit
measurement, the Cane and Dulk galactic background, and a power-law continuation
of the measurement anchored at 0.6 MHz with exponent {slope:.2f}.

Below 0.6 MHz, where the measurement does not reach, the candidates spread by up
to {float(np.nanmax(spread_db)):.1f} dB.  That spread is the honest uncertainty on
the noise reference at the bottom of the band, and it is reported here instead of
a best estimate.

The number that matters is the lower bound rather than the central value.  The
receive-efficiency cancellation holds only while the antenna temperature is
sky-dominated, which needs the sky to stay well above about {sky_dominance_floor_k:.0f} K.
The lowest candidate over the sub-megahertz range is {worst:.3g} K, which
{'clears' if worst > sky_dominance_floor_k else 'does not clear'} that floor.

This study bounds the answer.  It does not settle it.  WP12 settles it, and until
it does, any SNR quoted below 0.6 MHz carries this spread.
""".strip()

    return StudyResult(outputs=[path], summary=summary, metrics=metrics)
