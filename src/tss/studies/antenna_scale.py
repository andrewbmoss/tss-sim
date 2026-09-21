"""The antenna scale trade.

Layer:        study
WP:           WP3
Claims:       C1
Status:       implemented
Requirements: INS-01, INS-10, INS-14

A resonant dipole at the lower band edge is 454 m, which no nanosatellite will
deploy.  This study prices every practical length instead.  Three results come
out of it.  The dominant loss is the matching network rather than the conductor.
The penalty falls on transmit alone, because receive antenna efficiency cancels
from the SNR while the antenna temperature is sky-dominated.  And sky dominance
is itself a condition that has to be checked rather than assumed, since it fails
once the efficiency is low enough that the physical temperature of the structure
dominates the noise at the terminals.

OPEN ITEM.  The author's reported sky-dominance margin for the 5 m antenna is
0.5 dB.  Recomputing it from the defaults in the original small_antenna.py gives
about +10 dB at 0.6 MHz.  The two do not reconcile without roughly three times
lower coil Q or about 5 dB more matching-network loss.  This study reports both
numbers with the assumption each rests on rather than choosing one, and the
discrepancy is flagged in the summary.  It needs the author to say which antenna
assumption the 0.5 dB came from.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

#: Practical deployed scales, metres.  1 m is what a 6U bus manages with no
#: deployment mechanism worth arguing about; 25 m is an ambitious boom.
LENGTHS_M = (0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0)


@register_study(
    "antenna_scale",
    description="Price every practical antenna length against the band edge",
    work_packages=["WP3"],
    claims=["C1"],
    status="implemented",
    requirements=["INS-01", "INS-10", "INS-14"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    from tss.l2_instrument import antenna as ant
    from tss.l2_instrument import linkbudget as lb

    f_edge_hz = float(config.get("sweep", {}).get("f_start_hz", 3.3e5))
    f_mhz = np.array([f_edge_hz / 1e6, 0.6, 1.0, 5.0, 12.0, 20.0])

    resonant_m = float(np.asarray(ant.resonant_length_m(f_edge_hz)).ravel()[0])

    rows: list[tuple[float, float, float, float, float]] = []
    for length in LENGTHS_M:
        for f in f_mhz:
            eff_db = float(
                np.asarray(ant.short_dipole_efficiency(f * 1e6, length)).ravel()[0]
            )
            t_sky = float(np.asarray(lb.sky_temperature(f, "franklin")).ravel()[0])
            margin = float(
                np.asarray(ant.sky_dominance_margin_db(f * 1e6, length)).ravel()[0]
            )
            rows.append((length, f, eff_db, t_sky, margin))

    path = out_dir / "antenna_scale.csv"
    np.savetxt(
        path, np.array(rows), delimiter=",",
        header="length_m,f_mhz,radiation_efficiency_db,t_sky_k,sky_dominance_margin_db",
        comments="",
    )

    five = [r for r in rows if r[0] == 5.0]
    margin_at_0p6 = next((r[4] for r in five if abs(r[1] - 0.6) < 1e-9), float("nan"))
    margin_at_edge = five[0][4] if five else float("nan")

    metrics = {
        "band_edge_mhz": f_edge_hz / 1e6,
        "resonant_half_wave_length_m": resonant_m,
        "lengths_priced_m": list(LENGTHS_M),
        "sky_dominance_margin_db_5m_at_0p6mhz_recomputed": margin_at_0p6,
        "sky_dominance_margin_db_5m_at_band_edge_recomputed": margin_at_edge,
        "sky_dominance_margin_db_5m_as_reported_by_author": 0.5,
        "open_item": "recomputed and reported sky-dominance margins disagree",
    }

    summary = f"""
A half-wave dipole resonant at the {f_edge_hz / 1e6:.2f} MHz band edge is
{resonant_m:.0f} m.  That is not deployable from a nanosatellite, so the trade is
over electrically short antennas and what their inefficiency costs.

Radiation efficiency was evaluated for {len(LENGTHS_M)} deployed lengths at
{f_mhz.size} frequencies.  The loss is dominated by the matching coil rather than
by the conductor, and it falls on transmit alone, because receive efficiency
cancels from the SNR while the antenna temperature is sky-dominated.

OPEN ITEM, carried rather than resolved.  The sky-dominance margin for the 5 m
antenna recomputes to {margin_at_0p6:.1f} dB at 0.6 MHz from the defaults in the
original small_antenna.py, against the 0.5 dB reported previously.  Reconciling
them needs roughly three times lower coil Q or about 5 dB more matching-network
loss.  Both figures are in the metrics above with the assumption each rests on.
Which one is right changes how much margin the small-antenna argument has, so it
should be settled before this number is cited.
""".strip()

    return StudyResult(outputs=[path], summary=summary, metrics=metrics)
