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

THE 5 m ANTENNA.  With the Balanis short-dipole reactance (corrected 2026-10-05;
the earlier model used twice the terminal capacitance and was 3 dB optimistic)
the 5 m antenna's efficiency at 0.33 MHz is -46.1 dB.  Its sky-dominance margin
there is +0.5 dB against 2.6e7 K and -2.5 dB against 1.3e7 K, the bracket the
requirement set uses, so it is marginal to failing at the band edge.  The
earlier 'open item' set a 0.5 dB margin against a recomputed +10 dB; those were
the same model at different frequencies and skies (0.33 MHz on 1.3e7 K against
0.6 MHz on 2e7 K), not two antenna assumptions.
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
            eta, _ = ant.short_dipole_efficiency(f * 1e6, length)
            eta = float(np.asarray(eta).ravel()[0])
            eff_db = 10.0 * float(np.log10(eta))
            t_sky = float(np.asarray(lb.sky_temperature(f, "franklin")).ravel()[0])
            margin = float(np.asarray(ant.sky_dominance_margin_db(eta, t_sky)).ravel()[0])
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
    eta_edge = float(np.asarray(ant.short_dipole_efficiency(f_edge_hz, 5.0)[0]).ravel()[0])
    bracket = {
        "cane": float(np.asarray(ant.sky_dominance_margin_db(eta_edge, 2.6e7)).ravel()[0]),
        "low": float(np.asarray(ant.sky_dominance_margin_db(eta_edge, 1.3e7)).ravel()[0]),
    }

    metrics = {
        "band_edge_mhz": f_edge_hz / 1e6,
        "resonant_half_wave_length_m": resonant_m,
        "lengths_priced_m": list(LENGTHS_M),
        "radiation_efficiency_db_5m_at_band_edge": 10.0 * float(np.log10(eta_edge)),
        "sky_dominance_margin_db_5m_at_0p6mhz_franklin": margin_at_0p6,
        "sky_dominance_margin_db_5m_at_band_edge_franklin_extrapolated": margin_at_edge,
        "sky_dominance_margin_db_5m_at_band_edge_2p6e7k": bracket["cane"],
        "sky_dominance_margin_db_5m_at_band_edge_1p3e7k": bracket["low"],
    }

    summary = f"""
A half-wave dipole resonant at the {f_edge_hz / 1e6:.2f} MHz band edge is
{resonant_m:.0f} m.  That is not deployable from a nanosatellite, so the trade is
over electrically short antennas and what their inefficiency costs.

Radiation efficiency was evaluated for {len(LENGTHS_M)} deployed lengths at
{f_mhz.size} frequencies.  The loss is dominated by the matching coil rather than
by the conductor, and it falls on transmit alone, because receive efficiency
cancels from the SNR while the antenna temperature is sky-dominated.

The 5 m antenna is where sky dominance runs out.  Its efficiency at the band
edge is {10.0 * np.log10(eta_edge):.1f} dB, and its sky-dominance margin there is
{bracket["cane"]:+.1f} dB against 2.6e7 K and {bracket["low"]:+.1f} dB against 1.3e7 K, the
bracket the requirement set uses.  At 0.6 MHz on the Franklin and MacLean sky
it is {margin_at_0p6:+.1f} dB.  The model is an optimistic bound (no bus mutual
impedance, no plasma sheath), so the receive-efficiency cancellation should not
be assumed for a 5 m antenna at the band bottom.
""".strip()

    return StudyResult(outputs=[path], summary=summary, metrics=metrics)
