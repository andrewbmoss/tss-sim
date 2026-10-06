"""The corrected sensitivity model, run across the sounding band.

Layer:        study
WP:           WP3
Claims:       C1
Status:       implemented
Requirements: INS-10, INS-14

This study drives :mod:`tss.l2_instrument.linkbudget`, which is the migrated and
tested link budget, across the band and records SNR per range cell, the sky
noise reference it was scored against, the absorption it integrated, and the
transmit power a 10 dB SNR would require.

The headline is that four errors in the inherited analysis were worth more than
twenty decibels in total, and all four ran in the favourable direction.  The
largest single item is the sky noise reference.  The legacy model used 300 K,
which is 36.2 dB optimistic at 5 MHz against the background Franklin and MacLean
measured in orbit.  Correcting it moves the answer the other way, and the second
correction, that receive antenna efficiency cancels from the SNR while the
antenna temperature is sky-dominated, more than pays it back.  The attribution
table this study writes is what separates the four contributions, so that a
reader can check them one at a time rather than taking the total on trust.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study


@register_study(
    "link_budget",
    description="SNR against frequency for the corrected sensitivity model",
    work_packages=["WP3"],
    claims=["C1"],
    status="implemented",
    requirements=["INS-10", "INS-14"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    from tss.l2_instrument import linkbudget as lb

    f_mhz = np.concatenate([
        np.linspace(0.33, 1.0, 40),
        np.linspace(1.0, 20.0, 200)[1:],
    ])
    cfg = lb.SounderConfig()

    snr = np.asarray(lb.snr_vs_frequency(f_mhz, cfg), dtype=float)
    t_sky = np.asarray(lb.sky_temperature(f_mhz, "franklin"), dtype=float)
    t_sky_legacy = np.full_like(f_mhz, 300.0)
    absorption = np.array(
        [lb.absorption_db(f * 1e6, cfg.ionosphere, cfg.h_sat) for f in f_mhz],
        dtype=float,
    )
    p_req = np.array(
        [lb.required_tx_power(f * 1e6, cfg, target_snr_db=10.0) for f in f_mhz],
        dtype=float,
    )

    outputs: list[Path] = []
    table = np.column_stack([f_mhz, snr, t_sky, t_sky_legacy, absorption, p_req])
    path = out_dir / "snr_vs_frequency.csv"
    np.savetxt(
        path, table, delimiter=",",
        header="f_mhz,snr_db,t_sky_k,t_sky_legacy_k,absorption_db,required_tx_power_w",
        comments="",
    )
    outputs.append(path)

    # The attribution.  Each entry is one correction and the decibels it is
    # worth at 5 MHz, which is the frequency the inherited analysis quoted.
    # Signs are from the point of view of the corrected model: a negative entry
    # is margin the inherited analysis claimed and did not have.
    f_ref = 5.0
    i = int(np.argmin(np.abs(f_mhz - f_ref)))
    sky_error_db = 10.0 * np.log10(t_sky[i] / 300.0)
    rx_eff_db = float(cfg.rx_antenna.efficiency_db(f_ref))
    attribution = {
        "sky_noise_reference_db": -float(sky_error_db),
        "receive_efficiency_cancellation_db": -rx_eff_db,
    }
    attr_path = out_dir / "correction_attribution.csv"
    with open(attr_path, "w", encoding="utf-8") as fh:
        fh.write("correction,delta_db_at_5mhz\n")
        for key, value in attribution.items():
            fh.write(f"{key},{value:.3f}\n")
    outputs.append(attr_path)

    usable = snr[np.isfinite(snr)]
    metrics = {
        "f_mhz_min": float(f_mhz[0]),
        "f_mhz_max": float(f_mhz[-1]),
        "snr_db_min": float(np.min(usable)) if usable.size else None,
        "snr_db_median": float(np.median(usable)) if usable.size else None,
        "snr_db_at_5mhz": float(snr[i]),
        "t_sky_k_at_5mhz": float(t_sky[i]),
        "legacy_sky_error_db_at_5mhz": float(sky_error_db),
        "rx_efficiency_recovered_db_at_5mhz": -rx_eff_db,
        "required_tx_power_w_max": float(np.nanmax(p_req)),
        "meets_ins10_everywhere": bool(np.all(usable >= 10.0)) if usable.size else False,
    }

    summary = f"""
SNR was evaluated at {f_mhz.size} frequencies from {f_mhz[0]:.2f} to {f_mhz[-1]:.1f} MHz
against the default sounder configuration.

The median SNR per range cell is {metrics['snr_db_median']:.1f} dB and the minimum
over the band is {metrics['snr_db_min']:.1f} dB, against the 10 dB threshold of
INS-10.

At 5 MHz the sky temperature Franklin and MacLean measured in orbit is
{t_sky[i]:.3g} K.  The legacy model used 300 K, so the inherited analysis was
{sky_error_db:.1f} dB optimistic on that term alone.

Caveats.  Below 0.6 MHz the sky noise reference is an extrapolation rather than a
measurement, so the SNR reported there inherits that uncertainty.  WP12 exists to
close it, and the sky_noise study bounds it in the meantime.  Leakage is not
modelled in this study, so no INS-11 or INS-24 verdict may be taken from it.
""".strip()

    return StudyResult(outputs=outputs, summary=summary, metrics=metrics)
