"""What the magnetic field costs a topside inversion, and what correcting it buys.

Layer:        study
WP:           WP2, WP6
Claims:       C2
Status:       implemented
Requirements: SCI-05, SCI-07

The ordinary wave's group index near reflection is larger than the field-free
one by the cosecant of the angle between the field and the vertical (Budden
1961, eq. 12.12).  A lamination that ignores the field therefore places a
topside profile wrong, by an amount that grows toward the magnetic poles.  This
study puts numbers on that in two ways and holds the numbers the Technical
Volume quotes.

**The inversion alone.**  An exponential topside with a uniform field has an
exact ordinary-mode trace by quadrature of the magnetoionic group index.  That
trace is inverted with the field-free lamination and with the magnetoionic one,
and each is compared with the true depths.  Nothing else is in the comparison:
no noise, no instrument, no trace extraction.

**The whole chain.**  The pipeline is run from a Chapman truth through the
analytic propagation tier with the magnetoionic group index, the baseline
instrument and signal chain, and each of the two inversions, at a set of
latitudes on one meridian.  Only the ordinary echo is generated, so the result
isolates the inversion from the labelling of modes, which a scalar receiver
cannot do and which is INS-08.  Several noise seeds are run and the median is
reported.  A run in which the trace extractor lost the bottom of the trace is
dropped and counted, because that is a failure of the extractor and not of
either inversion.

What this does not show: the unlabelled trace of a scalar receiver with both
modes present, the error of a real field model, or a structured ionosphere.
"""

from __future__ import annotations

import io
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

#: The exponential topside of the inversion-only test.
F_SAT_HZ = 1.0e6
SCALE_HEIGHT_M = 1.5e5
ALT_SAT_M = 8.0e5
F_H_HZ = 0.9e6
ANGLES_DEG = (5.0, 10.0, 30.0, 50.0, 70.0, 90.0)

#: Latitudes of the whole-chain test, on the prime meridian.
LATITUDES_DEG = (0.0, 20.0, 40.0, 60.0, 80.0)


def exponential_trace_m(freq_hz: np.ndarray, theta_rad: float, *, n: int = 8000) -> np.ndarray:
    """Ordinary-mode virtual range of the exponential topside, by quadrature.

    ``f_N(z) = F exp(z / 2H)`` with ``z`` the depth below the spacecraft, so
    ``dz = 2H d(ln f_N)``.  The substitution ``f_N = f sin(phi)`` removes the
    singularity at reflection.
    """
    from tss.l1_propagation.magnetoionic import group_index

    out = np.empty(freq_hz.size, dtype=float)
    for i, f_k in enumerate(freq_hz):
        p0 = np.arcsin(F_SAT_HZ / f_k)
        phi = p0 + (np.pi / 2 - p0) * (np.arange(n) + 0.5) / n
        mu_g = group_index(np.sin(phi) ** 2, F_H_HZ / f_k, theta_rad, "O")
        out[i] = np.sum(mu_g * 2.0 * SCALE_HEIGHT_M / np.tan(phi)) * (np.pi / 2 - p0) / n
    return out


def inversion_only() -> np.ndarray:
    """Rows of ``(theta_deg, field_free_max_error_m, magnetoionic_max_error_m)``."""
    from tss.l4_retrieval.lamination import (
        constant_field,
        lamina_thicknesses_m,
        lamina_thicknesses_magnetoionic_m,
    )

    freq = F_SAT_HZ * 1.05 ** np.arange(1, 80)
    freq = freq[freq <= F_SAT_HZ * np.exp((ALT_SAT_M - 3.0e5) / (2.0 * SCALE_HEIGHT_M))]
    depth = 2.0 * SCALE_HEIGHT_M * np.log(freq / F_SAT_HZ)
    rows = []
    for theta_deg in ANGLES_DEG:
        theta = np.radians(theta_deg)
        trace = exponential_trace_m(freq, theta)
        plain = np.cumsum(lamina_thicknesses_m(trace, freq, F_SAT_HZ))
        corrected = np.cumsum(
            lamina_thicknesses_magnetoionic_m(
                trace, freq, F_SAT_HZ, ALT_SAT_M, constant_field(F_H_HZ, theta)
            )
        )
        rows.append(
            (theta_deg, float(np.max(np.abs(plain - depth))), float(np.max(np.abs(corrected - depth))))
        )
    return np.array(rows)


def whole_chain(config: dict[str, Any], seeds: tuple[int, ...]) -> np.ndarray:
    """Rows per latitude and inversion; see the module docstring."""
    from tss.interfaces import PlatformState
    from tss.pipeline import Pipeline, build_tiers
    from tss.studies import _grids, _sweep_from_config

    sweep, _ = _sweep_from_config(config, io.StringIO())
    alt_m, freq_hz = _grids(config)
    epoch = str(config.get("epoch", "2026-01-01T12:00:00Z"))

    rows = []
    for lat in LATITUDES_DEG:
        platform = PlatformState(
            alt_m=ALT_SAT_M,
            lat_deg=lat,
            lon_deg=0.0,
            velocity_ecef_mps=np.array([7600.0, 0.0, 0.0]),
            epoch=epoch,
        )
        for j, l4 in enumerate(("lamination", "lamination_magnetoionic")):
            kept: list[tuple[float, float, float]] = []
            theta = float("nan")
            lost = 0
            for seed in seeds:
                tiers = {
                    "l0": "chapman",
                    "l1": {"name": "analytic", "modes": ["O"]},
                    "l2": "baseline",
                    "l3": "baseline",
                    "l4": l4,
                    "l5": "standard",
                }
                pipeline = Pipeline(build_tiers(tiers), seed=int(seed), keep_intermediates=True)
                result = pipeline.run_configuration(
                    "field_inversion", platform=platform, sweep=sweep,
                    freq_hz=freq_hz, alt_m=alt_m, epoch=epoch,
                )
                truth_start = float(
                    pipeline.tiers["l1"].plasma_frequency_at(result.truth, ALT_SAT_M)
                )
                flags = result.retrieved.flags
                field = flags.get("field_model")
                if field:
                    theta = float(field["field_angle_from_vertical_deg"])
                if float(flags.get("start_plasma_hz", np.inf)) > 1.2 * truth_start:
                    lost += 1
                    continue
                errors = result.assessment.errors
                kept.append(
                    (
                        float(errors["ne_rms_fractional"]),
                        float(errors["height_error_rms_m"]),
                        float(errors["height_error_bias_m"]),
                    )
                )
            med = np.median(np.array(kept), axis=0) if kept else np.full(3, np.nan)
            rows.append((lat, j, theta, med[0], med[1], med[2], len(kept), lost))
    out = np.array(rows)
    # The field-free tier does not evaluate a field, so give its rows the angle
    # the magnetoionic tier found at the same latitude.
    for lat in LATITUDES_DEG:
        sel = out[:, 0] == lat
        out[sel, 2] = np.nanmax(out[sel, 2])
    return out


@register_study(
    "field_inversion",
    description="Field-free against magnetoionic lamination, alone and through the whole chain",
    work_packages=["WP2", "WP6"],
    claims=["C2"],
    status="implemented",
    requirements=["SCI-05", "SCI-07"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    opts = dict(config.get("field_inversion", {}))
    seeds = tuple(int(s) for s in opts.get("seeds", (1, 2, 3, 4, 5, 6)))

    alone = inversion_only()
    alone_path = out_dir / "field_inversion_alone.csv"
    np.savetxt(
        alone_path, alone, delimiter=",", comments="",
        header="theta_deg,field_free_max_error_m,magnetoionic_max_error_m",
    )

    chain = whole_chain(config, seeds)
    chain_path = out_dir / "field_inversion_chain.csv"
    np.savetxt(
        chain_path, chain, delimiter=",", comments="",
        header=("lat_deg,magnetoionic,theta_deg,ne_rms_fractional,height_error_rms_m,"
                "height_error_bias_m,n_runs,n_lost_low_band"),
    )

    plain, corrected = chain[chain[:, 1] == 0], chain[chain[:, 1] == 1]
    metrics = {
        "alone_field_free_max_error_km": [float(x) / 1e3 for x in alone[:, 1]],
        "alone_magnetoionic_max_error_km": [float(x) / 1e3 for x in alone[:, 2]],
        "alone_angles_deg": [float(x) for x in alone[:, 0]],
        "chain_latitudes_deg": [float(x) for x in corrected[:, 0]],
        "chain_theta_deg": [float(x) for x in corrected[:, 2]],
        "chain_field_free_height_rms_km": [float(x) / 1e3 for x in plain[:, 4]],
        "chain_magnetoionic_height_rms_km": [float(x) / 1e3 for x in corrected[:, 4]],
        "chain_field_free_ne_rms": [float(x) for x in plain[:, 3]],
        "chain_magnetoionic_ne_rms": [float(x) for x in corrected[:, 3]],
        "seeds": list(seeds),
    }
    summary = f"""
Inverting an exact ordinary-mode trace of an exponential topside, the field-free
lamination is in error by up to {alone[1, 1] / 1e3:.0f} km when the field is 10 degrees from
the vertical and {alone[3, 1] / 1e3:.0f} km at 50 degrees.  The magnetoionic lamination is
within {np.max(alone[:, 2]) / 1e3:.1f} km at every angle, which is the discretisation of the laminae.

Through the whole chain with the ordinary echo alone, the field-free inversion
gives a height error of {np.nanmin(plain[:, 4]) / 1e3:.1f} to {np.nanmax(plain[:, 4]) / 1e3:.1f} km rms across the latitudes run,
and the magnetoionic inversion {np.nanmin(corrected[:, 4]) / 1e3:.1f} to {np.nanmax(corrected[:, 4]) / 1e3:.1f} km.  SCI-07 asks for 20 km
at threshold and 10 km at goal.

A scalar receiver does not label the modes.  With both echoes present the
baseline trace extractor follows the nearer extraordinary echo above its cutoff
at the spacecraft, and that error is not in these numbers.
""".strip()
    return StudyResult(outputs=[alone_path, chain_path], summary=summary, metrics=metrics)
