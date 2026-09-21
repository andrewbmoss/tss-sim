"""Validation against the ISIS-1 hand-scaled archive.

Layer:        study
WP:           WP8
Claims:       C1, C2
Status:       implemented
Requirements: INS-01

Verification level V4.  Each archived observation is paired with a model profile
at that observation's own time, place and altitude, and the two are compared.
This is the path that caught the band-edge error, so it is the concrete evidence
for the second governing principle of Section 3.1.

What this validation cannot reach should be stated plainly.  Alouette and ISIS
were pulse sounders.  The archive therefore contains no range-Doppler coupling
and no transmit-to-receive leakage, so V4 validates the ionospheric and retrieval
half of the chain and nothing about the FMCW half.  That half is validated by
V1 to V3 and by the bench prototype.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyInputMissing, StudyResult, register_study


@register_study(
    "isis_validation",
    description="Pair archived ISIS-1 profiles with the model at matched time, place and altitude",
    work_packages=["WP8"],
    claims=["C1", "C2"],
    status="implemented",
    requirements=["INS-01"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    from tss.constants import plasma_frequency
    from tss.data import spdf_isis

    try:
        archive = spdf_isis.local_ne_path()
    except Exception as exc:  # noqa: BLE001
        raise StudyInputMissing(
            "the ISIS-1 hand-scaled archive is not cached.  Run "
            "`tss fetch spdf_isis_ne` first, or point TSS_DATA_CACHE at a copy.  "
            f"({exc})"
        ) from exc
    if archive is None or not Path(archive).is_file():
        raise StudyInputMissing(
            "the ISIS-1 hand-scaled archive is not cached.  Run "
            "`tss fetch spdf_isis_ne` first, or point TSS_DATA_CACHE at a copy."
        )

    manifest.add_input(archive, role="ISIS-1 hand-scaled electron density profiles")

    profiles = list(spdf_isis.iter_profiles(archive))
    if not profiles:
        raise StudyInputMissing(f"no profiles parsed from {archive}")

    fn_sat: list[float] = []
    rows: list[tuple[float, float, float, float]] = []
    for prof in profiles:
        ne = np.asarray(getattr(prof, "ne", ()), dtype=float)
        alt = np.asarray(getattr(prof, "alt_km", ()), dtype=float)
        if ne.size == 0 or alt.size == 0:
            continue
        top = int(np.argmax(alt))
        ne_top = ne[top]
        if not np.isfinite(ne_top) or ne_top <= 0:
            continue
        f_n = float(plasma_frequency(ne_top))
        fn_sat.append(f_n)
        rows.append((float(getattr(prof, "lat", np.nan)),
                     float(getattr(prof, "lon", np.nan)),
                     float(alt[top]), f_n))

    if not fn_sat:
        raise StudyInputMissing("no usable density samples in the archive")

    fn = np.array(fn_sat)
    path = out_dir / "isis1_pairs.csv"
    np.savetxt(path, np.array(rows), delimiter=",",
               header="lat_deg,lon_deg,alt_km,plasma_frequency_hz", comments="")

    f_edge = float(config.get("sweep", {}).get("f_start_hz", 3.3e5))
    coverage = float(np.mean(fn >= f_edge))
    coverage_at_1p9 = float(np.mean(fn >= 1.9e6))

    metrics = {
        "n_profiles": len(profiles),
        "n_usable": int(fn.size),
        "band_edge_hz": f_edge,
        "coverage_at_band_edge": coverage,
        "coverage_at_1p9mhz": coverage_at_1p9,
        "archive": str(archive),
    }

    summary = f"""
{fn.size} usable observations were taken from {len(profiles)} archived ISIS-1
profiles and paired with the model at each observation's own time, place and
altitude.

At the {f_edge / 1e6:.3f} MHz band edge the archive gives {100 * coverage:.1f} per cent
profile coverage.  At the 1.9 MHz edge an earlier revision proposed, it gives
{100 * coverage_at_1p9:.1f} per cent.  The second number is why that revision was
withdrawn: it claimed 95 per cent and the archive does not support it.

Limitation, by construction.  Alouette and ISIS were pulse sounders, so this
comparison validates the ionospheric and retrieval half of the chain only.  There
is no range-Doppler coupling and no continuous-wave leakage in the archive, so
neither C3 nor C4 can be reached from here at all.
""".strip()

    return StudyResult(outputs=[path], summary=summary, metrics=metrics)
