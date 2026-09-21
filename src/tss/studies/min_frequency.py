"""The lower edge of the sounding band.

Layer:        study
WP:           WP0, WP1
Claims:       C1, C2
Status:       implemented
Requirements: INS-01

A sounder cannot radiate below the local plasma frequency at its own position,
so the lower edge of the band is not a design choice.  It is a statistical
property of the ionosphere at the orbit altitude, and it is the parameter that
drives antenna size and therefore the whole small-platform argument.  Alouette I
began at 0.45 MHz, the Radio Plasma Imager at 3 kHz and the TOPside Automated
Sounder at 0.5 MHz, all round numbers with no supporting distribution behind
them.  This study derives the number instead.

THE COMBINATION DIRECTION IS THE WHOLE RESULT.  The per-model fifth percentiles
must be combined by taking the MINIMUM, not the maximum.  Taking the maximum
picks the most optimistic model at each point and produces 1.9 MHz, which
delivered 53.4 per cent profile coverage against the 95 per cent claimed.  The
minimum is the conservative envelope across a disagreeing ensemble and produces
0.33 MHz.  An earlier revision of this work took the maximum, passed every
internal check, and was caught only by the ISIS-1 archive comparison.  That is
the evidence for the second governing principle of Section 3.1, and
:func:`combine_ensemble` carries a test that pins the direction so the error
cannot come back.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyInputMissing, StudyResult, register_study

#: Coverage fractions the requirement is stated at.  Threshold then goal.
COVERAGE_LEVELS = (0.95, 0.99)


def combine_ensemble(per_model_percentiles: np.ndarray) -> float:
    """Combine per-model percentiles into one band edge.

    Takes the minimum across models, which is the conservative envelope.  Do not
    change this to a maximum.  See the module docstring: the maximum gave
    1.9 MHz and 53.4 per cent real coverage against 95 per cent claimed.
    """
    values = np.asarray(per_model_percentiles, dtype=float)
    if values.size == 0:
        raise ValueError("cannot combine an empty ensemble")
    return float(np.nanmin(values))


def coverage_fraction(fn_sat_hz: np.ndarray, f_min_hz: float) -> float:
    """Fraction of observations a sounder starting at ``f_min_hz`` can reach."""
    fn = np.asarray(fn_sat_hz, dtype=float)
    fn = fn[np.isfinite(fn)]
    if fn.size == 0:
        return float("nan")
    return float(np.mean(fn >= f_min_hz))


def f_min_for_coverage(fn_sat_hz: np.ndarray, coverage: float) -> float:
    """The lowest start frequency that reaches ``coverage`` of observations."""
    fn = np.asarray(fn_sat_hz, dtype=float)
    fn = fn[np.isfinite(fn)]
    if fn.size == 0:
        return float("nan")
    return float(np.percentile(fn, 100.0 * (1.0 - coverage)))


@register_study(
    "min_frequency",
    description="Derive the lower band edge from the modelled plasma frequency distribution",
    work_packages=["WP0", "WP1"],
    claims=["C1", "C2"],
    status="implemented",
    requirements=["INS-01"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    from tss.constants import plasma_frequency
    from tss.registry import UnknownTier, build

    altitudes_m = np.array([4.0e5, 5.0e5, 6.0e5, 8.0e5])
    lats = np.arange(-80.0, 81.0, 20.0)
    lons = np.arange(-180.0, 180.0, 45.0)
    epochs = [f"2014-{m:02d}-15T{h:02d}:00:00Z" for m in (1, 4, 7, 10) for h in (0, 6, 12, 18)]

    tier_names = ["chapman"]
    for candidate in ("iri2020", "analytic"):
        try:
            build("l0", candidate)
        except (UnknownTier, Exception):  # noqa: BLE001 - availability probe
            continue
        tier_names.append(candidate)

    per_model: dict[str, list[float]] = {}
    rows: list[tuple[float, float, float, float]] = []

    for tier_name in tier_names:
        try:
            tier = build("l0", tier_name)
        except Exception as exc:  # noqa: BLE001
            raise StudyInputMissing(
                f"L0 tier {tier_name!r} could not be built: {exc}"
            ) from exc
        for h_sat in altitudes_m:
            fn_samples: list[float] = []
            # state_at wants a grid, not a point.  A thin three-sample shell
            # centred on the satellite is enough: only the middle sample is
            # read, and a wider grid would cost time for nothing.
            grid = np.array([h_sat - 5.0e3, h_sat, h_sat + 5.0e3])
            for lat in lats:
                for lon in lons:
                    for epoch in epochs:
                        try:
                            state = tier.state_at(float(lat), float(lon), epoch, grid)
                        except NotImplementedError:
                            raise StudyInputMissing(
                                f"L0 tier {tier_name!r} is registered but not callable; "
                                "install the model it needs or select another tier."
                            ) from None
                        ne = float(np.asarray(state.ne).ravel()[1])
                        if np.isfinite(ne) and ne > 0:
                            fn_samples.append(float(plasma_frequency(ne)))
            if not fn_samples:
                continue
            fn = np.array(fn_samples)
            for coverage in COVERAGE_LEVELS:
                rows.append((h_sat, coverage, f_min_for_coverage(fn, coverage),
                             float(fn.size)))
            per_model.setdefault(f"{tier_name}", []).append(
                f_min_for_coverage(fn, COVERAGE_LEVELS[0])
            )

    if not rows:
        raise StudyInputMissing("no usable profiles were produced by any L0 tier")

    path = out_dir / "minfreq_summary.csv"
    np.savetxt(
        path, np.array(rows), delimiter=",",
        header="h_sat_m,coverage,f_min_hz,n_samples", comments="",
    )

    ensemble = np.array([v for vals in per_model.values() for v in vals])
    f_min_threshold = combine_ensemble(ensemble)
    unsafe = float(np.nanmax(ensemble))

    metrics = {
        "models": tier_names,
        "n_samples": int(sum(r[3] for r in rows) / len(COVERAGE_LEVELS)),
        "f_min_hz_at_95pc": f_min_threshold,
        "f_min_mhz_at_95pc": f_min_threshold / 1e6,
        "combination": "minimum of per-model 5th percentiles",
        "unsafe_maximum_combination_hz": unsafe,
        "requirement": "INS-01",
    }

    summary = f"""
The plasma frequency at the spacecraft was sampled over {len(tier_names)} model(s),
{altitudes_m.size} altitudes, {lats.size * lons.size} locations and {len(epochs)} epochs.

Combining the per-model fifth percentiles by taking the minimum gives a lower
band edge of {f_min_threshold / 1e6:.3f} MHz at 95 per cent profile coverage.

For comparison, combining them by taking the maximum would give
{unsafe / 1e6:.3f} MHz.  That is the error an earlier revision made.  It passed
every internal check and was caught only by comparison against hand-scaled ISIS-1
profiles, which showed the optimistic figure delivering 53.4 per cent real
coverage rather than the 95 per cent it claimed.  combine_ensemble() takes the
minimum and a unit test pins that direction.

Caveat.  This is a modelled distribution.  The archive comparison in the
isis_validation study is what makes it trustworthy, and the two should be read
together.
""".strip()

    return StudyResult(outputs=[path], summary=summary, metrics=metrics)
