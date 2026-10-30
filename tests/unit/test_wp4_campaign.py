"""WP4 Phase 0: the evaluation harness (tss.studies.wp4_campaign)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from tss.interfaces import EchoSet
from tss.studies import wp4_campaign as w


def _echoes(rows):
    f, m, r, fl = zip(*rows, strict=True)
    return EchoSet(
        freq_hz=np.asarray(f, float), mode=np.asarray(m, object),
        group_delay_s=2.0 * np.asarray(r, float) / w.C_LIGHT,
        amplitude_db=np.full(len(f), -100.0), absorption_db=np.zeros(len(f)),
        flags=np.asarray(fl, object),
    )


def _ionogram(freqs, ranges, res=1.0e3):
    return SimpleNamespace(freq_hz=np.asarray(freqs, float),
                           virtual_range_m=np.asarray(ranges, float)[:, None],
                           flags={"range_resolution_m_by_column": [res] * len(freqs)})


def test_trace_measures_classifies_o_x_false_missed_and_spurious():
    e = _echoes([
        (1e6, "O", 200e3, ""), (1e6, "X", 150e3, ""),
        (2e6, "O", 300e3, ""), (2e6, "X", 250e3, ""),
        (3e6, "O", 400e3, ""), (3e6, "X", 350e3, ""),
        (4e6, "O", 500e3, ""),
        (4.5e6, "O", 600e3, ""),  # the O cutoff step (G5c): not scored
        (6e6, "O", np.nan, "penetrating"), (6e6, "X", np.nan, "penetrating"),
    ])
    ion = _ionogram([1e6, 2e6, 3e6, 4e6, 4.5e6, 5e6, 6e6],
                   [195e3, 252e3, 600e3, np.nan, np.nan, 100e3, np.nan])
    m = w.trace_measures(ion, e)
    assert m["n_truth_o_columns"] == 4
    assert (m["n_o"], m["n_x"], m["n_false"], m["n_missed"], m["n_spurious"]) == (1, 1, 1, 1, 1)
    assert m["completeness_o"] == 0.25
    assert np.isclose(m["range_error_o_mean_m"], -5e3)
    assert np.isclose(m["range_error_x_mean_m"], 2e3)


def test_trace_measures_skips_the_cutoff_step():
    # 2 MHz is the O mode's last step below penetration (G5c): not scored.
    e = _echoes([(1e6, "O", 200e3, ""), (2e6, "O", 300e3, ""), (3e6, "O", np.nan, "penetrating")])
    m = w.trace_measures(_ionogram([1e6, 2e6, 3e6], [200e3, np.nan, np.nan]), e)
    assert m["n_truth_o_columns"] == 1 and m["completeness_o"] == 1.0


def test_local_time_sets_longitude():
    assert w._lon_for_local_time(12.0) == 0.0
    assert w._lon_for_local_time(0.0) in (180.0, -180.0)
    assert w._lon_for_local_time(21.0) == 135.0


def test_one_run_end_to_end_is_scored_and_reproducible():
    from tss.pipeline import build_tiers
    from tss.studies import _base_config

    cfg = _base_config(("grid.alt_step_m=20000",))
    shared = build_tiers({k: v for k, v in cfg["tiers"].items() if k not in ("l0", "l2")})
    a = w.run_one("smoke", 1, config=cfg, configuration="leo_500km_single", clutter=None, shared=shared)
    b = w.run_one("smoke", 1, config=cfg, configuration="leo_500km_single", clutter=None, shared=shared)
    assert a["clutter"] == "off" and a["clutter_rows"] == 0
    assert np.isfinite(a["sci05"]) and np.isfinite(a["sci07_m"]) and np.isfinite(a["ins10_db"])
    assert a["n_truth_o_columns"] > 0
    for k in ("sci05", "sci07_m", "ins10_db", "per_column"):
        assert a[k] == b[k]
