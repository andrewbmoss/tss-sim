# wp5_campaign

WP5 evaluation harness: C3 scenes by echo class and the route-one budget (Phase 0 baseline)

Status: **implemented**.  Work packages: WP5.  Claims: C3.  Requirements: SCI-15, INS-15, INS-26.

# wp5_campaign: WP5 Phase 0 baseline (parameter level), seed 20261124

117 scenes; budget 12 rows at 1000 trials.
Reproduced: **True**.

- A_pass: 56/56
- B_flagged: 0/32
- C_declined: 12/12
- E_pass: 2/2
- A_max_abs_drift_err_mps: 0.03003
- A_degeneracy_max_ratio: 2.223e-17
- B_err_over_predicted_min: 0.9994
- B_err_over_predicted_max: 1.001
- C_max_abs_slope_ratio: 0.0006038
- sheet_given_max_abs_err_free_mps: 0.1467
- sheet_given_max_abs_err_plasma_mps: 0.03393
- sheet_assume_err_psi15_phi5_mps: -166.5
- sheet_assume_err_psi45_phi5_mps: -464.9
- sheet_assume_abs_err_meridian_dip10_phi5_mps: 5.021
- sheet_assume_err_plasma_psi15_mps: -87.34
- sheet_assume_abs_err_plasma_meridian_mps: 2.557
- E_max_abs_drift_err_mps: 0.0157
- E_max_abs_group_path_dev_km: 0.0003792
- budget_max_abs_deviation: 0.004261
- D (slope ratio, class): {'0.5km': [1.0, 'compact_facet'], '1km': [0.999, 'compact_facet'], '2km': [0.985, 'compact_facet'], '4km': [0.766, 'extended'], '8km': [0.084, 'extended'], '16km': [-0.006, 'extended'], '32km': [-0.006, 'extended']}

| f MHz | P' km | fading | total m/s | TV m/s | deviation |
|---|---|---|---|---|---|
| 3 | 500 | none | 12.19 | 12.2 | -0.1 % |
| 3 | 500 | tau 0.5 s | 33.67 | 33.7 | -0.1 % |
| 3 | 1000 | none | 11.31 | 11.3 | +0.1 % |
| 3 | 1000 | tau 0.5 s | 30.86 | 30.9 | -0.1 % |
| 5 | 500 | none | 10.81 | 10.8 | +0.1 % |
| 5 | 500 | tau 0.5 s | 23.15 | 23.1 | +0.2 % |
| 5 | 1000 | none | 10.40 | 10.4 | -0.0 % |
| 5 | 1000 | tau 0.5 s | 20.11 | 20.1 | +0.1 % |
| 10 | 500 | none | 10.63 | 10.6 | +0.3 % |
| 10 | 500 | tau 0.5 s | 14.83 | 14.8 | +0.2 % |
| 10 | 1000 | none | 10.46 | 10.5 | -0.4 % |
| 10 | 1000 | tau 0.5 s | 14.47 | 14.5 | -0.2 % |

---

Run 2026-11-25T02:25:03+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp5_campaign` from the same commit and the same configuration reproduces it.
