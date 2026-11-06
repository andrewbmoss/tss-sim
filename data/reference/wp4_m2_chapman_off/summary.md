# wp4_campaign

WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2, C4.  Requirements: SCI-05, SCI-07, INS-05, INS-10.

# wp4_campaign: clutter off

6 profiles x 5 seeds on leo_500km_single (30 runs).  Tolerance for a trace match: 2 x column range resolution + 15 km.

| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | S/(N+C) dB med | O complete med (worst) | on X | false | spurious | O bias km | X bias km | label err | X samples | gap above hmF2 km med (worst) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mid_day_max | 0.005 (0.005) | 0.9 (0.9) | 14.1 (11.2) | 14.3 | 1.00 (1.00) | 0.00 | 0.00 | 0 | -3.3 | nan | 0.000 | 29 | 23.2 (23.3) |
| mid_day_min | 0.023 (0.025) | 2.1 (2.2) | 22.8 (22.6) | 22.8 | 1.00 (0.94) | 0.00 | 0.00 | 0 | -3.2 | nan | 0.000 | 31 | 21.6 (37.3) |
| mid_night_min | 0.017 (0.018) | 2.9 (5.9) | 11.3 (8.1) | 10.7 | 0.78 (0.78) | 0.00 | 0.11 | 0 | -7.2 | nan | 0.000 | 7 | 52.8 (53.3) |
| eq_evening_max | 0.029 (0.029) | 6.5 (6.6) | 23.3 (23.1) | 23.4 | 0.60 (0.60) | 0.00 | 0.00 | 0 | -6.5 | nan | 0.000 | 0 | 3.0 (3.1) |
| high_night | 0.033 (0.033) | 4.1 (4.1) | 21.2 (20.2) | 20.9 | 1.00 (1.00) | 0.00 | 0.00 | 0 | -4.7 | nan | 0.000 | 10 | 35.1 (35.2) |
| smoke | 0.011 (0.012) | 1.4 (1.5) | 25.5 (25.0) | 25.6 | 0.91 (0.82) | 0.00 | 0.00 | 0 | -3.2 | nan | 0.020 | 21 | 26.5 (46.1) |
| **all** | 0.021 (0.033) | 2.5 (6.6) | 22.2 (8.1) | 22.2 | 0.95 (0.60) | 0.00 | 0.00 | 0 | -3.9 | nan | 0.004 | 16 | 26.3 (53.3) |

Verdicts over all runs: SCI05 {'goal': 30}; SCI07 {'goal': 30}; INS10 {'fail': 2, 'threshold': 28}

---

Run 2026-11-06T07:42:58+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_campaign` from the same commit and the same configuration reproduces it.
