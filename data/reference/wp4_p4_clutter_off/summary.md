# wp4_campaign

WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2, C4.  Requirements: SCI-05, SCI-07, INS-05, INS-10.

# wp4_campaign: clutter off

6 profiles x 5 seeds on leo_500km_single (30 runs).  Tolerance for a trace match: 2 x column range resolution + 15 km.

| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | S/(N+C) dB med | O complete med (worst) | on X | false | spurious | O bias km | X bias km | label err | X samples | gap above hmF2 km med (worst) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mid_day_max | 0.010 (0.010) | 1.5 (1.5) | 14.1 (11.2) | 14.3 | 1.00 (1.00) | 0.00 | 0.00 | 0 | -3.3 | nan | 0.000 | 29 | 23.6 (23.7) |
| mid_day_min | 0.025 (0.026) | 2.2 (2.3) | 22.8 (22.6) | 22.8 | 1.00 (0.94) | 0.00 | 0.00 | 0 | -3.2 | nan | 0.000 | 31 | 21.7 (37.4) |
| mid_night_min | 0.060 (0.061) | 10.5 (12.2) | 11.3 (8.1) | 10.7 | 0.78 (0.78) | 0.00 | 0.11 | 0 | -7.2 | nan | 0.000 | 7 | 59.0 (59.2) |
| eq_evening_max | 0.036 (0.036) | 10.3 (10.4) | 23.3 (23.1) | 23.4 | 0.60 (0.60) | 0.00 | 0.00 | 0 | -6.5 | nan | 0.000 | 0 | 5.6 (5.6) |
| high_night | 0.065 (0.065) | 7.8 (7.8) | 21.2 (20.2) | 20.9 | 1.00 (1.00) | 0.00 | 0.00 | 0 | -4.7 | nan | 0.000 | 7 | 37.2 (37.3) |
| smoke | 0.021 (0.021) | 2.3 (2.5) | 25.5 (25.0) | 25.6 | 0.91 (0.82) | 0.00 | 0.00 | 0 | -3.2 | nan | 0.020 | 21 | 26.9 (46.6) |
| **all** | 0.031 (0.065) | 5.1 (12.2) | 22.2 (8.1) | 22.2 | 0.95 (0.60) | 0.00 | 0.00 | 0 | -3.9 | nan | 0.004 | 14 | 26.8 (59.2) |

Verdicts over all runs: SCI05 {'goal': 30}; SCI07 {'goal': 20, 'threshold': 10}; INS10 {'fail': 2, 'threshold': 28}

---

Run 2026-10-31T20:50:13+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_campaign` from the same commit and the same configuration reproduces it.
