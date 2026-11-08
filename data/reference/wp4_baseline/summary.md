# wp4_campaign

WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2, C4.  Requirements: SCI-05, SCI-07, INS-05, INS-10.

# wp4_campaign: clutter off

6 profiles x 5 seeds on leo_500km_single (30 runs).  Tolerance for a trace match: 2 x column range resolution + 15 km.

| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | S/(N+C) dB med | O complete med (worst) | on X | false | spurious | O bias km | X bias km | label err | X samples | gap above hmF2 km med (worst) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mid_day_max | 0.124 (0.125) | 20.1 (20.2) | 13.7 (13.4) | 13.7 | 0.10 (0.10) | 0.86 | 0.03 | 0 | -11.4 | -7.0 | 0.893 | 0 | 25.7 (48.2) |
| mid_day_min | 0.219 (0.220) | 18.7 (18.8) | 11.7 (11.0) | 11.7 | 0.15 (0.15) | 0.82 | 0.03 | 0 | -5.1 | -5.9 | 0.844 | 0 | 19.7 (27.4) |
| mid_night_min | 0.547 (0.548) | 93.8 (93.9) | -16.0 (-16.6) | -16.0 | 0.00 (0.00) | 0.11 | 0.33 | 0 | nan | -16.1 | 1.000 | 0 | 61.3 (62.0) |
| eq_evening_max | 0.046 (0.046) | 16.2 (16.3) | 18.5 (18.4) | 18.5 | 0.20 (0.20) | 0.00 | 0.80 | 0 | -1.7 | nan | 0.000 | 0 | 14.2 (14.5) |
| high_night | 1.170 (1.228) | 108.2 (111.7) | 17.2 (15.5) | 17.2 | 0.00 (0.00) | 0.17 | 0.17 | 0 | nan | -14.5 | 1.000 | 0 | 62.1 (64.4) |
| smoke | 0.116 (0.117) | 14.9 (15.0) | 20.9 (20.5) | 20.9 | 0.14 (0.09) | 0.86 | 0.00 | 0 | -13.4 | -7.4 | 0.872 | 0 | 15.9 (16.1) |
| **all** | 0.171 (1.228) | 19.4 (111.7) | 15.6 (-16.6) | 15.6 | 0.10 (0.00) | 0.49 | 0.11 | 0 | -8.4 | -7.4 | 0.862 | 0 | 25.6 (64.4) |

Verdicts over all runs: SCI05 {'fail': 15, 'goal': 5, 'threshold': 10}; SCI07 {'fail': 15, 'threshold': 15}; INS10 {'fail': 5, 'threshold': 25}

---

Run 2026-11-08T16:22:41+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_campaign` from the same commit and the same configuration reproduces it.
