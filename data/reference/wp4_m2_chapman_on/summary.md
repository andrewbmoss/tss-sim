# wp4_campaign

WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2, C4.  Requirements: SCI-05, SCI-07, INS-05, INS-10.

# wp4_campaign: clutter on

6 profiles x 5 seeds on leo_500km_single (30 runs).  Tolerance for a trace match: 2 x column range resolution + 15 km.

| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | S/(N+C) dB med | O complete med (worst) | on X | false | spurious | O bias km | X bias km | label err | X samples | gap above hmF2 km med (worst) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mid_day_max | 0.017 (0.024) | 3.0 (3.8) | 14.2 (13.4) | 13.1 | 0.93 (0.86) | 0.00 | 0.03 | 0 | -6.4 | 14.9 | 0.011 | 30 | 25.4 (26.0) |
| mid_day_min | 0.043 (0.048) | 4.1 (4.7) | 22.7 (20.7) | 6.4 | 0.61 (0.58) | 0.18 | 0.00 | 0 | -5.2 | 5.7 | 0.099 | 31 | 9.2 (24.5) |
| mid_night_min | 0.015 (0.015) | 3.0 (4.7) | -10.8 (-12.7) | -13.7 | 0.56 (0.56) | 0.00 | 0.11 | 0 | -7.7 | nan | 0.000 | 7 | 54.6 (54.7) |
| eq_evening_max | 0.002 (0.032) | 0.7 (8.3) | 23.4 (21.5) | -1.6 | 0.80 (0.40) | 0.00 | 0.00 | 0 | -6.6 | nan | 0.000 | 8 | 14.3 (41.5) |
| high_night | 0.011 (0.013) | 1.4 (1.5) | 23.3 (20.0) | 5.1 | 0.83 (0.58) | 0.00 | 0.00 | 0 | -4.8 | nan | 0.000 | 9 | 32.7 (110.7) |
| smoke | 0.029 (0.033) | 3.5 (3.8) | 24.2 (21.7) | 18.0 | 0.59 (0.59) | 0.14 | 0.00 | 0 | -4.3 | 4.0 | 0.121 | 22 | 27.6 (28.6) |
| **all** | 0.016 (0.048) | 3.0 (8.3) | 21.6 (-12.7) | 5.8 | 0.67 (0.40) | 0.00 | 0.00 | 0 | -5.4 | 5.3 | 0.057 | 18 | 25.8 (110.7) |

Verdicts over all runs: SCI05 {'goal': 30}; SCI07 {'goal': 30}; INS10 {'fail': 5, 'threshold': 25}

---

Run 2026-11-06T07:47:34+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_campaign` from the same commit and the same configuration reproduces it.
