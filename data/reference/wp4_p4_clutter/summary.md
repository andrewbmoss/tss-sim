# wp4_campaign

WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2, C4.  Requirements: SCI-05, SCI-07, INS-05, INS-10.

# wp4_campaign: clutter on

6 profiles x 5 seeds on leo_500km_single (30 runs).  Tolerance for a trace match: 2 x column range resolution + 15 km.

| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | S/(N+C) dB med | O complete med (worst) | on X | false | spurious | O bias km | X bias km | label err | X samples | gap above hmF2 km med (worst) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| mid_day_max | 0.034 (0.039) | 5.1 (5.7) | 14.2 (13.4) | 13.1 | 0.93 (0.86) | 0.00 | 0.03 | 0 | -6.4 | 14.9 | 0.011 | 30 | 26.6 (27.0) |
| mid_day_min | 0.107 (0.114) | 8.8 (9.5) | 22.7 (20.7) | 6.4 | 0.61 (0.58) | 0.18 | 0.00 | 0 | -5.2 | 5.7 | 0.099 | 31 | 11.1 (26.6) |
| mid_night_min | 0.063 (0.064) | 11.1 (11.9) | -10.8 (-12.7) | -13.7 | 0.56 (0.56) | 0.00 | 0.11 | 0 | -7.7 | nan | 0.000 | 7 | 59.6 (59.7) |
| eq_evening_max | 0.017 (0.039) | 5.8 (12.1) | 23.4 (21.5) | -1.6 | 0.80 (0.40) | 0.00 | 0.00 | 0 | -6.6 | nan | 0.000 | 8 | 16.5 (44.3) |
| high_night | 0.067 (0.076) | 7.9 (8.6) | 23.3 (20.0) | 5.1 | 0.83 (0.58) | 0.00 | 0.00 | 0 | -4.8 | nan | 0.000 | 9 | 36.6 (116.3) |
| smoke | 0.064 (0.069) | 6.9 (7.2) | 24.2 (21.7) | 18.0 | 0.59 (0.59) | 0.14 | 0.00 | 0 | -4.3 | 4.0 | 0.121 | 22 | 29.0 (30.1) |
| **all** | 0.063 (0.114) | 7.9 (12.1) | 21.6 (-12.7) | 5.8 | 0.67 (0.40) | 0.00 | 0.00 | 0 | -5.4 | 5.3 | 0.057 | 18 | 26.9 (116.3) |

Verdicts over all runs: SCI05 {'goal': 25, 'threshold': 5}; SCI07 {'goal': 24, 'threshold': 6}; INS10 {'fail': 5, 'threshold': 25}

---

Run 2026-10-31T20:52:48+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_campaign` from the same commit and the same configuration reproduces it.
