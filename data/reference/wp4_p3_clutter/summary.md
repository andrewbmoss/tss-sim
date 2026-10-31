# wp4_campaign

WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2, C4.  Requirements: SCI-05, SCI-07, INS-05, INS-10.

# wp4_campaign: clutter on

6 profiles x 5 seeds on leo_500km_single (30 runs).  Tolerance for a trace match: 2 x column range resolution + 15 km.

| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | S/(N+C) dB med | O complete med (worst) | on X | false | spurious | O bias km | X bias km |
|---|---|---|---|---|---|---|---|---|---|---|
| mid_day_max | 0.112 (1.180) | 18.4 (18.4) | 14.2 (13.4) | 13.1 | 0.10 (0.00) | 0.86 | 0.03 | 37 | -5.6 | -3.0 |
| mid_day_min | 0.819 (9.188) | 51.7 (51.7) | 22.7 (20.7) | 6.4 | 0.00 (0.00) | 0.67 | 0.00 | 21 | nan | -2.1 |
| mid_night_min | 12.267 (12.417) | 106.6 (106.6) | -10.8 (-12.7) | -13.9 | 0.00 (0.00) | 0.00 | 0.00 | 48 | nan | -6.6 |
| eq_evening_max | 0.820 (0.913) | nan (nan) | 23.4 (21.5) | -1.6 | 0.00 (0.00) | 0.00 | 0.00 | 26 | nan | nan |
| high_night | 9.015 (22.605) | nan (nan) | 17.8 (17.2) | 5.1 | 0.00 (0.00) | 0.00 | 0.00 | 51 | nan | nan |
| smoke | 6.728 (12.040) | nan (nan) | 21.7 (21.1) | 18.0 | 0.00 (0.00) | 0.00 | 0.00 | 30 | nan | nan |
| **all** | 2.337 (22.605) | 51.6 (106.6) | 20.1 (-12.7) | 5.8 | 0.00 (0.00) | 0.00 | 0.00 | 37 | -5.6 | -3.0 |

Verdicts over all runs: SCI05 {'fail': 26, 'threshold': 4}; SCI07 {'fail': 5, 'not_evaluated': 21, 'threshold': 4}; INS10 {'fail': 5, 'threshold': 25}

---

Run 2026-10-31T13:19:54+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_campaign` from the same commit and the same configuration reproduces it.
