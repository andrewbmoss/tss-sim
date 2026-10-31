# wp4_campaign

WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2, C4.  Requirements: SCI-05, SCI-07, INS-05, INS-10.

# wp4_campaign: clutter on

6 profiles x 5 seeds on leo_500km_single (30 runs).  Tolerance for a trace match: 2 x column range resolution + 15 km.

| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | S/(N+C) dB med | O complete med (worst) | on X | false | spurious | O bias km | X bias km |
|---|---|---|---|---|---|---|---|---|---|---|
| mid_day_max | 0.126 (0.135) | 20.2 (21.5) | 15.1 (12.0) | 13.1 | 0.07 (0.07) | 0.86 | 0.07 | 3 | -9.0 | -7.1 |
| mid_day_min | 0.323 (5.452) | 25.4 (25.4) | 21.8 (21.2) | 6.4 | 0.12 (0.00) | 0.79 | 0.09 | 6 | -5.7 | -6.0 |
| mid_night_min | 12.267 (32.323) | 108.3 (108.3) | -10.8 (-13.7) | -14.1 | 0.00 (0.00) | 0.00 | 0.00 | 31 | nan | -16.6 |
| eq_evening_max | 0.856 (0.913) | 17.4 (17.4) | 20.3 (17.3) | -1.6 | 0.00 (0.00) | 0.00 | 0.00 | 14 | -17.2 | nan |
| high_night | 22.383 (36.601) | nan (nan) | 16.6 (13.8) | 5.4 | 0.00 (0.00) | 0.00 | 0.00 | 32 | nan | nan |
| smoke | 4.617 (12.036) | 14.9 (14.9) | 21.8 (21.6) | 18.2 | 0.00 (0.00) | 0.00 | 0.00 | 27 | -13.6 | -7.4 |
| **all** | 0.874 (36.601) | 20.2 (108.3) | 18.1 (-13.7) | 5.9 | 0.00 (0.00) | 0.00 | 0.00 | 22 | -9.0 | -7.1 |

Verdicts over all runs: SCI05 {'fail': 22, 'goal': 1, 'threshold': 7}; SCI07 {'fail': 8, 'not_evaluated': 17, 'threshold': 5}; INS10 {'fail': 5, 'threshold': 25}

---

Run 2026-10-31T03:39:52+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_campaign` from the same commit and the same configuration reproduces it.
