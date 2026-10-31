# wp4_campaign

WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2, C4.  Requirements: SCI-05, SCI-07, INS-05, INS-10.

# wp4_campaign: clutter off

6 profiles x 5 seeds on leo_500km_single (30 runs).  Tolerance for a trace match: 2 x column range resolution + 15 km.

| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | S/(N+C) dB med | O complete med (worst) | on X | false | spurious | O bias km | X bias km |
|---|---|---|---|---|---|---|---|---|---|---|
| mid_day_max | 0.125 (0.135) | 20.2 (21.6) | 14.0 (11.8) | 13.7 | 0.10 (0.10) | 0.86 | 0.03 | 0 | -11.6 | -7.1 |
| mid_day_min | 0.219 (0.222) | 18.7 (18.8) | 11.8 (11.5) | 11.8 | 0.15 (0.15) | 0.82 | 0.03 | 0 | -5.2 | -6.0 |
| mid_night_min | 0.548 (0.549) | 93.9 (94.0) | -13.6 (-14.1) | -13.6 | 0.00 (0.00) | 0.11 | 0.33 | 0 | nan | -16.6 |
| eq_evening_max | 0.047 (0.049) | 16.3 (16.9) | 19.4 (18.4) | 18.9 | 0.20 (0.20) | 0.00 | 0.80 | 0 | -1.7 | nan |
| high_night | 1.171 (1.174) | 108.3 (108.4) | 16.5 (13.7) | 15.9 | 0.00 (0.00) | 0.17 | 0.17 | 0 | nan | -14.2 |
| smoke | 0.116 (0.127) | 14.9 (15.0) | 21.6 (20.2) | 21.2 | 0.14 (0.14) | 0.86 | 0.00 | 0 | -13.4 | -7.5 |
| **all** | 0.177 (1.174) | 18.8 (108.4) | 14.7 (-14.1) | 14.4 | 0.12 (0.00) | 0.48 | 0.12 | 0 | -8.3 | -7.5 |

Verdicts over all runs: SCI05 {'fail': 15, 'goal': 5, 'threshold': 10}; SCI07 {'fail': 13, 'threshold': 17}; INS10 {'fail': 5, 'threshold': 25}

---

Run 2026-10-31T03:40:23+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_campaign` from the same commit and the same configuration reproduces it.
