# wp4_campaign

WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2, C4.  Requirements: SCI-05, SCI-07, INS-05, INS-10.

# wp4_campaign: clutter off

6 profiles x 5 seeds on leo_500km_single (30 runs).  Tolerance for a trace match: 2 x column range resolution + 15 km.

| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | S/(N+C) dB med | O complete med (worst) | on X | false | spurious | O bias km | X bias km |
|---|---|---|---|---|---|---|---|---|---|---|
| mid_day_max | 0.109 (0.109) | 18.1 (18.1) | 14.1 (11.2) | 14.4 | 0.14 (0.14) | 0.86 | 0.00 | 0 | -4.5 | -2.5 |
| mid_day_min | 0.182 (0.188) | 16.1 (16.7) | 22.6 (22.5) | 22.6 | 0.18 (0.18) | 0.82 | 0.00 | 0 | -4.9 | -2.3 |
| mid_night_min | 0.186 (0.191) | 33.0 (34.2) | 9.7 (8.1) | 9.6 | 0.44 (0.44) | 0.44 | 0.11 | 0 | -2.7 | -5.4 |
| eq_evening_max | 0.039 (0.044) | 13.3 (14.7) | 23.3 (23.1) | 23.6 | 0.60 (0.60) | 0.40 | 0.00 | 0 | -6.6 | -13.7 |
| high_night | 1.131 (1.132) | 106.6 (106.6) | 17.4 (17.1) | 17.7 | 0.00 (0.00) | 0.33 | 0.00 | 0 | nan | -4.9 |
| smoke | 0.091 (0.091) | 12.2 (12.2) | 21.9 (21.8) | 22.1 | 0.14 (0.14) | 0.86 | 0.00 | 0 | -4.1 | -2.5 |
| **all** | 0.145 (1.132) | 17.4 (106.6) | 20.2 (8.1) | 20.2 | 0.16 (0.00) | 0.62 | 0.00 | 0 | -4.5 | -3.5 |

Verdicts over all runs: SCI05 {'fail': 15, 'goal': 10, 'threshold': 5}; SCI07 {'fail': 10, 'threshold': 20}; INS10 {'fail': 3, 'threshold': 27}

---

Run 2026-10-31T13:22:34+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_campaign` from the same commit and the same configuration reproduces it.
