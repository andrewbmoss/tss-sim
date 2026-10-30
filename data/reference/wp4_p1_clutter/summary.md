# wp4_campaign

WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2, C4.  Requirements: SCI-05, SCI-07, INS-05, INS-10.

# wp4_campaign: clutter on

6 profiles x 5 seeds on leo_500km_single (30 runs).  Tolerance for a trace match: 2 x column range resolution + 15 km.

| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | O complete med (worst) | on X | false | spurious | O bias km | X bias km |
|---|---|---|---|---|---|---|---|---|---|
| mid_day_max | 1.142 (1.811) | 20.8 (21.5) | 13.1 (11.3) | 0.00 (0.00) | 0.00 | 0.00 | 35 | -8.3 | -7.0 |
| mid_day_min | 0.240 (6.301) | 20.1 (20.2) | 6.4 (6.3) | 0.12 (0.00) | 0.79 | 0.03 | 51 | -4.7 | -5.9 |
| mid_night_min | 8.477 (20.856) | 65.7 (65.7) | -14.1 (-16.4) | 0.00 (0.00) | 0.00 | 0.00 | 48 | -12.6 | -15.9 |
| eq_evening_max | 0.320 (0.913) | nan (nan) | -1.6 (-1.6) | 0.00 (0.00) | 0.00 | 0.00 | 27 | nan | nan |
| high_night | 15.730 (36.887) | 108.0 (108.0) | 5.4 (5.3) | 0.00 (0.00) | 0.00 | 0.00 | 47 | nan | -11.7 |
| smoke | 2.133 (5.229) | nan (nan) | 18.2 (18.1) | 0.00 (0.00) | 0.00 | 0.00 | 70 | nan | nan |
| **all** | 1.477 (36.887) | 20.8 (108.0) | 5.9 (-16.4) | 0.00 (0.00) | 0.00 | 0.00 | 40 | -7.1 | -7.0 |

Verdicts over all runs: SCI05 {'fail': 28, 'threshold': 2}; SCI07 {'fail': 8, 'not_evaluated': 22}; INS10 {'fail': 20, 'threshold': 10}

---

Run 2026-10-30T14:12:30+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_campaign` from the same commit and the same configuration reproduces it.
