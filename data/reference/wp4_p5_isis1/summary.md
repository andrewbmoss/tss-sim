# isis1_retrieval

WP4 Phase 5: every usable ISIS-1 topside profile through the P4 chain, clutter on and off

Status: **implemented**.  Work packages: WP4, WP8.  Claims: C1, C3.  Requirements: SCI-05, SCI-07, INS-10.

7790 ISIS-1 records (iquality [4, 5, 6, 7, 8, 9]), clutter off, on; 38.95 h wall for this session (0 rows resumed from the checkpoint), 662.14 CPU h in all, on 17 workers.

| Clutter | N (failed) | SCI-05 median (p90, worst) | SCI-05 verdicts | SCI-07 km median (p90, worst) | SCI-07 verdicts | Span N_e error median (p90) | Span reached median (p10) | Coverage gap km median (p90) | O completeness median |
|---|---|---|---|---|---|---|---|---|---|
| off | 7790 (0) | 0.078 (0.571, 104.687) | {'fail': 2151, 'goal': 4680, 'not_evaluated': 5, 'threshold': 954} | 13.0 (56.0, 595.2) | {'fail': 2315, 'goal': 2879, 'not_evaluated': 5, 'threshold': 2591} | 0.064 (0.295) | 1.00 (0.90) | 27 (88) | 0.43 |
| on | 7790 (0) | 0.092 (0.510, 1404.083) | {'fail': 2118, 'goal': 4278, 'not_evaluated': 71, 'threshold': 1323} | 15.5 (51.1, 755.3) | {'fail': 2716, 'goal': 1977, 'not_evaluated': 96, 'threshold': 3001} | 0.082 (0.276) | 1.00 (0.81) | 18 (126) | 0.39 |

By stratum:

| Clutter | Group | N | SCI-05 median (p90) | SCI-07 km median (p90) | Span N_e error median | Coverage gap km median | SCI-05 fails |
|---|---|---|---|---|---|---|---|
| off | high/day/high_sun | 593 | 0.082 (0.771) | 17.6 (116.8) | 0.071 | 43 | 210 |
| off | high/day/low_sun | 41 | 0.145 (29.434) | 37.2 (410.9) | 0.156 | 70 | 19 |
| off | high/night/high_sun | 855 | 0.112 (1.204) | 20.9 (128.3) | 0.074 | 41 | 384 |
| off | high/night/low_sun | 34 | 0.482 (5.455) | 25.2 (234.4) | 0.283 | 57 | 17 |
| off | low/day/high_sun | 2136 | 0.062 (0.174) | 10.3 (25.0) | 0.049 | 24 | 304 |
| off | low/night/high_sun | 1626 | 0.072 (0.474) | 10.9 (28.9) | 0.054 | 26 | 402 |
| off | mid/day/high_sun | 1351 | 0.075 (0.362) | 12.9 (47.0) | 0.072 | 23 | 304 |
| off | mid/day/low_sun | 20 | 0.125 (0.566) | 19.2 (53.1) | 0.090 | 12 | 7 |
| off | mid/night/high_sun | 1083 | 0.121 (0.843) | 16.7 (78.4) | 0.079 | 27 | 476 |
| off | mid/night/low_sun | 51 | 0.167 (1.994) | 15.2 (134.5) | 0.095 | 34 | 28 |
| on | high/day/high_sun | 593 | 0.107 (0.672) | 23.3 (85.7) | 0.100 | 60 | 196 |
| on | high/day/low_sun | 41 | 0.136 (2.414) | 30.2 (198.2) | 0.129 | 76 | 13 |
| on | high/night/high_sun | 855 | 0.122 (2.479) | 23.2 (212.6) | 0.105 | 62 | 375 |
| on | high/night/low_sun | 34 | 0.263 (8.092) | 69.4 (294.4) | 0.248 | 90 | 18 |
| on | low/day/high_sun | 2136 | 0.063 (0.200) | 10.9 (30.5) | 0.051 | 0 | 316 |
| on | low/night/high_sun | 1626 | 0.089 (0.468) | 13.3 (31.6) | 0.068 | 17 | 426 |
| on | mid/day/high_sun | 1351 | 0.094 (0.234) | 16.8 (41.4) | 0.091 | 7 | 299 |
| on | mid/day/low_sun | 20 | 0.098 (0.528) | 21.9 (79.7) | 0.098 | 38 | 8 |
| on | mid/night/high_sun | 1083 | 0.121 (0.852) | 18.0 (67.8) | 0.100 | 28 | 434 |
| on | mid/night/low_sun | 51 | 0.425 (2.838) | 31.2 (168.6) | 0.166 | 56 | 33 |

By quality index:

| Clutter | Group | N | SCI-05 median (p90) | SCI-07 km median (p90) | Span N_e error median | Coverage gap km median | SCI-05 fails |
|---|---|---|---|---|---|---|---|
| off | 4 | 2989 | 0.068 (0.239) | 11.1 (27.6) | 0.058 | 23 | 540 |
| off | 5 | 1937 | 0.077 (0.395) | 12.3 (42.4) | 0.062 | 25 | 515 |
| off | 6 | 1491 | 0.090 (0.927) | 15.1 (100.8) | 0.070 | 33 | 535 |
| off | 7 | 826 | 0.091 (0.932) | 17.0 (114.4) | 0.071 | 35 | 310 |
| off | 8 | 447 | 0.097 (1.263) | 17.7 (130.3) | 0.071 | 33 | 193 |
| off | 9 | 100 | 0.312 (1.456) | 26.7 (128.1) | 0.093 | 31 | 58 |
| on | 4 | 2989 | 0.078 (0.225) | 12.5 (29.1) | 0.068 | 3 | 523 |
| on | 5 | 1937 | 0.095 (0.324) | 15.8 (44.6) | 0.083 | 17 | 547 |
| on | 6 | 1491 | 0.104 (1.274) | 18.3 (129.9) | 0.093 | 31 | 531 |
| on | 7 | 826 | 0.109 (1.430) | 21.4 (167.4) | 0.100 | 46 | 296 |
| on | 8 | 447 | 0.112 (2.051) | 20.9 (170.3) | 0.101 | 48 | 174 |
| on | 9 | 100 | 0.127 (1.640) | 20.9 (130.8) | 0.090 | 46 | 47 |

---

Run 2026-11-05T13:25:05+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study isis1_retrieval` from the same commit and the same configuration reproduces it.
