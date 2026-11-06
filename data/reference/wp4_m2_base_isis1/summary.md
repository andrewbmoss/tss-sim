# isis1_retrieval

WP4 Phase 5: every usable ISIS-1 topside profile through the P4 chain, clutter on and off

Status: **implemented**.  Work packages: WP4, WP8.  Claims: C1, C3.  Requirements: SCI-05, SCI-07, INS-10.

7790 ISIS-1 records (iquality [4, 5, 6, 7, 8, 9]), clutter off, on; 4.27 h wall for this session (0 rows resumed from the checkpoint), 72.56 CPU h in all, on 17 workers.

| Clutter | N (failed) | SCI-05 median (p90, worst) | SCI-05 verdicts | SCI-07 km median (p90, worst) | SCI-07 verdicts | Span N_e error median (p90) | Span reached median (p10) | Coverage gap km median (p90) | O completeness median |
|---|---|---|---|---|---|---|---|---|---|
| off | 7790 (0) | 0.078 (0.571, 104.687) | {'fail': 2151, 'goal': 4680, 'not_evaluated': 5, 'threshold': 954} | 13.0 (56.0, 595.2) | {'fail': 2315, 'goal': 2879, 'not_evaluated': 5, 'threshold': 2591} | 0.064 (0.295) | 1.00 (0.90) | 27 (88) | 0.43 |
| on | 7790 (0) | 0.092 (0.510, 1404.083) | {'fail': 2118, 'goal': 4278, 'not_evaluated': 71, 'threshold': 1323} | 15.5 (51.1, 755.3) | {'fail': 2716, 'goal': 1977, 'not_evaluated': 96, 'threshold': 3001} | 0.082 (0.276) | 1.00 (0.81) | 18 (126) | 0.39 |

M2 rule for the measured set: median at threshold or better and at least 80% of runs at threshold or better; the measured span is the gate, the full profile (Chapman cap below the scaled heights) is shown beside it.

| Clutter | Scope | Requirement | Median | At threshold or better | Pass |
|---|---|---|---|---|---|
| off | span | SCI-05 | 0.074 | 73.1% | no |
| off | span | SCI-07 | 12.0 km | 72.8% | no |
| off | full | SCI-05 | 0.078 | 72.3% | no |
| off | full | SCI-07 | 13.0 km | 70.2% | no |
| off | **gate (span, both)** | | | | **fail** |
| on | span | SCI-05 | 0.088 | 73.1% | no |
| on | span | SCI-07 | 15.1 km | 65.3% | no |
| on | full | SCI-05 | 0.092 | 71.9% | no |
| on | full | SCI-07 | 15.5 km | 63.9% | no |
| on | **gate (span, both)** | | | | **fail** |

By plasma frequency at the spacecraft:

| Clutter | Group | N | SCI-05 median (p90) | SCI-05 span median (p90) | SCI-07 km median (p90) | SCI-07 span km median | Span N_e error median | Coverage gap km median | SCI-05 fails | SCI-05 span fails |
|---|---|---|---|---|---|---|---|---|---|---|
| off | 0-1.5 MHz | 2295 | 0.142 (1.320) | 0.141 (1.266) | 19.3 (138.6) | 16.2 | 0.081 | 32 | 1118 | 1103 |
| off | 1.5-2.5 MHz | 2570 | 0.086 (0.459) | 0.084 (0.427) | 14.2 (53.7) | 13.4 | 0.071 | 27 | 724 | 701 |
| off | 2.5-4 MHz | 1795 | 0.071 (0.196) | 0.068 (0.185) | 11.8 (26.1) | 11.4 | 0.063 | 24 | 282 | 263 |
| off | >= 4 MHz | 1130 | 0.022 (0.078) | 0.020 (0.073) | 4.2 (13.0) | 3.5 | 0.019 | 22 | 27 | 20 |
| on | 0-1.5 MHz | 2295 | 0.164 (2.065) | 0.158 (2.001) | 24.1 (188.5) | 21.8 | 0.115 | 34 | 1184 | 1171 |
| on | 1.5-2.5 MHz | 2570 | 0.102 (0.279) | 0.099 (0.262) | 17.7 (40.1) | 17.6 | 0.092 | 21 | 705 | 656 |
| on | 2.5-4 MHz | 1795 | 0.080 (0.162) | 0.076 (0.153) | 12.9 (25.6) | 13.0 | 0.073 | 11 | 213 | 190 |
| on | >= 4 MHz | 1130 | 0.035 (0.071) | 0.032 (0.068) | 6.4 (12.9) | 5.7 | 0.031 | 0 | 16 | 11 |

By stratum:

| Clutter | Group | N | SCI-05 median (p90) | SCI-05 span median (p90) | SCI-07 km median (p90) | SCI-07 span km median | Span N_e error median | Coverage gap km median | SCI-05 fails | SCI-05 span fails |
|---|---|---|---|---|---|---|---|---|---|---|
| off | high/day/high_sun | 593 | 0.082 (0.771) | 0.078 (0.736) | 17.6 (116.8) | 16.2 | 0.071 | 43 | 210 | 209 |
| off | high/day/low_sun | 41 | 0.145 (29.434) | 0.145 (26.008) | 37.2 (410.9) | 26.7 | 0.156 | 70 | 19 | 19 |
| off | high/night/high_sun | 855 | 0.112 (1.204) | 0.109 (1.150) | 20.9 (128.3) | 14.9 | 0.074 | 41 | 384 | 378 |
| off | high/night/low_sun | 34 | 0.482 (5.455) | 0.419 (5.363) | 25.2 (234.4) | 23.0 | 0.283 | 57 | 17 | 17 |
| off | low/day/high_sun | 2136 | 0.062 (0.174) | 0.057 (0.170) | 10.3 (25.0) | 10.0 | 0.049 | 24 | 304 | 287 |
| off | low/night/high_sun | 1626 | 0.072 (0.474) | 0.071 (0.402) | 10.9 (28.9) | 9.9 | 0.054 | 26 | 402 | 382 |
| off | mid/day/high_sun | 1351 | 0.075 (0.362) | 0.074 (0.358) | 12.9 (47.0) | 12.7 | 0.072 | 23 | 304 | 301 |
| off | mid/day/low_sun | 20 | 0.125 (0.566) | 0.114 (0.575) | 19.2 (53.1) | 18.1 | 0.090 | 12 | 7 | 7 |
| off | mid/night/high_sun | 1083 | 0.121 (0.843) | 0.115 (0.824) | 16.7 (78.4) | 13.8 | 0.079 | 27 | 476 | 462 |
| off | mid/night/low_sun | 51 | 0.167 (1.994) | 0.146 (2.018) | 15.2 (134.5) | 14.6 | 0.095 | 34 | 28 | 25 |
| on | high/day/high_sun | 593 | 0.107 (0.672) | 0.103 (0.666) | 23.3 (85.7) | 22.1 | 0.100 | 60 | 196 | 193 |
| on | high/day/low_sun | 41 | 0.136 (2.414) | 0.138 (2.289) | 30.2 (198.2) | 29.9 | 0.129 | 76 | 13 | 12 |
| on | high/night/high_sun | 855 | 0.122 (2.479) | 0.118 (2.344) | 23.2 (212.6) | 19.7 | 0.105 | 62 | 375 | 371 |
| on | high/night/low_sun | 34 | 0.263 (8.092) | 0.266 (7.208) | 69.4 (294.4) | 58.4 | 0.248 | 90 | 18 | 17 |
| on | low/day/high_sun | 2136 | 0.063 (0.200) | 0.058 (0.189) | 10.9 (30.5) | 10.5 | 0.051 | 0 | 316 | 296 |
| on | low/night/high_sun | 1626 | 0.089 (0.468) | 0.085 (0.413) | 13.3 (31.6) | 13.2 | 0.068 | 17 | 426 | 384 |
| on | mid/day/high_sun | 1351 | 0.094 (0.234) | 0.091 (0.222) | 16.8 (41.4) | 16.9 | 0.091 | 7 | 299 | 282 |
| on | mid/day/low_sun | 20 | 0.098 (0.528) | 0.096 (0.522) | 21.9 (79.7) | 22.3 | 0.098 | 38 | 8 | 8 |
| on | mid/night/high_sun | 1083 | 0.121 (0.852) | 0.117 (0.814) | 18.0 (67.8) | 16.6 | 0.100 | 28 | 434 | 431 |
| on | mid/night/low_sun | 51 | 0.425 (2.838) | 0.402 (2.729) | 31.2 (168.6) | 23.8 | 0.166 | 56 | 33 | 34 |

By quality index:

| Clutter | Group | N | SCI-05 median (p90) | SCI-05 span median (p90) | SCI-07 km median (p90) | SCI-07 span km median | Span N_e error median | Coverage gap km median | SCI-05 fails | SCI-05 span fails |
|---|---|---|---|---|---|---|---|---|---|---|
| off | 4 | 2989 | 0.068 (0.239) | 0.065 (0.226) | 11.1 (27.6) | 10.5 | 0.058 | 23 | 540 | 510 |
| off | 5 | 1937 | 0.077 (0.395) | 0.074 (0.372) | 12.3 (42.4) | 11.4 | 0.062 | 25 | 515 | 499 |
| off | 6 | 1491 | 0.090 (0.927) | 0.088 (0.882) | 15.1 (100.8) | 13.8 | 0.070 | 33 | 535 | 524 |
| off | 7 | 826 | 0.091 (0.932) | 0.087 (0.915) | 17.0 (114.4) | 14.7 | 0.071 | 35 | 310 | 306 |
| off | 8 | 447 | 0.097 (1.263) | 0.093 (1.227) | 17.7 (130.3) | 14.4 | 0.071 | 33 | 193 | 190 |
| off | 9 | 100 | 0.312 (1.456) | 0.285 (1.321) | 26.7 (128.1) | 20.5 | 0.093 | 31 | 58 | 58 |
| on | 4 | 2989 | 0.078 (0.225) | 0.074 (0.212) | 12.5 (29.1) | 12.3 | 0.068 | 3 | 523 | 472 |
| on | 5 | 1937 | 0.095 (0.324) | 0.093 (0.303) | 15.8 (44.6) | 15.4 | 0.083 | 17 | 547 | 525 |
| on | 6 | 1491 | 0.104 (1.274) | 0.100 (1.238) | 18.3 (129.9) | 17.6 | 0.093 | 31 | 531 | 522 |
| on | 7 | 826 | 0.109 (1.430) | 0.107 (1.379) | 21.4 (167.4) | 19.6 | 0.100 | 46 | 296 | 289 |
| on | 8 | 447 | 0.112 (2.051) | 0.108 (1.974) | 20.9 (170.3) | 19.4 | 0.101 | 48 | 174 | 173 |
| on | 9 | 100 | 0.127 (1.640) | 0.122 (1.509) | 20.9 (130.8) | 19.0 | 0.090 | 46 | 47 | 47 |

---

Run 2026-11-05T18:29:18+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study isis1_retrieval` from the same commit and the same configuration reproduces it.
