# isis1_retrieval

WP4 Phase 5: every usable ISIS-1 topside profile through the P4 chain, clutter on and off

Status: **implemented**.  Work packages: WP4, WP8.  Claims: C1, C3.  Requirements: SCI-05, SCI-07, INS-10.

7790 ISIS-1 records (iquality [4, 5, 6, 7, 8, 9]), clutter off, on; 4.21 h wall for this session (0 rows resumed from the checkpoint), 71.47 CPU h in all, on 17 workers.

| Clutter | N (failed) | SCI-05 median (p90, worst) | SCI-05 verdicts | SCI-07 km median (p90, worst) | SCI-07 verdicts | Span N_e error median (p90) | Span reached median (p10) | Coverage gap km median (p90) | O completeness median |
|---|---|---|---|---|---|---|---|---|---|
| off | 7790 (0) | 0.050 (0.364, 115.562) | {'fail': 1480, 'goal': 5718, 'threshold': 592} | 8.3 (31.0, 690.9) | {'fail': 1298, 'goal': 4566, 'threshold': 1926} | 0.038 (0.127) | 1.00 (0.92) | 24 (81) | 0.47 |
| on | 7790 (0) | 0.039 (0.315, 101.661) | {'fail': 1263, 'goal': 6155, 'not_evaluated': 9, 'threshold': 363} | 6.7 (27.3, 592.3) | {'fail': 1085, 'goal': 5318, 'not_evaluated': 9, 'threshold': 1378} | 0.030 (0.101) | 1.00 (0.86) | 14 (106) | 0.40 |

M2 rule for the measured set: median at threshold or better and at least 80% of runs at threshold or better; the measured span is the gate, the full profile (Chapman cap below the scaled heights) is shown beside it.

| Clutter | Scope | Requirement | Median | At threshold or better | Pass |
|---|---|---|---|---|---|
| off | span | SCI-05 | 0.048 | 81.8% | yes |
| off | span | SCI-07 | 7.2 km | 87.6% | yes |
| off | full | SCI-05 | 0.050 | 81.0% | yes |
| off | full | SCI-07 | 8.3 km | 83.3% | yes |
| off | **gate (span, both)** | | | | **pass** |
| on | span | SCI-05 | 0.037 | 84.6% | yes |
| on | span | SCI-07 | 5.7 km | 89.3% | yes |
| on | full | SCI-05 | 0.039 | 83.7% | yes |
| on | full | SCI-07 | 6.7 km | 86.0% | yes |
| on | **gate (span, both)** | | | | **pass** |

By plasma frequency at the spacecraft:

| Clutter | Group | N | SCI-05 median (p90) | SCI-05 span median (p90) | SCI-07 km median (p90) | SCI-07 span km median | Span N_e error median | Coverage gap km median | SCI-05 fails | SCI-05 span fails |
|---|---|---|---|---|---|---|---|---|---|---|
| off | 0-1.5 MHz | 2295 | 0.107 (1.200) | 0.102 (1.154) | 14.5 (108.0) | 11.8 | 0.060 | 29 | 976 | 943 |
| off | 1.5-2.5 MHz | 2570 | 0.057 (0.215) | 0.055 (0.202) | 9.3 (21.8) | 8.3 | 0.045 | 23 | 387 | 376 |
| off | 2.5-4 MHz | 1795 | 0.039 (0.114) | 0.036 (0.109) | 6.4 (15.1) | 5.4 | 0.030 | 21 | 100 | 88 |
| off | >= 4 MHz | 1130 | 0.015 (0.041) | 0.013 (0.036) | 3.0 (7.7) | 2.2 | 0.011 | 21 | 17 | 9 |
| on | 0-1.5 MHz | 2295 | 0.083 (1.240) | 0.079 (1.218) | 12.6 (103.3) | 10.0 | 0.044 | 28 | 859 | 839 |
| on | 1.5-2.5 MHz | 2570 | 0.044 (0.194) | 0.042 (0.167) | 7.4 (19.7) | 6.6 | 0.035 | 17 | 322 | 288 |
| on | 2.5-4 MHz | 1795 | 0.029 (0.069) | 0.028 (0.065) | 4.9 (11.2) | 4.2 | 0.023 | 7 | 71 | 58 |
| on | >= 4 MHz | 1130 | 0.018 (0.033) | 0.015 (0.030) | 3.7 (7.7) | 2.6 | 0.012 | 0 | 11 | 8 |

By stratum:

| Clutter | Group | N | SCI-05 median (p90) | SCI-05 span median (p90) | SCI-07 km median (p90) | SCI-07 span km median | Span N_e error median | Coverage gap km median | SCI-05 fails | SCI-05 span fails |
|---|---|---|---|---|---|---|---|---|---|---|
| off | high/day/high_sun | 593 | 0.050 (0.428) | 0.048 (0.423) | 10.8 (44.1) | 9.6 | 0.043 | 31 | 122 | 120 |
| off | high/day/low_sun | 41 | 0.125 (40.819) | 0.127 (40.819) | 24.2 (435.3) | 21.2 | 0.102 | 84 | 19 | 19 |
| off | high/night/high_sun | 855 | 0.066 (0.941) | 0.060 (0.893) | 11.3 (63.7) | 8.9 | 0.045 | 31 | 263 | 256 |
| off | high/night/low_sun | 34 | 0.075 (5.611) | 0.076 (5.257) | 16.9 (241.5) | 13.8 | 0.068 | 54 | 14 | 14 |
| off | low/day/high_sun | 2136 | 0.037 (0.144) | 0.037 (0.137) | 6.4 (17.7) | 5.5 | 0.027 | 22 | 190 | 182 |
| off | low/night/high_sun | 1626 | 0.051 (0.401) | 0.050 (0.342) | 7.7 (22.7) | 6.7 | 0.037 | 24 | 325 | 307 |
| off | mid/day/high_sun | 1351 | 0.045 (0.179) | 0.043 (0.171) | 7.7 (20.6) | 6.8 | 0.038 | 20 | 160 | 154 |
| off | mid/day/low_sun | 20 | 0.071 (0.323) | 0.061 (0.293) | 10.6 (33.8) | 8.7 | 0.060 | 7 | 3 | 3 |
| off | mid/night/high_sun | 1083 | 0.067 (0.682) | 0.064 (0.629) | 10.0 (49.7) | 8.0 | 0.045 | 25 | 359 | 340 |
| off | mid/night/low_sun | 51 | 0.135 (0.690) | 0.110 (0.801) | 13.3 (55.6) | 11.8 | 0.069 | 21 | 25 | 21 |
| on | high/day/high_sun | 593 | 0.048 (0.389) | 0.045 (0.394) | 9.9 (42.4) | 8.5 | 0.037 | 52 | 106 | 101 |
| on | high/day/low_sun | 41 | 0.093 (41.558) | 0.094 (40.382) | 21.0 (435.4) | 15.2 | 0.070 | 79 | 16 | 17 |
| on | high/night/high_sun | 855 | 0.056 (1.892) | 0.053 (1.841) | 10.2 (165.6) | 8.4 | 0.040 | 51 | 235 | 231 |
| on | high/night/low_sun | 34 | 0.123 (7.563) | 0.096 (6.864) | 18.2 (314.7) | 16.2 | 0.068 | 66 | 14 | 14 |
| on | low/day/high_sun | 2136 | 0.024 (0.104) | 0.023 (0.103) | 4.7 (13.8) | 3.8 | 0.019 | 0 | 129 | 116 |
| on | low/night/high_sun | 1626 | 0.038 (0.379) | 0.037 (0.300) | 5.9 (21.5) | 5.1 | 0.028 | 14 | 310 | 284 |
| on | mid/day/high_sun | 1351 | 0.037 (0.114) | 0.036 (0.111) | 6.6 (17.4) | 5.6 | 0.031 | 3 | 109 | 105 |
| on | mid/day/low_sun | 20 | 0.049 (0.470) | 0.050 (0.465) | 10.2 (48.8) | 9.5 | 0.049 | 36 | 4 | 4 |
| on | mid/night/high_sun | 1083 | 0.055 (0.579) | 0.052 (0.553) | 8.4 (44.7) | 6.7 | 0.037 | 23 | 321 | 302 |
| on | mid/night/low_sun | 51 | 0.102 (1.091) | 0.105 (1.140) | 10.5 (71.0) | 10.3 | 0.056 | 42 | 19 | 19 |

By quality index:

| Clutter | Group | N | SCI-05 median (p90) | SCI-05 span median (p90) | SCI-07 km median (p90) | SCI-07 span km median | Span N_e error median | Coverage gap km median | SCI-05 fails | SCI-05 span fails |
|---|---|---|---|---|---|---|---|---|---|---|
| off | 4 | 2989 | 0.040 (0.167) | 0.039 (0.154) | 6.7 (18.8) | 5.6 | 0.030 | 22 | 339 | 318 |
| off | 5 | 1937 | 0.053 (0.278) | 0.051 (0.258) | 8.3 (25.2) | 7.2 | 0.039 | 22 | 351 | 335 |
| off | 6 | 1491 | 0.060 (0.642) | 0.059 (0.634) | 10.0 (51.1) | 8.7 | 0.045 | 28 | 411 | 390 |
| off | 7 | 826 | 0.057 (0.655) | 0.053 (0.614) | 10.3 (44.2) | 8.8 | 0.042 | 29 | 203 | 197 |
| off | 8 | 447 | 0.066 (1.127) | 0.060 (1.051) | 11.1 (71.7) | 8.9 | 0.044 | 27 | 140 | 141 |
| off | 9 | 100 | 0.069 (1.259) | 0.059 (1.165) | 11.6 (70.4) | 8.5 | 0.046 | 23 | 36 | 35 |
| on | 4 | 2989 | 0.031 (0.141) | 0.029 (0.127) | 5.2 (15.3) | 4.4 | 0.024 | 0 | 284 | 250 |
| on | 5 | 1937 | 0.040 (0.242) | 0.039 (0.214) | 6.6 (23.2) | 5.7 | 0.030 | 14 | 302 | 288 |
| on | 6 | 1491 | 0.047 (0.615) | 0.045 (0.623) | 8.2 (49.9) | 6.9 | 0.034 | 26 | 339 | 331 |
| on | 7 | 826 | 0.051 (0.719) | 0.049 (0.646) | 9.7 (53.3) | 8.1 | 0.037 | 39 | 186 | 177 |
| on | 8 | 447 | 0.052 (1.566) | 0.049 (1.593) | 9.2 (114.8) | 7.5 | 0.036 | 43 | 120 | 115 |
| on | 9 | 100 | 0.052 (1.200) | 0.047 (1.046) | 8.4 (84.3) | 7.6 | 0.033 | 43 | 32 | 32 |

---

Run 2026-11-06T11:49:43+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study isis1_retrieval` from the same commit and the same configuration reproduces it.
