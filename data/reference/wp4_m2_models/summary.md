# IRI-2020 against ISIS-1, and the coverage of the data (WP4 P6)

7790 usable ISIS-1 records (quality 4-9, 500-1100 km), the records of `isis1_retrieval`.  Each IRI-2020 profile (epoch indices, foF2 storm model on) is scored against the measured profile over the measured span on SCI-05's metric (rms fractional N_e) and SCI-07's (rms height at matched density), with the Rev G limits (0.15 / 0.10; 20 / 10 km).  The plasma frequency at the spacecraft is a ratio model/measured in dB (20 log10).

## overall

| Group | N | Topside | SCI-05 metric median (p90) | at threshold or better | SCI-07 metric km median (p90) | at threshold or better | f_N at spacecraft dB median [p5, p95] |
|---|---|---|---|---|---|---|---|
| all | 7790 | NeQuick | 0.417 (1.241) | 12% | 70.2 (162.0) | 8% | 0.58 [-4.09, 4.81] |
| all | 7790 | IRI2001 | 0.524 (2.465) | 6% | 85.9 (225.1) | 4% | 2.43 [-2.92, 9.23] |
| all | 7790 | IRI2001corr | 0.446 (1.792) | 8% | 73.4 (193.4) | 6% | 0.87 [-3.91, 7.34] |

## by fn band

| Group | N | Topside | SCI-05 metric median (p90) | at threshold or better | SCI-07 metric km median (p90) | at threshold or better | f_N at spacecraft dB median [p5, p95] |
|---|---|---|---|---|---|---|---|
| 0-1.5 MHz | 2295 | NeQuick | 0.646 (2.217) | 3% | 96.2 (207.9) | 1% | 2.02 [-2.32, 6.74] |
| 0-1.5 MHz | 2295 | IRI2001 | 1.764 (4.987) | 0% | 177.4 (321.0) | 0% | 6.41 [2.00, 12.64] |
| 0-1.5 MHz | 2295 | IRI2001corr | 1.097 (3.751) | 1% | 132.0 (277.3) | 0% | 3.78 [-0.27, 10.61] |
| 1.5-2.5 MHz | 2570 | NeQuick | 0.393 (1.014) | 12% | 64.4 (144.7) | 7% | 0.82 [-3.79, 3.90] |
| 1.5-2.5 MHz | 2570 | IRI2001 | 0.516 (1.294) | 5% | 78.6 (152.6) | 1% | 3.01 [0.15, 5.62] |
| 1.5-2.5 MHz | 2570 | IRI2001corr | 0.369 (1.123) | 8% | 61.2 (137.4) | 4% | 1.15 [-1.79, 4.40] |
| 2.5-4 MHz | 1795 | NeQuick | 0.322 (0.740) | 20% | 50.6 (138.2) | 16% | -0.12 [-4.15, 3.30] |
| 2.5-4 MHz | 1795 | IRI2001 | 0.302 (0.619) | 14% | 49.1 (108.6) | 10% | 0.67 [-1.74, 2.89] |
| 2.5-4 MHz | 1795 | IRI2001corr | 0.312 (0.648) | 17% | 50.7 (133.4) | 14% | -0.39 [-3.03, 2.36] |
| >= 4 MHz | 1130 | NeQuick | 0.376 (0.565) | 13% | 64.9 (130.7) | 12% | -1.84 [-5.03, 1.21] |
| >= 4 MHz | 1130 | IRI2001 | 0.377 (0.583) | 9% | 65.2 (133.5) | 9% | -2.06 [-5.08, 0.69] |
| >= 4 MHz | 1130 | IRI2001corr | 0.409 (0.622) | 9% | 69.2 (135.9) | 8% | -2.76 [-6.00, 0.31] |

## by solar

| Group | N | Topside | SCI-05 metric median (p90) | at threshold or better | SCI-07 metric km median (p90) | at threshold or better | f_N at spacecraft dB median [p5, p95] |
|---|---|---|---|---|---|---|---|
| solar_min | 95 | NeQuick | 0.365 (3.378) | 8% | 52.9 (152.7) | 8% | 0.63 [-2.74, 10.17] |
| solar_min | 95 | IRI2001 | 0.370 (7.729) | 5% | 54.0 (203.5) | 5% | 1.25 [-1.39, 13.80] |
| solar_min | 95 | IRI2001corr | 0.374 (5.323) | 5% | 52.5 (163.2) | 4% | 0.53 [-2.46, 11.73] |
| intermediate | 445 | NeQuick | 0.414 (0.843) | 4% | 82.9 (163.0) | 3% | -1.50 [-4.79, 4.71] |
| intermediate | 445 | IRI2001 | 0.431 (1.324) | 3% | 87.9 (170.2) | 1% | 1.05 [-3.89, 8.24] |
| intermediate | 445 | IRI2001corr | 0.433 (0.780) | 4% | 84.4 (169.7) | 2% | -1.12 [-5.66, 4.58] |
| solar_max | 7250 | NeQuick | 0.418 (1.253) | 12% | 69.5 (162.0) | 9% | 0.67 [-4.02, 4.76] |
| solar_max | 7250 | IRI2001 | 0.537 (2.511) | 7% | 85.9 (228.4) | 4% | 2.58 [-2.82, 9.23] |
| solar_max | 7250 | IRI2001corr | 0.447 (1.852) | 8% | 72.8 (195.7) | 6% | 0.97 [-3.80, 7.43] |

## by stratum

| Group | N | Topside | SCI-05 metric median (p90) | at threshold or better | SCI-07 metric km median (p90) | at threshold or better | f_N at spacecraft dB median [p5, p95] |
|---|---|---|---|---|---|---|---|
| high/day/high_sun | 593 | NeQuick | 0.411 (0.766) | 7% | 96.7 (192.8) | 3% | -1.88 [-5.52, 3.31] |
| high/day/high_sun | 593 | IRI2001 | 0.800 (3.920) | 2% | 139.3 (298.4) | 0% | 4.17 [-0.51, 11.80] |
| high/day/high_sun | 593 | IRI2001corr | 0.679 (3.258) | 4% | 126.2 (296.6) | 1% | 3.28 [-0.69, 10.49] |
| high/day/low_sun | 41 | NeQuick | 0.452 (5.651) | 0% | 133.4 (275.9) | 0% | -3.07 [-4.70, 10.03] |
| high/day/low_sun | 41 | IRI2001 | 0.466 (10.858) | 0% | 122.1 (360.9) | 0% | -0.29 [-2.85, 14.05] |
| high/day/low_sun | 41 | IRI2001corr | 0.406 (8.544) | 0% | 125.0 (336.5) | 0% | -2.12 [-3.27, 11.76] |
| high/night/high_sun | 855 | NeQuick | 0.455 (1.234) | 5% | 87.3 (192.7) | 2% | -0.01 [-4.87, 5.62] |
| high/night/high_sun | 855 | IRI2001 | 1.413 (6.607) | 1% | 150.2 (345.1) | 0% | 6.06 [-0.26, 13.99] |
| high/night/high_sun | 855 | IRI2001corr | 0.960 (4.589) | 2% | 129.8 (317.5) | 0% | 4.42 [-0.91, 11.64] |
| high/night/low_sun | 34 | NeQuick | 0.517 (2.798) | 9% | 102.2 (172.3) | 9% | 0.93 [-3.58, 11.84] |
| high/night/low_sun | 34 | IRI2001 | 0.594 (3.978) | 0% | 103.2 (174.2) | 0% | 1.94 [-2.66, 14.63] |
| high/night/low_sun | 34 | IRI2001corr | 0.510 (3.036) | 3% | 96.7 (154.6) | 0% | 1.28 [-2.82, 12.08] |
| low/day/high_sun | 2136 | NeQuick | 0.414 (0.954) | 12% | 72.0 (150.1) | 9% | -0.02 [-4.26, 3.21] |
| low/day/high_sun | 2136 | IRI2001 | 0.420 (1.205) | 6% | 75.4 (171.9) | 5% | 0.85 [-4.28, 5.39] |
| low/day/high_sun | 2136 | IRI2001corr | 0.408 (0.722) | 9% | 66.6 (145.0) | 8% | -0.81 [-5.26, 2.21] |
| low/night/high_sun | 1626 | NeQuick | 0.477 (1.687) | 7% | 66.6 (168.6) | 6% | 1.85 [-2.39, 5.73] |
| low/night/high_sun | 1626 | IRI2001 | 0.855 (2.871) | 3% | 92.0 (228.5) | 2% | 4.35 [-0.24, 8.90] |
| low/night/high_sun | 1626 | IRI2001corr | 0.482 (1.829) | 5% | 68.9 (184.4) | 4% | 2.00 [-2.68, 5.92] |
| mid/day/high_sun | 1351 | NeQuick | 0.292 (0.850) | 23% | 49.5 (121.1) | 16% | 0.28 [-2.58, 3.40] |
| mid/day/high_sun | 1351 | IRI2001 | 0.324 (0.788) | 15% | 55.9 (125.9) | 8% | 0.45 [-3.23, 4.19] |
| mid/day/high_sun | 1351 | IRI2001corr | 0.303 (0.700) | 15% | 50.5 (115.4) | 10% | -0.64 [-3.51, 2.58] |
| mid/day/low_sun | 20 | NeQuick | 0.290 (0.448) | 5% | 59.0 (86.9) | 5% | -1.06 [-3.34, 0.76] |
| mid/day/low_sun | 20 | IRI2001 | 0.302 (0.454) | 10% | 60.6 (107.5) | 0% | -0.43 [-2.53, 1.27] |
| mid/day/low_sun | 20 | IRI2001corr | 0.346 (0.466) | 10% | 63.4 (105.6) | 5% | -1.34 [-3.67, -0.27] |
| mid/night/high_sun | 1083 | NeQuick | 0.514 (1.783) | 12% | 73.4 (172.9) | 9% | 1.68 [-2.16, 5.64] |
| mid/night/high_sun | 1083 | IRI2001 | 0.821 (2.534) | 9% | 99.0 (227.6) | 5% | 3.42 [-0.47, 8.06] |
| mid/night/high_sun | 1083 | IRI2001corr | 0.627 (2.070) | 10% | 83.9 (204.0) | 6% | 2.17 [-1.69, 6.59] |
| mid/night/low_sun | 51 | NeQuick | 0.299 (0.723) | 8% | 38.7 (94.8) | 8% | -0.49 [-2.67, 2.24] |
| mid/night/low_sun | 51 | IRI2001 | 0.338 (0.988) | 8% | 44.1 (112.4) | 10% | 0.43 [-1.45, 3.14] |
| mid/night/low_sun | 51 | IRI2001corr | 0.339 (0.955) | 4% | 42.1 (110.3) | 6% | -0.45 [-2.34, 2.24] |

## by quality

| Group | N | Topside | SCI-05 metric median (p90) | at threshold or better | SCI-07 metric km median (p90) | at threshold or better | f_N at spacecraft dB median [p5, p95] |
|---|---|---|---|---|---|---|---|
| 4 | 2989 | NeQuick | 0.360 (0.945) | 16% | 56.4 (132.8) | 13% | 0.48 [-3.72, 3.66] |
| 4 | 2989 | IRI2001 | 0.406 (1.199) | 9% | 64.8 (147.0) | 7% | 1.02 [-3.90, 5.78] |
| 4 | 2989 | IRI2001corr | 0.363 (0.864) | 12% | 56.8 (135.8) | 9% | -0.18 [-4.65, 3.54] |
| 5 | 1937 | NeQuick | 0.437 (1.399) | 12% | 69.9 (159.4) | 8% | 1.14 [-3.50, 4.97] |
| 5 | 1937 | IRI2001 | 0.577 (2.143) | 7% | 88.6 (202.5) | 3% | 2.75 [-2.66, 7.86] |
| 5 | 1937 | IRI2001corr | 0.438 (1.491) | 9% | 70.1 (167.8) | 6% | 0.87 [-3.76, 5.41] |
| 6 | 1491 | NeQuick | 0.458 (1.572) | 6% | 81.9 (178.0) | 4% | 0.79 [-3.85, 5.33] |
| 6 | 1491 | IRI2001 | 0.832 (3.541) | 3% | 110.3 (264.8) | 1% | 3.90 [-1.69, 10.90] |
| 6 | 1491 | IRI2001corr | 0.544 (2.525) | 4% | 91.4 (228.8) | 3% | 1.96 [-2.91, 8.41] |
| 7 | 826 | NeQuick | 0.449 (1.219) | 7% | 92.7 (193.6) | 4% | -0.61 [-5.05, 5.37] |
| 7 | 826 | IRI2001 | 0.934 (4.112) | 3% | 129.1 (307.7) | 1% | 4.49 [-0.68, 11.91] |
| 7 | 826 | IRI2001corr | 0.696 (3.136) | 3% | 113.1 (270.2) | 1% | 3.11 [-1.58, 10.06] |
| 8 | 447 | NeQuick | 0.464 (1.700) | 5% | 85.6 (217.5) | 2% | -0.56 [-5.57, 7.91] |
| 8 | 447 | IRI2001 | 0.911 (5.699) | 4% | 136.0 (337.1) | 2% | 4.43 [-1.06, 13.89] |
| 8 | 447 | IRI2001corr | 0.715 (4.305) | 5% | 118.6 (296.4) | 2% | 3.11 [-2.17, 12.04] |
| 9 | 100 | NeQuick | 0.445 (1.869) | 8% | 93.2 (207.2) | 5% | -0.70 [-4.41, 5.84] |
| 9 | 100 | IRI2001 | 0.615 (3.897) | 1% | 95.5 (255.1) | 0% | 3.64 [0.31, 11.23] |
| 9 | 100 | IRI2001corr | 0.620 (3.121) | 1% | 107.0 (228.0) | 0% | 2.70 [-1.86, 10.94] |

## Coverage: archive (38953 records, 1969-1980, spacecraft to 3550 km)

| Solar activity | Records |
|---|---|
| solar_min | 577 |
| intermediate | 4724 |
| solar_max | 33652 |

| Spacecraft height | Records |
|---|---|
| 0-500 km | 0 |
| 500-700 km | 2949 |
| 700-900 km | 2538 |
| 900-1100 km | 2303 |
| 1100-1400 km | 3355 |
| 1400-2000 km | 6087 |
| 2000-4000 km | 21721 |

| Quality index | Records |
|---|---|
| -1 | 8 |
| 0 | 2 |
| 1 | 1 |
| 2 | 1 |
| 3 | 1 |
| 4 | 4280 |
| 5 | 7477 |
| 6 | 14171 |
| 7 | 7640 |
| 8 | 4591 |
| 9 | 781 |

Empty latitude x local-time cells (15 deg x 3 h, of 96): solar_min 70; intermediate 20; solar_max 0

## Coverage: usable (7790 records, 1969-1979, spacecraft to 1100 km)

| Solar activity | Records |
|---|---|
| solar_min | 95 |
| intermediate | 445 |
| solar_max | 7250 |

| Spacecraft height | Records |
|---|---|
| 0-500 km | 0 |
| 500-700 km | 2949 |
| 700-900 km | 2538 |
| 900-1100 km | 2301 |
| 1100-1400 km | 2 |
| 1400-2000 km | 0 |
| 2000-4000 km | 0 |

| Quality index | Records |
|---|---|
| 4 | 2989 |
| 5 | 1937 |
| 6 | 1491 |
| 7 | 826 |
| 8 | 447 |
| 9 | 100 |

Empty latitude x local-time cells (15 deg x 3 h, of 96): solar_min 85; intermediate 52; solar_max 11

`coverage.png`: the usable records by latitude and local time, per solar-activity class.
