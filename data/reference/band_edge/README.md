# Reference: the band-edge run that set INS-01a (Requirements Rev G)

Frozen on 2026-10-19 (WP0 closure, Phase 0) so that the band-edge study in this
repository (`tss study min_frequency`, plan `band_edge_reference`) can be held to
the run the requirement was set from.

| File | What it is |
|---|---|
| `minfreq_summary.csv` | 60 cases: O1/O2 x 500-900 km x solar min (R12 10) / max (R12 120) x NeQuick / IRI-2001 / IRI-2001 corrected, 20,000 samples each. Written by `Dissertation/MinFreqStudy/run_minfreq_study.py` on 2026-09-07 10:30 (MDT). |
| `completeness_curves.csv` | The completeness curves of the same run. |
| `run.log` | The run's console log. |
| `orbit_sample_fixture.json` | Five samples per orbit class (900 km, n = 20,000, seed 17) drawn with `MinFreqStudy/orbits.py`, against which `tss.platform.orbits` is tested. |

The governing value is O1, 900 km, solar minimum, NeQuick: 0.334 MHz at 95 per
cent coverage, the minimum of the three topsides (IRI-2001 0.477, IRI-2001
corrected 0.354).

SHA-256 of the sources at the time of freezing:

| Source | SHA-256 |
|---|---|
| `MinFreqStudy/minfreq_study.py` | 5733c47214ee19159fd6db45d52ad74f1178034759143bb720cfb0e1c6e1371d |
| `MinFreqStudy/run_minfreq_study.py` | e71dce4618be9a455818dde28572589c219d100794a42c017cd058c0ac7f29f5 |
| `MinFreqStudy/orbits.py` | 9c485cd1169e65f9d1d07b3a936a2e9c2e0c6e870ee89b5d27733bb60301a02a |
| `MinFreqStudy/iri_adapter.py` | 3766a180bd6330f58f8fa5ca6ba75ac3517457ebd281a22f1da648185a69faec |
| `MinFreqStudy/matlab/iri_batch.m` | 9378bf1eb6e7f8b600bbe5b650bc3711a7fa0f685b458cd1f083ee13e868072c |
| `MinFreqStudy/out/minfreq_summary.csv` | ef402591407a89338f39207684918ebbf26c79edb223e9a2ee3470b5890535e5 |
| `MinFreqStudy/out/completeness_curves.csv` | 1410392a77ab7ff90ea3e7c68962292d04498e837c624a8fcaa35cc7b9fd5ad4 |
| `Dissertation/TSS_SRM_RevG.xlsx` | fa7286861d19ea2abc4539839fd3be2c3255aeae32f9477869a9ab33d3716a18 |
| `Dissertation/ISISValidation/t1_results.csv` | f4f696224c3da216812977499e2fc99e3856be00eeee2f4ccf46bdc04fb13858 |

`minfreq_study.py` was last modified at 10:59 on 2026-09-07, 29 minutes after
the run; the regression test in `tests/unit/test_band_edge.py` is what shows the
arithmetic did not change.

PHaRLAP 4.7.4 index files: `ig_rz.dat` covers January 1958 to November 2026 and
`apf107.dat` ends on 2024-12-09, so an epoch-driven run is limited to 1958-01-01
to 2024-12-09 (`tss.l0_ionosphere.backends.iri_index_span`).
