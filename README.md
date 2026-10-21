# tss-sim

Modelling and simulation codebase for a low size, weight and power FMCW topside
ionospheric sounder intended for a CubeSat.  This repository supports the
doctoral dissertation of Andrew B. Moss at Utah State University.

The instrument under study transmits under 1 W across HF, 0.33 MHz to 20 MHz,
from electrically small antennas on a platform of 6U or less.  The lower band
edge of 0.33 MHz is a derived number and not a design choice; it comes out of
the minimum sounding frequency study described in `docs/work-packages.md`.  The
simulator exists to answer whether such an instrument can return usable topside
electron density profiles, and to quantify what it cannot do.

The Python package is `tss`.  MATLAB functions carry a `tss_` prefix.  The two
implementations are peers, bridged by files rather than by the MATLAB Engine
API; see `docs/matlab-python-bridge.md`.

## What the simulator is

It is a truth-in, profile-out chain of six layers.  An ionosphere is generated
at L0, echoes are propagated at L1, the instrument observes them at L2, signal
processing forms an ionogram at L3, retrieval inverts it back to a density
profile at L4, and L5 scores that profile against the L0 truth it started from
and against the machine-readable requirement set.

Because truth is known exactly at L0, every error at L5 can be attributed to a
layer.  That is the property the whole design is built around.  The drawback is
that the chain is only as trustworthy as its weakest tier, so each layer carries
more than one fidelity implementation and the verification plan in
`docs/verification.md` compares them against each other and against archived
real data.

## Architecture

```
                 campaign YAML  ──────────────┐
                 (config/campaigns/*.yaml)    │  names one tier per layer,
                                              │  a seed, and the configurations
                                              ▼
  ┌───────────────────────────────────────────────────────────────────────┐
  │                          tss.pipeline.Pipeline                        │
  │            sequencing only, no physics, one RNG per config            │
  └───────────────────────────────────────────────────────────────────────┘
        │
        │  L0  Ionosphere state            src/tss/l0_ionosphere/
        │      truth Ne, collision frequency, geomagnetic field
        │      IRI-2020, NeQuick 2, alpha-Chapman, tabulated
        ▼      ── IonosphereState ──────────────────────────────────┐
        │                                                           │
        │  L1  Propagation                 src/tss/l1_propagation/  │
        │      group delay, amplitude, absorption, O/X splitting,   │
        │      Doppler.  Analytic tier and PHaRLAP ray-traced tier. │
        ▼      ── EchoSet ───────────────────────────────┐          │
        │                                                │          │
        │  L2  Instrument                  src/tss/l2_instrument/   │
        │      sweep generation, electrically small antenna, front  │
        │      end, noise, phase noise, transmit-to-receive leakage │
        ▼      ── ReceivedSignal (+ SweepSchedule) ──┐    │          │
        │                                           │    │          │
        │  L3  Signal processing           src/tss/l3_signal/       │
        │      de-chirp, transform, bin index to apparent range,    │
        │      trace extraction, Doppler and range decoupling       │
        ▼      ── Ionogram ──────────────────┐      │    │          │
        │                                    │      │    │          │
        │  L4  Retrieval                   src/tss/l4_retrieval/    │
        │      virtual height to true height by lamination          │
        ▼      ── DensityProfile ──┐         │      │    │          │
        │                          │         │      │    │          │
        │  L5  Assessment                  src/tss/l5_assessment/   │
        │      truth comparison and requirement scoring             │
        ▼                          ▼         ▼      ▼    ▼          ▼
     Assessment              results/<study>/  manifest.json  figures/
                             per configuration, per requirement ID
```

Every tier registers itself under a layer and a name with `tss.registry`, and a
campaign selects tiers by those names.  Adding a fidelity tier does not mean
editing the pipeline.  The contracts the tiers must honour are in
`src/tss/interfaces.py` and are documented in `docs/architecture.md`.

## Quickstart

```bash
# 1.  Environment.  The pinned environment carries numpy, scipy, pyyaml,
#     pytest, matplotlib and the optional ionosphere model bindings.
conda env create -f environment.yml
conda activate tss

# 2.  Install the package in editable mode.
pip install -e .

# 3.  See which tier implementations are registered.
tss tiers

# 4.  Show the requirement set and check that every entry is traceable.
#     Exit status 2 means at least one entry cannot be used to score.
tss requirements
tss requirements --claim C3

# 5.  Resolve a campaign without running it.
tss config baseline_6u

# 6.  Run the test suite.  No test touches the network.
pytest -q
pytest -q tests/l2_instrument          # one layer
pytest -q -m "not slow"                # skip the long convergence cases
```

MATLAB users start from `matlab/startup_tss.m`, which puts the `tss_` functions
and the PHaRLAP installation on the path.  See `docs/matlab-python-bridge.md`.

Running a full campaign and writing results plus a manifest:

```bash
tss run baseline_6u --out results/baseline_6u --strict
tss verify --level V1
tss figure band_edge_cdf
```

`--strict` refuses to run from a dirty working tree.  Use it for anything
destined for the dissertation, and omit it for exploratory work, in which case
the manifest records honestly that the tree was dirty.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/tss/interfaces.py` | Layer contracts and the objects that pass between tiers.  Read this first. |
| `src/tss/registry.py` | Tier registration and lookup by layer and name. |
| `src/tss/config.py` | Campaign loading, inheritance, merging and validation. |
| `src/tss/requirements.py` | Requirement set loading and scoring. |
| `src/tss/provenance.py` | Run manifests, config hashing, git state, input checksums. |
| `src/tss/pipeline.py` | The L0 to L5 orchestrator.  Sequencing only. |
| `src/tss/cli.py` | The `tss` command line driver. |
| `src/tss/units.py`, `constants.py` | SI conventions, conversions, physical constants. |
| `src/tss/l0_ionosphere/` | Ionosphere truth models and the perturbation library. |
| `src/tss/l1_propagation/` | Echo generation, absorption, magnetoionic splitting. |
| `src/tss/l2_instrument/` | Waveform, antenna, receive chain, link budget, leakage. |
| `src/tss/l3_signal/` | De-chirp, transform, ionogram formation, novel algorithms. |
| `src/tss/l4_retrieval/` | Virtual to true height inversion and the error budget. |
| `src/tss/l5_assessment/` | Truth comparison, requirement scoring, feasibility maps. |
| `src/tss/data/` | External data fetchers and the cache under `data/external/`. |
| `src/tss/studies/` | Named studies that produce dissertation tables and figures. |
| `src/tss/verify/` | The V1 to V4 verification drivers. |
| `matlab/` | `tss_*` functions, `startup_tss.m`, MATLAB-side tiers. |
| `config/default.yaml` | Settings every campaign inherits. |
| `config/campaigns/` | Campaign definitions, one file per run. |
| `config/requirements/requirements.yaml` | The single machine-readable requirement set. |
| `tests/` | Test suite, mirroring the layer packages. |
| `tests/fixtures/` | Checked-in small data, used in place of any network access. |
| `data/external/` | Downloaded external data.  Not checked in.  See `docs/data-sources.md`. |
| `results/` | Run outputs, one directory per study, each with a manifest. |
| `docs/` | Architecture, work packages, traceability, data, bridge, verification, ADRs. |
| `scripts/` | One-off maintenance utilities, for example `import_srm.py`. |

## Claims

The dissertation advances five claims.  The simulator exists to close them.

| Claim | Statement | Status |
| --- | --- | --- |
| C1 | Sub-watt FMCW from electrically small antennas yields usable topside ionograms across the band. | Resolved affirmative where a measured noise reference exists.  Open below 0.6 MHz, which is why WP12 was created. |
| C2 | Range, density precision and along-track resolution are simultaneously satisfiable. | Open and quantified. |
| C3 | Spacecraft motion separates Doppler from range using one transmit and receive pair. | Open.  Principal contribution.  Unstarted. |
| C4 | Leakage is suppressible with commercial software-defined radio hardware. | Open.  This is now the claim the concept turns on. |
| C5 | A vector receive configuration delivers arrival angle and polarization mode separation at this scale. | Open, goal level. |

C1 is stated as resolved only where a measured noise reference exists.  Below
about 0.6 MHz no such reference has been found in the literature, so the claim
is open there and the honest statement is that the band edge derived at
0.33 MHz rests on extrapolated noise.  WP12 exists to remove that extrapolation.

## Work package status

| WP | Title | Status |
| --- | --- | --- |
| WP0 | Close the science requirements | Complete (2026-10-20).  Requirement set at Rev G as amended 2026-10-20 (72 requirements). |
| WP1 | Ionosphere truth model and perturbation library | Partial.  IRI-2020 in use; NeQuick 2 and the depletion generator outstanding. |
| WP2 | Propagation and echo model | Partial.  Absorption delivered and tested; ray-traced tier and off-vertical echoes outstanding. |
| WP3 | Instrument model and link budget | In progress.  Link budget complete; Phase 0 (reference echo, INS-10 trace) done 2026-10-20.  Leakage and de-chirped residual next.  **Critical path.** |
| WP4 | Baseline signal processing chain | Not started.  Next package to start. |
| WP5 | Novel algorithms | Not started.  Motion-based Doppler and range decoupling, leakage suppression. |
| WP6 | Retrieval and error budget | Not started. |
| WP7 | Trade study and design closure | Not started.  Feasibility maps per configuration. |
| WP8 | Verification and validation | Partial.  Link budget suite and archive comparison delivered. |
| WP9 | Prototype specification | Not started.  Bench build. |
| WP10 | Plasma bubble modelling | Not started, goal level. |
| WP11 | Vector receive and arrival angle | Not started, goal level. |
| WP12 | Sub-megahertz sky noise reference | Not started, newly created. |

`docs/work-packages.md` carries the full table with claims served, implementing
modules and the next action for each, plus the sequencing argument.

## Adding a new tier

A tier is a class that satisfies the protocol for its layer in
`src/tss/interfaces.py` and registers itself.  The steps are:

1. Create a module under the layer package, for example
   `src/tss/l1_propagation/raytrace.py`.
2. Write the docstring metadata block required by `CONTRIBUTING.md`, naming the
   layer, work package, claims, status and requirements the module touches.
3. Decorate the class:

   ```python
   from tss.registry import register

   @register("l1", "raytrace")
   class RayTracedPropagation:
       def __init__(self, *, max_hops: int = 1, step_km: float = 1.0) -> None:
           ...

       def describe(self) -> dict:
           return {"name": self.name, "max_hops": self.max_hops,
                   "step_km": self.step_km, "engine": "PHaRLAP 4.7.4"}

       def echoes(self, state, platform, freq_hz, rng):
           ...
   ```

4. Import the module from the layer package `__init__.py`, which is what makes
   registration happen when `tss` is imported.  Keep any heavy or optional
   dependency, PHaRLAP and the ionosphere model bindings included, inside the
   method that needs it rather than at module scope.
5. Add tests under `tests/<layer package>/` covering at least one closed-form
   limiting case.
6. Select it from a campaign:

   ```yaml
   tiers:
     l1: {name: raytrace, max_hops: 1, step_km: 0.5}
   ```

`tss tiers` will then list it, and `tss config <campaign>` will validate that
the campaign can be built before anything is run.

## Results, manifests and reproducibility

Every run writes its outputs under `results/<study or campaign>/` together with
a `manifest.json`.  The manifest records the git commit and whether the tree was
dirty, the fully resolved configuration and its SHA-256, the tier
implementations selected with their settings, the RNG seed, the package
versions and host, and every external input file with its checksum.

A result whose manifest says the tree was dirty is a result that cannot be
reproduced.  That is stated in the manifest rather than hidden, and
`--strict` refuses to produce one.

Figures are regenerated from stored results by `tss figure <name>`, using the
same driver that produces the corresponding table.  Nothing in `results/` is
checked in; the manifest plus the configuration is what is checked in, and the
result is rebuilt from those.

## Licence

Code in this repository is released under the MIT licence; see `LICENSE`.

External data sets are not covered by that licence and carry their own terms.
The NASA SPDF Alouette and ISIS holdings, the wspr.live spot data and the
data.gov absorption measurements each have attribution requirements set out in
`docs/data-sources.md`.  IRI-2020, NeQuick 2 and PHaRLAP are obtained and
licensed separately by the user and are not redistributed here.

## Citation

If you use this code, cite the dissertation.  A machine-readable record is in
`CITATION.cff`.

```
Moss, A. B.  A low size, weight and power FMCW topside ionospheric sounder for
CubeSat constellations.  PhD dissertation, Utah State University.  In
preparation.  Software: tss-sim.
```
