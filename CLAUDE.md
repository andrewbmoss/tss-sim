# tss-sim: the dissertation's research code

Modelling and simulation for the FMCW topside ionosphere sounder (TSS), the
core of Andrew B. Moss's PhD dissertation (ECE, Utah State University).  The
claims C1-C5 are closed here, and the numbers in the Technical Volume and the
dissertation come from studies run here.

The written side is `C:\Users\Andrew\PhD\Dissertation` (its own repo, open to
every session through `.claude/settings.json`).  Read its `CLAUDE.md` before
touching documents there.

## Read first, every session

1. `RESUME.md`: what the code does today, the phase log, what comes next, and
   the decisions waiting on Andrew.  It is the source of truth for status.
2. `CONTRIBUTING.md`: the rules below in full.
3. `README.md` and `docs/architecture.md` for the L0-L5 layout;
   `docs/work-packages.md` for WP status.

## Environment

- `conda activate tss` before anything.  Plain `python` in PowerShell is a bare
  3.13 without the dependencies.  Ask before installing or upgrading packages
  (environment changes go in `environment.yml`).
- MATLAB and PHaRLAP back the default `iri` tier (`model: pharlap`); `smoke`
  and the contract suite need neither.  MATLAB starts from
  `matlab/startup_tss.m`.  Bridge rules: `docs/matlab-python-bridge.md`.

## Commands

```
pytest -q -m "not slow"          # fast suite, under about a minute
pytest -q                        # everything available on this machine
ruff check src tests scripts     # must stay clean
tss verify                       # what verification has and has not run
tss run <campaign> --out results/<name> --strict
tss study <name>                 # studies that produce dissertation tables
tss figure <name>
```

## Rules that protect the dissertation

1. **Results for the dissertation are run with `--strict`** from a clean,
   committed tree, so the manifest records a reproducible commit.  Exploratory
   runs without it are fine but must never be cited.
2. **A commit that changes a number cited in the dissertation says so** in its
   body and names the TV table or figure (e.g. `c3r1budget`,
   `tab:fieldinversion`).  Then update the Technical Volume in
   `..\Dissertation\TechnicalVolume\` in the same piece of work, or list the
   mismatch for Andrew.
3. **Reference results in `data/reference/` are checked in and pinned.**  Do not
   regenerate or overwrite one without being asked; a rerun goes to a new
   folder or is compared against the pinned one.
4. **Requirements:** `config/requirements/srm_catalog.yaml` is generated from
   `..\Dissertation\TSS_SRM_RevG.xlsx` by `scripts/import_srm.py`; never edit it
   by hand.  Threshold changes follow `docs/requirements-traceability.md` so the
   YAML, the Rev G document and the SRM stay in step.
5. **Code conventions** (CONTRIBUTING.md): SI with unit suffixes (`alt_km`,
   never `alt`); one tier per module, registered by decorator; no heavy
   imports at module scope; the docstring metadata block on every module;
   tests mirror the layers; no test touches the network; never loosen a
   tolerance to fix a Python/MATLAB discrepancy.
6. **Commits:** one logical change, `scope: imperative subject` under 72
   characters, body ending `Refs: WPx, Cx, REQ-IDs`.  Docs move with the code
   (`docs/architecture.md`, `docs/work-packages.md` and the README table
   together, ADRs for lasting decisions).
7. **Never push**, tag, rewrite history or delete branches without asking.
   Local commits are fine.
8. **Ask before long runs.**  Many studies take minutes to hours; state the
   expected run time and run a smoke or reduced case first.
9. **Keep RESUME.md current**: update it at the end of each phase or piece of
   work, in the existing style.

## Decisions and project notes

Plans and decisions recorded in the claude.ai Project are exported to
`..\Dissertation\Plans\project-notes\` under their original file names, one
subfolder per work package (`README.md` there is the index).  A note cited here
as `claude/<name>.md` keeps that file name but may be in any subfolder
(`WP5\C3-research\`, `requirements\`, ...): find it by file name anywhere under
`..\Dissertation\Plans\project-notes\`.  If it is missing, say so and ask
Andrew to export it.  Andrew decides at the gates (G1, G2, ...); propose, do not
decide.  Copies of the WP plans are in `..\Dissertation\Plans\`.

## Using helpers (to save plan usage)

- Long test output, logs, docs or papers: the `digest` subagent.
- Running pytest, ruff or a smoke campaign: the `test-runner` subagent.
- Finding where a tier, requirement or study lives: `Explore`.
