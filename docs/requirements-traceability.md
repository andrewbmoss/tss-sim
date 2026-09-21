# Requirements traceability

`config/requirements/requirements.yaml` is the single machine-readable source of
the requirement set.  L5 scores against that file and against nothing else.  No
threshold is written into a tier, a study or a figure script.

The reason is the first governing principle of Section 3.1 of the proposal: no
simulated performance figure is scored pass or fail until the number it is
scored against exists and is traced to a source.  Putting the numbers in one
file makes that principle checkable by a program rather than by a reader.

The loader and the scorer are `src/tss/requirements.py`.  The command line view
is `tss requirements`, which prints the set and validates it, and exits with
status 2 if any entry cannot be used to score.

## Schema of an entry

```yaml
revision: E

requirements:

  INS-01:
    title: Sounding band and profile coverage
    sense: range
    unit: MHz
    threshold: [0.33, 20.0]
    goal: [0.30, 20.0]
    claims: [C1]
    work_packages: [WP0, WP1, WP12]
    source: TSS_Science_Requirements_RevE.docx section 4.2
    derivation: ""
    notes: >
      95% profile coverage across the band.  The lower edge is derived from
      the minimum sounding frequency study over 1.2 million IRI-2020 profiles.
      The external noise reference below 0.6 MHz is extrapolated, so a verdict
      on coverage below that frequency carries the extrapolated-noise flag.

  INS-10:
    title: Sensitivity, signal to noise per range cell
    sense: min
    unit: dB
    threshold: 10.0
    goal: 13.0
    claims: [C1]
    work_packages: [WP3, WP8]
    source: derived
    derivation: docs/derivations/ins10_snr_threshold.md
    notes: >
      Set by the trace detection probability needed to meet the INS-01
      coverage figure at the worst-case absorption condition.

  INS-18:
    title: Vector receive capability
    sense: bool
    unit: ""
    threshold: true
    goal: null
    claims: [C5]
    work_packages: [WP11]
    source: TSS_Science_Requirements_RevE.docx section 4.6
    derivation: ""
    notes: Goal level.  Not required for the baseline configuration.
```

The fields:

| Field | Meaning |
| --- | --- |
| `id` | The mapping key, for example `INS-01` or `SCI-09`.  `INS-` for instrument requirements, `SCI-` for science requirements.  IDs are never reused.  A retired requirement keeps its ID and gains a `notes` entry saying it was retired and in which revision. |
| `title` | One short line.  This is what appears in the assessment table, so keep it under about 60 characters. |
| `sense` | One of `min`, `max`, `range`, `bool`.  See below. |
| `unit` | The unit the threshold and goal are expressed in, as a human-readable string.  Empty for dimensionless and for `bool`. |
| `threshold` | The number that must be met.  A two-element list for `sense: range`.  `true` for `bool`.  Required for every sense except `bool`. |
| `goal` | The number one would like to meet.  May be `null`, in which case there is no goal verdict and the best achievable outcome is `threshold`. |
| `claims` | The claims this requirement contributes to, from C1 to C5.  May be empty for a requirement that is a constraint rather than evidence, though in practice most map to at least one. |
| `work_packages` | The work packages that produce an achieved value for this requirement, or that validate it. |
| `source` | Where the number comes from.  Either a citation into the Rev E document, an external standard, a literature reference, or the literal string `derived`. |
| `derivation` | Required when `source: derived`.  A path to the document or notebook that derives the number.  Ignored otherwise. |
| `notes` | Free text.  Caveats, conditions, and anything a reader needs in order not to misread the verdict.  This is where conditional results are recorded. |

## The four senses

`sense` says which direction is better, and it is what `Requirement.evaluate`
dispatches on.

**`min`** means the achieved value must be at or above the stated number.  Used
for sensitivity, isolation, coverage and coherent processing interval.  INS-10
at 10 dB SNR per range cell, INS-15 at 1.0 s coherent processing and INS-23 at
40 dB transmit-receive isolation are all `min`.  The verdict is `goal` if the
achieved value reaches the goal, `threshold` if it reaches the threshold, and
`fail` otherwise.  The margin is achieved minus the level that was met, or
achieved minus threshold when it failed, so a negative margin says how far
short it fell.

**`max`** means the achieved value must be at or below the stated number.  Used
for limits on error, resolution, mass and power.  INS-05 at 10 km virtual range
resolution, SCI-09 at 50 km along-track sampling, SCI-23 at 50 m/s zonal drift
resolution and INS-26 at 2e-9 Allan deviation are `max`.

**`range`** means the achieved value must lie inside `[lower, upper]`.  INS-01's
sounding band is the case this exists for.  `threshold` and `goal` are both
two-element lists, and `validate()` rejects an entry whose `range` threshold is
not a pair.  The margin is the distance to the nearer edge.

**`bool`** is a capability that is either demonstrated or not.  INS-18 vector
receive is the example.  There is no margin, and the verdict is `goal` when the
capability is demonstrated and `fail` when it is not.  There is no intermediate
`threshold` verdict for a boolean, which is a slight asymmetry with the other
senses and is accepted because a half-demonstrated capability is not a useful
thing to report.

A requirement that is not present in the achieved mapping L5 hands to
`RequirementSet.score` is reported as `not_evaluated`, not as a failure.  That
distinction matters.  A campaign that did not exercise Doppler should not be
recorded as having failed INS-09.

## A requirement without a source cannot be used to score

`RequirementSet.validate()` returns a problem for any entry with an empty
`source`, and `tss requirements` exits with status 2 when there are problems.
It also returns a problem for an entry whose `source` is `derived` but which
names no `derivation`.

The rule is not a formality.  The difference between a number that came from the
mission science case and a number that someone picked because it looked
reasonable is the difference between a defensible result and an unfalsifiable
one.  If a threshold is invented in order to make a result pass, that fact is
visible in the diff of this file rather than buried in a script.

The practical consequence when adding a requirement is that you cannot add it
and fill in the source later.  Either cite the Rev E document, or write the
derivation note first and cite that.  A derivation note is short: what the
number is, what it is derived from, the arithmetic, and what would change it.
They live under `docs/derivations/`.

The same rule is why `threshold` is required for every non-boolean sense.  An
entry with a goal and no threshold would silently score everything below the
goal as a failure.

## Keeping the YAML in step with the documents

Three artefacts have to agree:

1. `TSS_Science_Requirements_RevE.docx`, the prose requirement document.  It is
   the authority.  It carries the rationale, the negotiation history and the
   sign-off.
2. `TSS_SRM_RevE.xlsx`, the science requirements matrix.  It carries the
   traceability from science objective to measurement to instrument
   requirement, in a form that can be tabulated.
3. `config/requirements/requirements.yaml`, which is what actually scores
   results.

The YAML is downstream of the other two.  It is never edited first.  The
procedure for a change is:

1. Change the Rev E document.  If the change alters a threshold, a goal or a
   sense, the revision letter goes up; Rev E is not edited in place once
   results have been scored against it.  Add the new revision letter to the
   `revision` key in the YAML at the same time.
2. Update the SRM spreadsheet so the traceability chain still closes.
3. Run `python scripts/import_srm.py TSS_SRM_RevE.xlsx` to regenerate a
   candidate YAML from the spreadsheet, and review the diff against the
   committed file.  The script exists to make this a diff review rather than a
   retype.  It does not commit anything; it writes to a scratch path and prints
   the diff.
4. Carry across by hand anything the spreadsheet does not hold: `derivation`
   paths, `notes`, and the `claims` mapping where it is richer in the YAML than
   in the SRM.
5. Run `tss requirements`.  A non-zero exit means the set is not scorable.
6. Run `pytest -q tests/test_metadata.py`, which checks that every requirement
   ID named in a module docstring metadata block exists in the set.  This is
   what catches a renumbering that was applied to the document but not to the
   code.
7. Commit the document change, the spreadsheet change and the YAML change
   together, with a `Refs:` line naming the requirement IDs.

Results already scored against an earlier revision are not silently
invalidated.  Each run manifest records the requirement set revision it scored
against, so an old result remains interpretable; it just says Rev D on it.  When
a threshold changes, the affected campaigns are rerun and the old results are
kept until the new ones exist, rather than deleted first.

## Reading the traceability in the other direction

`RequirementSet` supports both directions, which is what the dissertation
chapters need.

```python
from tss.requirements import load_requirements

reqs = load_requirements()
reqs.for_claim("C3")          # every requirement C3 is scored against
reqs.for_work_package("WP3")  # every requirement WP3 has to produce
reqs["INS-23"].evaluate(42.1) # score one achieved value
```

From the command line:

```bash
tss requirements --claim C4
tss requirements --json > results/requirements_revE.json
```

The module docstring metadata block described in `CONTRIBUTING.md` closes the
last link, from a requirement ID back to the module that produces an achieved
value for it.  With the YAML, the metadata blocks and the manifests, the chain
runs: science objective, SRM row, Rev E requirement, YAML entry, implementing
module, campaign run, manifest, scored assessment.  Any break in that chain
should be visible to one of the checks above rather than to a reader.
