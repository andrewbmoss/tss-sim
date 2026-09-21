# ADR 0004  Requirements as machine-readable data

Status: accepted
Date: 2025
Supersedes: nothing

## Context

The first governing principle of Section 3.1 of the proposal is that no
simulated performance figure is scored pass or fail until the number it is
scored against exists and is traced to a source.  Stated as a principle it is
uncontroversial.  Stated as a practice it is difficult, because the natural way
to write a scoring script is to put the threshold in the script.

Thresholds in scripts fail in three ways.  They drift apart when the same
requirement is scored in two places.  They get quietly adjusted to make a
borderline result pass, and the adjustment is invisible in a diff full of other
changes.  And they lose their provenance, so a number that came from the mission
science case becomes indistinguishable from a number someone picked because it
looked reasonable.

The requirement set is also not small and not static.  It runs to dozens of
entries across instrument and science requirements, it is at Rev E, and it
exists in two upstream documents that have to stay consistent with each other:
`TSS_Science_Requirements_RevE.docx` and `TSS_SRM_RevE.xlsx`.

The requirements also have to be readable in two directions.  From a
requirement, which work package produces an achieved value for it.  From a
claim, which requirements it is scored against.  The dissertation chapters need
both.

## Decision

**One machine-readable file**, `config/requirements/requirements.yaml`, is the
single source.  L5 scores against it and against nothing else.  No threshold
appears in a tier, a study or a figure script.

**Each entry carries** `id`, `title`, `sense`, `unit`, `threshold`, `goal`,
`claims`, `work_packages`, `source`, `derivation` and `notes`.  The loader and
scorer are `src/tss/requirements.py`.

**Four senses**: `min`, `max`, `range` and `bool`.  `sense` says which direction
is better and is what the scorer dispatches on.  `range` exists for INS-01's
sounding band, where the achieved value has to lie inside an interval.  `bool`
exists for capabilities such as INS-18 vector receive, which are demonstrated or
not.

**A requirement without a source cannot be used to score.**
`RequirementSet.validate()` reports a problem for an empty `source`, and for a
`source: derived` with no `derivation` naming where the number comes from.
`tss requirements` exits with status 2 when there are problems.  A threshold is
required for every non-boolean sense, because an entry with only a goal would
silently score everything below the goal as a failure.

**A requirement not evaluated is reported as `not_evaluated`, not as a
failure**, and requirements absent from the achieved mapping are still shown.
An assessment always carries the full set.

**The YAML is downstream** of the Rev E document and the SRM spreadsheet and is
never edited first.  `scripts/import_srm.py` regenerates a candidate from the
spreadsheet so that keeping them in step is a diff review rather than a retype.
The procedure is in `docs/requirements-traceability.md`.

## Consequences

Changing a threshold is a one-line diff in one file, in a commit of its own,
with the source field visible in the same hunk.  Adjusting a number to make a
result pass is still possible, and it is now conspicuous, which is the most that
can be asked of a process.

Scoring is uniform.  Every requirement goes through the same evaluator and
produces the same record shape, so the assessment table in a results chapter is
generated rather than assembled.  Adding a requirement adds a row everywhere it
should appear, with no code change.

Traceability runs in both directions through `for_claim` and `for_work_package`,
and the docstring metadata blocks required by `CONTRIBUTING.md` close the last
link from a requirement ID back to the module that produces an achieved value
for it.  `tests/test_metadata.py` fails when a module names a requirement ID
that is not in the set, which is what catches a renumbering applied to the
document but not to the code.

The source rule has a cost that is felt when adding a requirement.  It cannot be
added with the intention of filling in the source later, because validation will
refuse it.  Either the Rev E document is cited or a derivation note is written
first.  That is friction at exactly the moment when friction is least welcome,
and it is the point of the rule.

Three artefacts now have to be kept consistent instead of two, and the YAML is
the one a reader of the dissertation will never see.  If it drifts from Rev E,
every scored result is wrong in a way that no test will catch, because the tests
check that the YAML is well formed and not that it matches the document.  The
import script reduces this risk and does not remove it.  It is the largest
unmitigated weakness in this decision and it is worth restating whenever the
revision letter changes.

Results scored against an earlier revision are not invalidated silently, because
the manifest records the requirement set revision.  An old result stays
interpretable; it just says Rev D on it.  The consequence is that a results
chapter mixing revisions has to say so, and the affected campaigns are rerun
before the old results are deleted rather than after.

Reporting `not_evaluated` entries makes an assessment more verbose than it
needs to be, since a campaign that only exercises L0 and L1 still prints the
full instrument requirement table.  That verbosity is accepted because a
silently shortened requirement table reads as a clean sweep, and a clean sweep
is exactly the thing an examiner should be able to check.

A side benefit is that `tss requirements --claim C3` gives a direct answer to
what evidence a claim needs, which turns out to be the most useful single
command for planning work.
