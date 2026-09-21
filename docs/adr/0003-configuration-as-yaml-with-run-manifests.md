# ADR 0003  Configuration as YAML campaigns with run manifests

Status: accepted
Date: 2025
Supersedes: nothing

## Context

The dissertation promises that every figure is regenerable from a stored
configuration by the same driver that produces the corresponding table.  That
promise fails in the ordinary way: settings end up as literals in driver
scripts, a script is edited to make a new figure, the old figure can no longer
be reproduced, and six months later nobody can say which version of the script
made which plot.

The specific things that have to be pinned for a result to be reproducible here
are the instrument configuration, the platform configuration, which tier
implementation was selected at each of the six layers and with what settings,
the RNG seed, the requirement set revision it was scored against, and the
external data files that were read.  That is a lot of state, and it changes
often during a study.

There is also a real need to vary one thing at a time.  A trade study is a
family of runs differing in one parameter, and writing out a full configuration
per member of that family invites copy-paste divergence in the parts that were
supposed to be identical.

## Decision

**A run is defined by a YAML campaign file** under `config/campaigns/`.  The
campaign names the instrument configuration, the platform configurations to
sweep, the tier implementations for each layer, the seed and the requirement set
to score against.  Nothing is hard-coded in a driver script.

**Resolution order**, each overriding the previous: `config/default.yaml`, then
any file named in the campaign's `inherits` list in order, then the campaign
file itself, then command-line overrides of the form `--set a.b.c=value`.
Mappings merge key by key.  Every other type, lists included, replaces
wholesale, because a campaign that lists three configurations means those three
and not those three appended to the defaults.

**A campaign is validated before it is run.**  `tss config <campaign>` resolves
and prints it.  Validation rejects a campaign that selects no tier for a layer,
that names a tier which is not registered, or that names no configurations at
all.  The last one is refused because Section 1.3 of the proposal reports
results per configuration and a result without its configuration is not a
result.

**Every run writes a manifest** next to its outputs, through
`src/tss/provenance.py`.  The manifest records the git commit and whether the
working tree was dirty, the fully resolved configuration and its SHA-256, the
tier implementations with their `describe()` output, the seed, the Python
version and host, and every external input file with its SHA-256 and size.

**A publication run refuses a dirty working tree.**  `check_clean(strict=True)`
raises, and the `--strict` flag on `tss run` and `tss study` sets it.
Exploratory runs omit the flag and get a manifest that honestly records the tree
as dirty.

**Results are not checked in.**  `results/` is ignored.  What is checked in is
the configuration and the code, and the result is rebuilt from those.

## Consequences

A figure can be traced back to the exact settings that produced it, including
which fidelity tier was selected at each layer, which is the piece that would
otherwise be hardest to recover.  The config hash makes it cheap to tell whether
two results came from the same settings without comparing the files by eye.

The inheritance chain keeps a trade study family consistent.  A campaign that
inherits a baseline and overrides one key differs from the baseline in exactly
that key, and the diff shows it.

Command-line overrides make a one-off variation possible without editing a file,
and because the override is applied before hashing, the manifest records the
resolved value rather than the file's value.  A run made with `--set` is
therefore still reproducible, though only from its manifest and not from the
campaign file alone.  That is a genuine sharp edge: a reader who looks at the
campaign file and not at the manifest will see the wrong number.

The dirty-tree rule is the part most likely to be found annoying in practice.
It means a publication run cannot be made from a working tree with an uncommitted
experiment in it, which during active development is most of the time.  The
alternative, producing results that claim a commit they were not built from, is
worse.  The escape hatch is to omit `--strict`, and the cost of using the escape
hatch is a manifest that says the result cannot be reproduced.  That is the
right trade: the fact is recorded rather than prevented.

YAML brings its own hazards.  Unquoted values parse in surprising ways, deep
merge can hide where a value came from, and a typo in a key name is silently
accepted as a new key rather than rejected.  Validation catches the structural
cases, tier names and missing layers, but it does not catch a misspelled tuning
parameter that a tier then never reads.  This is a known weakness and the
mitigation is that `describe()` records what the tier actually used, so the
manifest shows the effective value rather than the requested one.  Comparing
those two is how a silently ignored key is found.

Not checking in results means a figure in the dissertation cannot be
regenerated without rerunning the campaign, which for the large sweeps takes
real time.  That is accepted because the alternative is a repository full of
large binary outputs whose provenance still has to be trusted.

A side benefit to hashing the resolved configuration is that a rerun can be
skipped when the hash already exists in `results/`, which makes a partially
completed campaign resumable without any extra bookkeeping.
