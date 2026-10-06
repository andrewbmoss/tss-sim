"""Layer: support
WP: WP0, WP3, WP5, WP7, WP8, WP12
Claims: C1, C2, C3, C4, C5
Status: partial
Requirements: -

The study registry, the study runner and the campaign runner.

A study is a named, reproducible campaign that produces a result the
dissertation cites.  It is not a script that happens to live in the repository.
The difference is that a study runs from a stored configuration, writes its
outputs and a run manifest into one timestamped directory, and is addressable by
name from the command line:

    tss study list
    tss study link_budget
    tss study min_frequency --set grid.alt_step_m=10000

Section 4.2 of the proposal promises that every figure is regenerable from a
stored configuration by the same driver that produces the corresponding table.
That promise is why the study is the driver: the table and the figure read the
same stored result, and :mod:`tss.figures` loads what a study wrote rather than
recomputing anything.  A number recomputed at plotting time is not the number
that was scored.

WHAT EVERY STUDY WRITES
-----------------------
Into ``results/<study>/<timestamp>/``:

    manifest.json   provenance, via :class:`tss.provenance.Manifest`
    summary.md      the short prose result, with its caveats
    *.csv, *.npz    the tables the dissertation cites

and ``results/<study>/latest.json``, a small pointer record naming the most
recent run directory and its headline numbers.  The pointer exists because the
figure registry declares its inputs as ``results/<study>/*.json`` and a reader
looking for the current answer should not have to sort timestamps by hand.

STATUS IS PART OF THE DESCRIPTOR
--------------------------------
``implemented`` means the study runs and its result is citable.  ``partial``
means it runs but bounds the answer rather than settling it, which is the honest
description of :mod:`tss.studies.sky_noise`.  ``stub`` means it raises
:class:`NotImplementedError` naming the work package that will fill it in.  Two
studies are stubs and both are named in the risk section of the proposal, so
they are registered rather than omitted: a reader who asks for
``doppler_separation`` should be told that WP5 owns it, not that no such study
exists.

EXIT STATUS
-----------
:func:`run_study` and :func:`run_campaign` return a process exit status rather
than raising, because :mod:`tss.cli` hands it straight to the shell:

    0   the study ran and wrote its outputs
    1   a campaign scored ``fail`` on at least one threshold requirement
    2   the study is a stub; the message names the work package
    3   the study could not run for want of an input, such as an uncached
        archive; the message names the file and the command that fetches it
    4   no study of that name is registered
    5   the working tree is dirty and ``--strict`` was given

Those are distinguished rather than collapsed because they call for different
actions, and a nightly job that reruns the study set should be able to tell a
missing archive from a failed requirement.
"""

from __future__ import annotations

import json
import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TextIO

__all__ = [
    "Study",
    "StudyResult",
    "STUDIES",
    "register_study",
    "available_studies",
    "run_study",
    "run_campaign",
    "results_root",
    "StudyInputMissing",
]


class StudyInputMissing(RuntimeError):
    """A study cannot run because an input it needs is not on this machine.

    Raised rather than returned so that a study body can give up at the point it
    discovers the problem.  :func:`run_study` turns it into exit status 3 and
    prints the message, which must name the missing file and the command that
    produces it.
    """


@dataclass
class StudyResult:
    """What one study produced.

    Attributes
    ----------
    outputs
        Files written, in the order they were written.  Copied into the
        manifest.
    summary
        Markdown body for ``summary.md``.  The runner adds the heading and the
        provenance footer, so a study writes the result and its caveats and
        nothing else.
    metrics
        Headline numbers, JSON-serialisable.  These go into ``latest.json`` and
        are what a later figure or a citation check reads.
    exit_code
        Non-zero when the study ran but its result should be treated as a
        failure.  A study that simply could not run raises instead.
    """

    outputs: list[Path] = field(default_factory=list)
    summary: str = ""
    metrics: dict[str, Any] = field(default_factory=dict)
    exit_code: int = 0


#: What a study callable is handed and what it must return.  The manifest is
#: passed in rather than returned so that a study can record the external files
#: it consumed at the moment it opens them.
StudyCallable = Callable[..., StudyResult]


@dataclass(frozen=True)
class Study:
    """One named study: what it answers, what it serves, and whether it works."""

    name: str
    description: str
    work_packages: tuple[str, ...]
    claims: tuple[str, ...]
    status: str                       # "implemented" | "partial" | "stub"
    run: StudyCallable
    requirements: tuple[str, ...] = ()
    blocked_on: str = ""              # the work package a stub waits on

    def as_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "work_packages": list(self.work_packages),
            "claims": list(self.claims),
            "status": self.status,
            "requirements": list(self.requirements),
            "blocked_on": self.blocked_on,
        }


STUDIES: dict[str, Study] = {}

_STATUSES = ("implemented", "partial", "stub")


def register_study(
    name: str,
    *,
    description: str,
    work_packages: Sequence[str],
    claims: Sequence[str] = (),
    status: str = "stub",
    requirements: Sequence[str] = (),
    blocked_on: str = "",
) -> Callable[[StudyCallable], StudyCallable]:
    """Register a study.  Used as a decorator on the study's ``run`` function.

    ``work_packages`` and ``status`` are required rather than optional.  A study
    whose work package is not recorded cannot be placed against the plan in
    ``docs/work-packages.md``, and a study whose status is not recorded invites a
    reader to take a stub for a result.
    """
    if status not in _STATUSES:
        raise ValueError(f"status must be one of {_STATUSES}, got {status!r}")
    if not work_packages:
        raise ValueError(f"study {name!r} must name at least one work package")

    def _decorate(fn: StudyCallable) -> StudyCallable:
        if name in STUDIES:
            raise ValueError(f"study {name!r} is already registered")
        STUDIES[name] = Study(
            name=name,
            description=description,
            work_packages=tuple(work_packages),
            claims=tuple(claims),
            status=status,
            run=fn,
            requirements=tuple(requirements),
            blocked_on=blocked_on,
        )
        return fn

    return _decorate


def available_studies() -> dict[str, dict[str, Any]]:
    """Every registered study, for ``tss study list`` and for tests."""
    return {name: study.as_dict() for name, study in sorted(STUDIES.items())}


def results_root() -> Path:
    """``results/`` under the repository root."""
    from tss.config import repo_root

    return repo_root() / "results"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _base_config(overrides: Sequence[str]) -> dict[str, Any]:
    """``config/default.yaml`` with the command-line overrides applied.

    A study is not a campaign: it does not name a tier for every layer and it is
    not validated against the campaign schema.  It still resolves its numbers
    from the stored defaults rather than from literals in its own source, so that
    ``--set`` reaches a study the same way it reaches a campaign and the resolved
    mapping can be hashed into the manifest.
    """
    from tss.config import apply_overrides, config_dir, load_yaml

    path = config_dir() / "default.yaml"
    config = load_yaml(path) if path.exists() else {}
    return apply_overrides(config, overrides)


def _print_study_table(stream: TextIO) -> None:
    for name, record in available_studies().items():
        marker = {"implemented": "  ", "partial": "~ ", "stub": "* "}[record["status"]]
        wps = ",".join(record["work_packages"])
        claims = ",".join(record["claims"]) or "-"
        print(f"{marker}{name:20s} {wps:16s} {claims:10s} {record['description']}",
              file=stream)
    print(
        "\n~ bounds the answer rather than settling it.  "
        "* is a stub; the message names the work package that owns it.",
        file=stream,
    )


def run_study(
    name: str,
    out: str | Path | None = None,
    overrides: Sequence[str] = (),
    strict: bool = False,
    *,
    stream: TextIO | None = None,
) -> int:
    """Run one named study and return a process exit status.

    ``name`` of ``"list"`` prints the registry and returns 0, which is how a
    reader finds out what exists without reading this package.
    """
    out_stream = stream if stream is not None else sys.stdout

    if name == "list":
        _print_study_table(out_stream)
        return 0

    study = STUDIES.get(name)
    if study is None:
        known = ", ".join(sorted(STUDIES)) or "nothing"
        print(
            f"no study named {name!r}; registered: {known}.  "
            "Run `tss study list` for what each one answers.",
            file=out_stream,
        )
        return 4

    from tss.provenance import Manifest, check_clean

    try:
        check_clean(strict=strict)
    except RuntimeError as exc:
        print(str(exc), file=out_stream)
        return 5

    config = _base_config(overrides)
    out_dir = Path(out) if out is not None else results_root() / name / _timestamp()
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = Manifest(
        study=name,
        config=config,
        notes=study.description,
        tiers={
            "status": study.status,
            "work_packages": list(study.work_packages),
            "claims": list(study.claims),
        },
        seed=int(config.get("seed", 0)),
    )

    try:
        result = study.run(out_dir, config=config, manifest=manifest)
    except NotImplementedError as exc:
        _remove_if_empty(out_dir)
        print(f"study {name!r} is a stub.\n  {exc}", file=out_stream)
        return 2
    except StudyInputMissing as exc:
        _remove_if_empty(out_dir)
        print(f"study {name!r} cannot run here.\n  {exc}", file=out_stream)
        return 3

    summary_path = _write_summary(out_dir, study, result)
    result.outputs.append(summary_path)
    manifest.outputs = [str(p) for p in result.outputs]
    manifest_path = manifest.write(out_dir)
    _write_pointer(name, out_dir, study, result, manifest_path)

    print(
        f"{study.name} [{study.status}] wrote {len(result.outputs)} file(s) to {out_dir}",
        file=out_stream,
    )
    for path in result.outputs:
        print(f"  {path}", file=out_stream)
    return result.exit_code


def _remove_if_empty(directory: Path) -> None:
    """Take back the directory a study never wrote into.

    A stub that left an empty timestamped directory behind would make the
    results tree look like it holds a result.
    """
    try:
        if directory.is_dir() and not any(directory.iterdir()):
            directory.rmdir()
    except OSError:
        pass


def _write_summary(out_dir: Path, study: Study, result: StudyResult) -> Path:
    """Write ``summary.md``: the heading, the study's own prose, the footer."""
    created = datetime.now(timezone.utc).isoformat(timespec="seconds")
    lines = [
        f"# {study.name}",
        "",
        study.description,
        "",
        f"Status: **{study.status}**.  "
        f"Work packages: {', '.join(study.work_packages)}.  "
        f"Claims: {', '.join(study.claims) or 'none stated'}.  "
        f"Requirements: {', '.join(study.requirements) or 'none stated'}.",
        "",
        result.summary.strip(),
        "",
        "---",
        "",
        f"Run {created}.  Provenance is in `manifest.json` next to this file.  "
        f"Rerunning `tss study {study.name}` from the same commit and the same "
        "configuration reproduces it.",
        "",
    ]
    path = out_dir / "summary.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def _write_pointer(
    name: str,
    out_dir: Path,
    study: Study,
    result: StudyResult,
    manifest_path: Path,
) -> Path:
    """Write ``results/<study>/latest.json``.

    Only written when the run directory sits under ``results/<study>/``, which is
    the default.  A run directed elsewhere with ``--out`` is somebody's scratch
    run and must not move the pointer the figures read.
    """
    payload = {
        "study": name,
        "status": study.status,
        "run_directory": str(out_dir),
        "manifest": str(manifest_path),
        "outputs": [str(p) for p in result.outputs],
        "metrics": result.metrics,
        "work_packages": list(study.work_packages),
        "claims": list(study.claims),
    }
    default_parent = results_root() / name
    pointer = (
        default_parent / "latest.json"
        if _same_directory(out_dir.parent, default_parent)
        else out_dir / "latest.json"
    )
    pointer.parent.mkdir(parents=True, exist_ok=True)
    pointer.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    return pointer


def _same_directory(left: Path, right: Path) -> bool:
    """True when two paths name the same directory, without requiring both exist."""
    try:
        return left.resolve() == right.resolve()
    except OSError:
        return str(left) == str(right)


# ---------------------------------------------------------------------------
# Campaigns
# ---------------------------------------------------------------------------


def _sweep_from_config(
    config: Mapping[str, Any], stream: TextIO
) -> tuple[Any, dict[str, Any]]:
    """Build the :class:`~tss.interfaces.SweepSchedule` a campaign asks for.

    The default block is the A1 stepped sweep, whose duration follows from its
    step law (159 steps of 26.8 ms, 4.26 s) and whose period of zero means that
    duration.  A stated block can still be internally inconsistent: a uniform
    sweep of 0.33 to 20 MHz at 5.6 MHz/s takes 3.5125 s, and the earlier default
    stated a period of 3.3 s against it, which is the INS-07 goal rather than a
    sweep duration.
    :func:`tss.l3_signal.base.sweep_duration_s` refuses such a schedule, because
    the alternative is a silently truncated sweep and a silently wrong range
    scale.

    Rather than refuse the run, the driver lengthens the period to the duration
    the band and the rate actually require, prints what it did, and records both
    numbers in the manifest.  The period sets how often a sweep repeats and not
    the range scale, so lengthening it changes the ionogram cadence and nothing
    else, and the cadence is then scored against INS-07 as the number the run
    actually achieved.  A campaign that wants the stated cadence has to state a
    faster ``k_f_hz_per_s`` or a narrower band, which is a design decision and
    not something a driver should make silently.
    """
    from tss.interfaces import SweepSchedule

    from tss import sweep_law

    block = dict(config.get("sweep", {}))
    f_start = float(block.get("f_start_hz", 3.3e5))
    f_stop = float(block.get("f_stop_hz", 2.0e7))
    period = float(block.get("sweep_period_s", 0.0) or 0.0)
    stepped = float(block.get("step_dwell_s", 0.0) or 0.0) > 0.0
    if stepped:
        # The A1 stepped law of Requirements Rev G.  Its rate is the peak step
        # rate and its duration is n_steps * step_dwell_s; see tss.sweep_law.
        from tss.l2_instrument.waveform import schedule_from_config

        law = schedule_from_config({**block, "f_start_hz": f_start, "f_stop_hz": f_stop,
                                    "sweep_period_s": 0.0})
        k_f = float(law.k_f_hz_per_s)
        required = sweep_law.active_duration_s(law)
    else:
        k_f = float(block.get("k_f_hz_per_s", 5.6e6))
        required = abs(f_stop - f_start) / abs(k_f)

    notes: dict[str, Any] = {
        "sweep_law": "stepped" if stepped else "uniform",
        "sweep_duration_s": required,
        "sweep_period_requested_s": period,
        "sweep_period_used_s": period,
        "sweep_period_extended": False,
    }
    if period > 0.0 and required > period * (1.0 + 1e-9):
        print(
            f"  sweep.sweep_period_s is {period:g} s but the band needs "
            f"{required:.4g} s at {k_f:g} Hz/s; the driver is using "
            f"{required:.4g} s and recording both.",
            file=stream,
        )
        period = required
        notes["sweep_period_used_s"] = period
        notes["sweep_period_extended"] = True

    if period <= 0.0:
        period = required
        notes["sweep_period_used_s"] = period
    schedule = SweepSchedule(
        f_start_hz=f_start,
        f_stop_hz=f_stop,
        k_f_hz_per_s=k_f,
        sweep_period_s=period,
        direction=int(block.get("direction", 1)),
        nonlinearity_ppm=float(block.get("nonlinearity_ppm", 0.0)),
        step_dwell_s=float(block.get("step_dwell_s", 0.0) or 0.0) if stepped else 0.0,
        step_fraction=float(block.get("step_fraction", 0.05)) if stepped else 0.0,
        step_cap_hz=float(block.get("step_cap_hz", 150.0e3)) if stepped else 0.0,
    )
    if stepped:
        notes["n_steps"] = sweep_law.step_count(schedule)
        notes["peak_rate_hz_per_s"] = k_f
    return schedule, notes


def _grids(config: Mapping[str, Any]) -> tuple[Any, Any]:
    """The altitude and sounding-frequency grids, from the ``grid`` block.

    ``n_freq`` is read from the same block and defaults to 128.  It is the one
    quantity the campaign schema does not already carry, and a coarse grid is
    what makes the smoke campaign fast enough to run inside a test.
    """
    import numpy as np

    grid = dict(config.get("grid", {}))
    alt = np.arange(
        float(grid.get("alt_min_m", 9.0e4)),
        float(grid.get("alt_max_m", 2.0e6)) + 1.0,
        float(grid.get("alt_step_m", 5.0e3)),
    )
    sweep = dict(config.get("sweep", {}))
    if float(sweep.get("step_dwell_s", 0.0) or 0.0) > 0.0:
        # A stepped sweep sounds one frequency per step, at the step centre,
        # so the echo grid is the step centres and n_freq does not apply.
        from tss import sweep_law
        from tss.l2_instrument.waveform import schedule_from_config

        law = schedule_from_config({
            **sweep,
            "f_start_hz": float(sweep.get("f_start_hz", 3.3e5)),
            "f_stop_hz": float(sweep.get("f_stop_hz", 2.0e7)),
            "sweep_period_s": 0.0,
        })
        return alt, sweep_law.step_centres_hz(law)
    freq = np.linspace(
        float(sweep.get("f_start_hz", 3.3e5)),
        float(sweep.get("f_stop_hz", 2.0e7)),
        max(int(grid.get("n_freq", 128)), 2),
    )
    return alt, freq


def _platform_for(configuration: Mapping[str, Any], epoch: str) -> Any:
    """Platform state for one configuration, velocity included.

    The velocity comes from :func:`tss.platform.orbits.platform_state_at`, which
    computes it analytically from the orbit rather than by differencing two
    positions.  C3 rests on the platform motion being large and precisely known,
    so a campaign that handed L1 a state with no velocity would be quietly
    removing the thing the principal contribution is about.
    """
    from tss.platform.orbits import CircularOrbit, platform_state_at

    orbit_block = dict(configuration.get("orbit", {}))
    orbit = CircularOrbit(
        altitude=float(orbit_block.get("altitude_m", 5.0e5)),
        inclination_deg=float(orbit_block.get("inclination_deg", 97.4)),
        name=str(configuration.get("name", "")),
    )
    return platform_state_at(orbit, 0.0, epoch=epoch)


def run_campaign(
    campaign: str,
    out: str | Path | None = None,
    overrides: Sequence[str] = (),
    strict: bool = False,
    *,
    stream: TextIO | None = None,
) -> int:
    """Run every configuration of a campaign through L0 to L5.

    Loads the campaign with :func:`tss.config.load_campaign`, builds the tiers it
    names with :func:`tss.pipeline.build_tiers`, runs each named configuration
    through a :class:`tss.pipeline.Pipeline`, renders each
    :class:`~tss.interfaces.Assessment` with :mod:`tss.l5_assessment.report`, and
    writes one combined manifest for the campaign alongside the per-configuration
    outputs.

    Returns non-zero if any configuration scored ``fail`` on a threshold
    requirement.  A requirement that was not evaluated is not a failure and is
    not counted as one; it is reported in the per-configuration table, which is
    where ``docs/architecture.md`` says an unevaluated requirement belongs.

    ``strict`` refuses a dirty working tree.  A campaign whose own configuration
    sets ``run.strict_provenance`` is treated as though ``--strict`` had been
    given, because that flag exists to mark a campaign whose output is destined
    for the dissertation.
    """
    import numpy as np

    from tss.config import load_campaign, load_yaml, resolve_path
    from tss.l5_assessment.report import to_dict, write_report
    from tss.pipeline import Pipeline, build_tiers
    from tss.provenance import Manifest, check_clean

    out_stream = stream if stream is not None else sys.stdout

    loaded = load_campaign(campaign, overrides)
    config = loaded.config
    publication = bool(config.get("run", {}).get("strict_provenance", False))
    try:
        check_clean(strict=strict or publication)
    except RuntimeError as exc:
        print(str(exc), file=out_stream)
        return 5

    out_dir = Path(out) if out is not None else results_root() / loaded.name / _timestamp()
    out_dir.mkdir(parents=True, exist_ok=True)

    pipeline = Pipeline(
        build_tiers(loaded.tiers),
        seed=loaded.seed,
        keep_intermediates=bool(config.get("run", {}).get("keep_intermediates", True)),
    )

    sweep, sweep_notes = _sweep_from_config(config, out_stream)
    alt_m, freq_hz = _grids(config)
    epoch = str(config.get("epoch", "2026-01-01T12:00:00Z"))

    outputs: list[Path] = []
    failures: dict[str, list[str]] = {}
    summaries: dict[str, Any] = {}

    for name in loaded.configurations:
        configuration = load_yaml(resolve_path(name, subdir="configurations"))
        result = pipeline.run_configuration(
            name,
            platform=_platform_for(configuration, epoch),
            sweep=sweep,
            freq_hz=freq_hz,
            alt_m=alt_m,
            epoch=epoch,
        )
        assessment = result.assessment
        if assessment is None:
            failures[name] = ["the pipeline returned no assessment"]
            continue

        text_path, json_path = write_report(assessment, out_dir, stem=name)
        outputs.extend([text_path, json_path])

        if result.ionogram is not None:
            npz = out_dir / f"{name}_ionogram.npz"
            np.savez_compressed(
                npz,
                freq_hz=result.ionogram.freq_hz,
                virtual_range_m=result.ionogram.virtual_range_m,
                snr_db=result.ionogram.snr_db,
                modes=np.asarray([str(m) for m in result.ionogram.modes]),
            )
            outputs.append(npz)

        record = to_dict(assessment)
        summaries[name] = {
            "summary": record["summary"],
            "timing_s": result.timing_s,
            "source": assessment.source,
        }
        failed = sorted(
            rid for rid, entry in assessment.requirements.items()
            if entry.get("verdict") == "fail"
        )
        if failed:
            failures[name] = failed

        counts = record["summary"]
        print(
            f"  {name:28s} {counts['goal']} goal, {counts['threshold']} threshold, "
            f"{counts['fail']} fail, {counts['not_evaluated']} not evaluated",
            file=out_stream,
        )

    tiers_record = dict(pipeline.describe())
    tiers_record.update(sweep_notes)
    manifest = Manifest(
        study=f"campaign/{loaded.name}",
        configuration=", ".join(loaded.configurations),
        config=dict(config),
        tiers=tiers_record,
        seed=loaded.seed,
        notes=(
            f"Campaign {loaded.name} from {loaded.source}.  "
            f"{len(loaded.configurations)} configuration(s), "
            f"{freq_hz.size} sounding frequencies, {alt_m.size} altitudes."
        ),
    )

    campaign_json = out_dir / "campaign.json"
    campaign_json.write_text(
        json.dumps(
            {
                "campaign": loaded.name,
                "configurations": list(loaded.configurations),
                "per_configuration": summaries,
                "failures": failures,
            },
            indent=2,
            default=str,
        ),
        encoding="utf-8",
    )
    outputs.append(campaign_json)
    manifest.outputs = [str(p) for p in outputs]
    manifest.write(out_dir)

    print(f"campaign {loaded.name} wrote to {out_dir}", file=out_stream)
    if failures:
        for name, ids in sorted(failures.items()):
            print(f"  {name}: failed {', '.join(ids)}", file=out_stream)
        return 1
    return 0


# Importing the study modules registers them.  Each one is numpy-only at module
# scope; pandas, scipy, matplotlib, docx and openpyxl are imported inside the
# functions that need them, so `import tss.studies` stays cheap.
from tss.studies import (  # noqa: E402,F401
    antenna_scale,
    doppler_separation,
    isis_validation,
    link_budget,
    min_frequency,
    route_one,
    sky_noise,
    trade_space,
    transcription,
)
