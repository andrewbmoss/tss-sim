"""Layer: data
WP: WP1
Claims: -
Status: implemented
Requirements: -

Discovery of locally installed ionospheric models and of MATLAB.

Nothing here downloads anything.  PHaRLAP is not redistributable, NeQuick is
distributed under its own terms, MATLAB is licensed, and iri2016 needs a Fortran
compiler.  What this module does is answer one question honestly for each of
them: can L0 actually call it right now on this machine.  Having one place that
answers that means a tier does not have to carry its own detection code, and a
campaign can fail at configuration time with a message a person can act on
rather than at hour three of a run with an ImportError.

Availability is decided by a trial call, not by the presence of a file.  That is
deliberate and it is the whole point of this module.  An earlier version of the
author's MATLAB tooling tested for PHaRLAP with ``exist('raytrace_3d')`` and
rejected a working installation, because at the moment of the test the PHaRLAP
directories had not yet been added to the MATLAB path and ``exist`` answers a
question about the path, not a question about the disk.  The reverse failure is
just as easy: an ``iri2016`` package can be importable and still be unusable
because its compiled extension was built against a different numpy, and a
``matlab`` on PATH can be a stub from an uninstalled version.  A file test
answers "is something there".  A trial call answers "does it work", which is the
question a campaign actually needs answered.  Trials are therefore the default,
they are cached for the life of the process, and ``trial=False`` exists only for
the case where the caller wants a fast inventory and accepts a weaker answer.

Usage::

    from tss.data.local_models import probe, require

    status = probe()                  # what is available, with versions
    require("pharlap")                # raises with instructions if it is not
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = [
    "ModelStatus",
    "probe",
    "require",
    "pharlap_path",
    "nequick_path",
    "matlab_executable",
    "MissingModel",
    "DEFAULT_PHARLAP_PATH",
    "TRIAL_TIMEOUT_S",
]


DEFAULT_PHARLAP_PATH = r"C:\Users\Andrew\PhD\Chartier Verification\pharlap_4.7.4"
"""Where PHaRLAP lives on the author's machine.

``TSS_PHARLAP_PATH`` overrides it, and should be set on any other machine.
"""

PHARLAP_ENV = "TSS_PHARLAP_PATH"
NEQUICK_ENV = "TSS_NEQUICK_PATH"
MATLAB_ENV = "TSS_MATLAB"

TRIAL_TIMEOUT_S = 120.0
"""How long a trial call may take.

MATLAB's start-up alone is tens of seconds on a cold cache, so this is generous.
A trial that times out is reported as unavailable with the timeout as the
reason, which is a true statement about usability even when the install is fine.
"""

_PHARLAP_MARKERS: tuple[str, ...] = ("raytrace_3d", "iri2016", "gen_iono_grid_2d")
"""Functions a working PHaRLAP provides.  Used as the subject of the trial call."""

_probe_cache: dict[tuple[str, bool], ModelStatus] = {}


class MissingModel(RuntimeError):
    """Raised by :func:`require` when a model is not usable here."""


@dataclass(frozen=True)
class ModelStatus:
    """What was learned about one model.

    ``available`` is the answer to "can L0 call this now".  ``verified_by_trial``
    says whether that answer came from actually running something, so a caller
    can tell a confident yes from a hopeful one.  ``reason`` is always populated
    when ``available`` is False and is written to be read by a person.
    """

    name: str
    available: bool
    version: str = ""
    path: str = ""
    reason: str = ""
    verified_by_trial: bool = False
    detail: Mapping[str, Any] = field(default_factory=dict)

    def describe(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "available": self.available,
            "version": self.version,
            "path": self.path,
            "reason": self.reason,
            "verified_by_trial": self.verified_by_trial,
            "detail": dict(self.detail),
        }


# ---------------------------------------------------------------------------
# Locating things
# ---------------------------------------------------------------------------


def pharlap_path() -> Path | None:
    """The PHaRLAP root, from the environment or the default, if it exists."""
    candidates: list[str] = []
    env = os.environ.get(PHARLAP_ENV)
    if env:
        candidates.append(env)
    candidates.append(DEFAULT_PHARLAP_PATH)
    for c in candidates:
        p = Path(c).expanduser()
        if p.is_dir():
            return p
    return None


def nequick_path() -> Path | None:
    """The NeQuick 2 root, if one can be found.

    NeQuick 2 has no single canonical install layout.  What is looked for is a
    directory named by ``TSS_NEQUICK_PATH``, or a ``NeQuick`` directory beside
    PHaRLAP, that contains the modip and CCIR coefficient files the model cannot
    run without.
    """
    candidates: list[Path] = []
    env = os.environ.get(NEQUICK_ENV)
    if env:
        candidates.append(Path(env).expanduser())
    root = pharlap_path()
    if root is not None:
        candidates.extend(
            [root / "nequick", root / "NeQuick", root.parent / "NeQuick2", root.parent / "nequick2"]
        )
    for p in candidates:
        if p.is_dir():
            return p
    return None


def matlab_executable() -> str | None:
    """The MATLAB executable, from ``TSS_MATLAB`` or from PATH."""
    env = os.environ.get(MATLAB_ENV)
    if env:
        p = Path(env).expanduser()
        if p.is_file() and os.access(p, os.X_OK):
            return str(p)
        found = shutil.which(env)
        if found:
            return found
    return shutil.which("matlab")


# ---------------------------------------------------------------------------
# Trials
# ---------------------------------------------------------------------------


def _run(cmd: Sequence[str], *, timeout: float) -> tuple[int, str, str]:
    """Run a command and return ``(returncode, stdout, stderr)``.

    Never raises for a non-zero exit.  A missing executable, a timeout and a
    crash all come back as a return code and a message, because all three mean
    the same thing to the caller: not usable.
    """
    try:
        proc = subprocess.run(
            list(cmd),
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except FileNotFoundError as exc:
        return 127, "", str(exc)
    except subprocess.TimeoutExpired:
        return 124, "", f"timed out after {timeout:.0f} s"
    except OSError as exc:
        return 126, "", str(exc)
    return proc.returncode, proc.stdout or "", proc.stderr or ""


def _probe_matlab(trial: bool) -> ModelStatus:
    exe = matlab_executable()
    if exe is None:
        return ModelStatus(
            name="matlab",
            available=False,
            reason=(
                "no 'matlab' on PATH.  Add MATLAB's bin directory to PATH, or set "
                f"{MATLAB_ENV} to the executable."
            ),
        )
    if not trial:
        return ModelStatus(name="matlab", available=True, path=exe,
                           reason="", verified_by_trial=False)

    code, out, err = _run(
        [exe, "-batch", "fprintf('TSS_MATLAB_VERSION %s\\n', version);"],
        timeout=TRIAL_TIMEOUT_S,
    )
    if code != 0:
        return ModelStatus(
            name="matlab",
            available=False,
            path=exe,
            reason=f"'{exe} -batch' exited {code}: {(err or out).strip()[:300]}",
            verified_by_trial=True,
        )
    version = ""
    for line in out.splitlines():
        if line.startswith("TSS_MATLAB_VERSION"):
            version = line.split(None, 1)[-1].strip()
            break
    return ModelStatus(
        name="matlab", available=True, version=version, path=exe, verified_by_trial=True
    )


def _probe_pharlap(trial: bool) -> ModelStatus:
    root = pharlap_path()
    if root is None:
        return ModelStatus(
            name="pharlap",
            available=False,
            reason=(
                f"no PHaRLAP directory found.  Set {PHARLAP_ENV} to the unpacked "
                f"PHaRLAP root; the default looked for was {DEFAULT_PHARLAP_PATH}.  "
                "PHaRLAP is requested from DST Group and is not redistributable."
            ),
        )

    if not trial:
        return ModelStatus(
            name="pharlap",
            available=True,
            path=str(root),
            reason="directory present; not verified by trial",
            verified_by_trial=False,
            detail={"note": "a directory is not a working install, see the module docstring"},
        )

    matlab = _probe_matlab(trial=True)
    if not matlab.available:
        return ModelStatus(
            name="pharlap",
            available=False,
            path=str(root),
            reason=(
                f"PHaRLAP is at {root} but MATLAB is not usable, and PHaRLAP is a "
                f"MATLAB toolbox.  {matlab.reason}"
            ),
            verified_by_trial=True,
        )

    # The trial adds the PHaRLAP tree to the path first and only then asks
    # whether the functions resolve.  Asking before the addpath is what produced
    # the false negative described in the module docstring.
    marker = _PHARLAP_MARKERS[0]
    root_literal = str(root).replace("'", "''")
    script = (
        f"addpath(genpath('{root_literal}'));"
        f"if isempty(which('{marker}')), error('TSS:missing', "
        f"'{marker} did not resolve after addpath'); end;"
        f"fprintf('TSS_PHARLAP_OK %s\\n', which('{marker}'));"
    )
    env = dict(os.environ)
    exe = matlab.path
    try:
        proc = subprocess.run(
            [exe, "-batch", script],
            capture_output=True,
            text=True,
            timeout=TRIAL_TIMEOUT_S,
            check=False,
            env=env,
        )
        code, out, err = proc.returncode, proc.stdout or "", proc.stderr or ""
    except (OSError, subprocess.TimeoutExpired) as exc:
        code, out, err = 124, "", str(exc)

    if code != 0 or "TSS_PHARLAP_OK" not in out:
        return ModelStatus(
            name="pharlap",
            available=False,
            path=str(root),
            reason=(
                f"PHaRLAP at {root} did not answer a trial call: {(err or out).strip()[:300]}.  "
                "Check that the mex files match this platform and MATLAB version."
            ),
            verified_by_trial=True,
        )

    return ModelStatus(
        name="pharlap",
        available=True,
        version=_pharlap_version_from_path(root),
        path=str(root),
        verified_by_trial=True,
        detail={"matlab": matlab.version, "marker": marker},
    )


def _pharlap_version_from_path(root: Path) -> str:
    """Read the version out of the directory name, e.g. ``pharlap_4.7.4``.

    The directory name is the only version string PHaRLAP reliably exposes
    without running it, and it is only a label.  It is reported for the manifest
    and is not used to decide anything.
    """
    name = root.name
    if "_" in name:
        tail = name.rsplit("_", 1)[-1]
        if tail and tail[0].isdigit():
            return tail
    return ""


def _probe_iri2016(trial: bool) -> ModelStatus:
    try:
        import iri2016  # type: ignore
    except Exception as exc:
        return ModelStatus(
            name="iri2016",
            available=False,
            reason=(
                f"the iri2016 package is not importable: {exc}.  "
                "Install it with `pip install iri2016`, which needs a Fortran "
                "compiler, or use the MATLAB IRI that PHaRLAP bundles."
            ),
        )

    version = str(getattr(iri2016, "__version__", "") or "")
    path = str(getattr(iri2016, "__file__", "") or "")
    if not trial:
        return ModelStatus(
            name="iri2016", available=True, version=version, path=path,
            reason="importable; not verified by trial", verified_by_trial=False
        )

    # An importable iri2016 is not a working iri2016.  The package wraps a
    # compiled Fortran extension which can be absent, built for another platform,
    # or built against a different numpy, and none of that shows up until it is
    # called.  So call it.
    try:
        runner = getattr(iri2016, "IRI", None) or getattr(iri2016, "iri", None)
        if runner is None:
            raise AttributeError("no IRI entry point on the iri2016 module")
        result = runner("2020-01-01T12", (300, 300, 1), 40.0, -111.0)
        ok = result is not None
    except Exception as exc:
        return ModelStatus(
            name="iri2016",
            available=False,
            version=version,
            path=path,
            reason=(
                f"iri2016 imports but a trial call failed: {type(exc).__name__}: {exc}.  "
                "The compiled extension is usually the cause; rebuild it against this "
                "Python and numpy."
            ),
            verified_by_trial=True,
        )

    return ModelStatus(
        name="iri2016", available=bool(ok), version=version, path=path,
        verified_by_trial=True,
        reason="" if ok else "the trial call returned nothing",
    )


def _probe_nequick(trial: bool) -> ModelStatus:
    """Detect NeQuick 2.

    NeQuick is checked three ways because it arrives three ways: as a Python
    package, as a standalone Fortran build, and bundled inside PHaRLAP as the
    topside option of its IRI.  The third is the one the author actually uses,
    and it is reported as available only when PHaRLAP itself passed its trial,
    which keeps the answer honest.
    """
    try:
        import nequick  # type: ignore

        version = str(getattr(nequick, "__version__", "") or "")
        return ModelStatus(
            name="nequick2",
            available=True,
            version=version,
            path=str(getattr(nequick, "__file__", "") or ""),
            verified_by_trial=False,
            reason="python package importable; no standard trial entry point",
            detail={"flavour": "python package"},
        )
    except Exception:
        pass

    root = nequick_path()
    if root is not None:
        data_files = [p.name for p in root.glob("*") if p.suffix.lower() in {".asc", ".dat", ".txt"}]
        return ModelStatus(
            name="nequick2",
            available=bool(data_files),
            path=str(root),
            verified_by_trial=False,
            reason="" if data_files else f"{root} has no coefficient files in it",
            detail={"flavour": "standalone", "coefficient_files": len(data_files)},
        )

    pharlap = _cached(_probe_pharlap, trial)
    if pharlap.available:
        return ModelStatus(
            name="nequick2",
            available=True,
            path=pharlap.path,
            version=pharlap.version,
            verified_by_trial=pharlap.verified_by_trial,
            reason="",
            detail={
                "flavour": "via PHaRLAP",
                "note": (
                    "PHaRLAP's IRI offers NeQuick as a topside option; this is the "
                    "route the T1 comparison used."
                ),
            },
        )

    return ModelStatus(
        name="nequick2",
        available=False,
        reason=(
            f"no NeQuick found.  Set {NEQUICK_ENV} to a NeQuick 2 directory holding "
            "its modip and CCIR coefficient files, install a nequick Python package, "
            "or make PHaRLAP available and use its NeQuick topside option."
        ),
    )


_PROBES = {
    "pharlap": _probe_pharlap,
    "iri2016": _probe_iri2016,
    "nequick2": _probe_nequick,
    "matlab": _probe_matlab,
}


def _cached(fn: Any, trial: bool) -> ModelStatus:
    key = (fn.__name__, bool(trial))
    hit = _probe_cache.get(key)
    if hit is None:
        hit = fn(trial)
        _probe_cache[key] = hit
    return hit


def probe(
    names: Sequence[str] | None = None,
    *,
    trial: bool = True,
    refresh: bool = False,
) -> dict[str, dict[str, Any]]:
    """Report what is available, with versions where they can be had.

    ``trial=True`` runs a real call for each model, which is the default and is
    what the module docstring argues for.  It is slow, tens of seconds when
    MATLAB is involved, so results are cached for the life of the process;
    ``refresh=True`` clears that cache first.

    The returned mapping is JSON-serialisable and goes straight into a run
    manifest, so a result records not only which model was used but that the
    model answered a trial call on the machine that produced it.
    """
    if refresh:
        _probe_cache.clear()
    wanted = list(names) if names is not None else list(_PROBES)
    unknown = [n for n in wanted if n not in _PROBES]
    if unknown:
        raise KeyError(
            f"unknown model(s) {', '.join(unknown)}; known: {', '.join(sorted(_PROBES))}"
        )
    out: dict[str, dict[str, Any]] = {}
    for name in wanted:
        out[name] = _cached(_PROBES[name], trial).describe()
    out["_environment"] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "trial": trial,
    }
    return out


def require(name: str, *, trial: bool = True) -> ModelStatus:
    """Return the status of ``name``, or raise with something actionable.

    Called at configuration time by any tier that needs an external model, so
    that a campaign fails before it starts rather than partway through.
    """
    if name not in _PROBES:
        raise KeyError(
            f"unknown model {name!r}; known: {', '.join(sorted(_PROBES))}"
        )
    status = _cached(_PROBES[name], trial)
    if not status.available:
        raise MissingModel(
            f"{name} is not usable on this machine.  {status.reason}"
        )
    return status
