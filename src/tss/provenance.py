"""Run manifests.

The proposal promises that every figure is regenerable from a stored
configuration by the same driver that produces the corresponding table.  That
promise is only worth anything if each result carries enough information to
reproduce it, so every run writes a manifest next to its outputs.

A manifest records the git commit, whether the working tree was dirty, the
resolved configuration and its hash, the tier implementations selected, the RNG
seed, the package versions, and the external data files consumed with their
checksums.  A result whose manifest says the tree was dirty is a result that
cannot be reproduced, and :func:`check_clean` exists so a campaign intended for
publication can refuse to run in that state.
"""

from __future__ import annotations

import getpass
import hashlib
import json
import platform
import socket
import subprocess
import sys
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

__all__ = ["Manifest", "config_hash", "file_digest", "git_state", "check_clean"]


def config_hash(config: Mapping[str, Any]) -> str:
    """Stable SHA-256 over a configuration mapping."""
    blob = json.dumps(config, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def file_digest(path: str | Path, chunk: int = 1 << 20) -> str:
    """SHA-256 of a file, for recording which external data a run consumed."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def git_state(repo: str | Path | None = None) -> dict[str, Any]:
    """Return the commit, branch and dirty flag, or a reason they are unknown."""
    cwd = str(repo) if repo is not None else str(Path(__file__).resolve().parents[2])
    def _run(*args: str) -> str | None:
        try:
            return subprocess.check_output(
                ["git", *args], cwd=cwd, stderr=subprocess.DEVNULL, text=True
            ).strip()
        except (OSError, subprocess.CalledProcessError):
            return None

    commit = _run("rev-parse", "HEAD")
    if commit is None:
        return {"commit": None, "branch": None, "dirty": None,
                "note": "not a git working tree, or git is not on PATH"}
    status = _run("status", "--porcelain")
    return {
        "commit": commit,
        "branch": _run("rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": bool(status),
        "dirty_files": [line[3:] for line in (status or "").splitlines()][:50],
    }


@dataclass
class Manifest:
    """Everything needed to reproduce one run."""

    study: str
    configuration: str = ""
    config: Mapping[str, Any] = field(default_factory=dict)
    tiers: Mapping[str, Any] = field(default_factory=dict)
    seed: int | None = None
    inputs: list[dict[str, str]] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)
    notes: str = ""

    created_utc: str = field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(timespec="seconds")
    )
    git: Mapping[str, Any] = field(default_factory=git_state)
    environment: Mapping[str, Any] = field(
        default_factory=lambda: {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "host": socket.gethostname(),
            "user": getpass.getuser(),
        }
    )

    def add_input(self, path: str | Path, role: str = "") -> None:
        """Record an external data file and its checksum."""
        p = Path(path)
        self.inputs.append(
            {"path": str(p), "role": role,
             "sha256": file_digest(p) if p.is_file() else "",
             "bytes": str(p.stat().st_size) if p.is_file() else ""}
        )

    def write(self, directory: str | Path, filename: str = "manifest.json") -> Path:
        """Write the manifest as JSON and return the path."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        payload = asdict(self)
        payload["config_sha256"] = config_hash(self.config)
        out = directory / filename
        out.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str))
        return out


def check_clean(repo: str | Path | None = None, *, strict: bool = True) -> None:
    """Raise if the working tree is dirty.

    Campaigns whose output is destined for the dissertation call this with
    ``strict=True``.  Exploratory runs pass ``strict=False`` and get a manifest
    that honestly records the tree as dirty.
    """
    state = git_state(repo)
    if state.get("dirty") and strict:
        files = ", ".join(state.get("dirty_files", [])[:5])
        raise RuntimeError(
            "refusing to run a publication campaign from a dirty working tree; "
            f"uncommitted: {files}"
        )
