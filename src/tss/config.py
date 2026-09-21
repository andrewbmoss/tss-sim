"""Configuration loading, merging and validation.

A run is defined by a YAML campaign file.  The campaign names an instrument
configuration, a set of platform configurations to sweep, the tier
implementations to use for each layer, and the requirement set to score
against.  Nothing is hard-coded in a driver script.

Resolution order, later overriding earlier:

    1. ``config/default.yaml``
    2. any file named in the campaign's ``inherits`` list, in order
    3. the campaign file itself
    4. command-line overrides of the form ``--set a.b.c=value``

The resolved mapping is hashed into the run manifest, so a result can always be
traced back to the exact settings that produced it.
"""

from __future__ import annotations

import copy
import json
import re
from collections.abc import Iterable, Mapping, MutableMapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

__all__ = [
    "ConfigError",
    "load_campaign",
    "load_yaml",
    "deep_merge",
    "apply_overrides",
    "resolve_path",
    "repo_root",
    "config_dir",
    "Campaign",
]


class ConfigError(ValueError):
    """Raised for a malformed or self-contradictory configuration."""


def repo_root() -> Path:
    """The repository root, located from this file rather than from cwd."""
    return Path(__file__).resolve().parents[2]


def config_dir() -> Path:
    return repo_root() / "config"


def resolve_path(name: str | Path, *, subdir: str = "") -> Path:
    """Resolve a config reference.

    Accepts an absolute path, a path relative to the repository root, or a bare
    name to be found under ``config/<subdir>/`` with or without a ``.yaml``
    suffix.  This is what lets a campaign say ``instrument: baseline_6u``.
    """
    p = Path(name)
    if p.is_absolute() and p.exists():
        return p
    for candidate in (
        repo_root() / p,
        config_dir() / subdir / p,
        config_dir() / subdir / f"{p}.yaml",
        config_dir() / p,
        config_dir() / f"{p}.yaml",
    ):
        if candidate.exists():
            return candidate
    raise ConfigError(f"cannot resolve configuration reference {name!r}")


class _Loader(yaml.SafeLoader):
    """SafeLoader with a float resolver that accepts an unsigned exponent.

    PyYAML implements YAML 1.1, whose float resolver requires a sign on the
    exponent.  ``3.3e+5`` is a float and ``3.3e5`` is a string.  That is a trap
    for a configuration file full of physical quantities, and it fails late and
    confusingly: the value flows all the way to a comparison in the requirement
    scorer before anything complains, and the error names a type rather than
    the file and key that caused it.  Rather than require a sign on every
    exponent in every config file, this loader resolves the YAML 1.2 form as
    well.  The drawback is that the files are no longer strictly YAML 1.1, but
    they are read only by this package and the 1.2 behaviour is what anyone
    writing ``5.0e5`` intends.
    """


_FLOAT_1_2 = re.compile(
    r"""^(?:[-+]?(?:[0-9][0-9_]*)\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |\.[0-9][0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?[0-9][0-9_]*(?:[eE][-+]?[0-9]+)
        |[-+]?\.(?:inf|Inf|INF)
        |\.(?:nan|NaN|NAN))$""",
    re.X,
)
_Loader.add_implicit_resolver("tag:yaml.org,2002:float", _FLOAT_1_2, list("-+0123456789."))


def load_yaml(path: str | Path) -> dict[str, Any]:
    """Load a YAML file into a plain dict, with an empty file giving ``{}``."""
    text = Path(path).read_text(encoding="utf-8")
    data = yaml.load(text, Loader=_Loader)  # noqa: S506 - _Loader derives from SafeLoader
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"{path} must contain a mapping at the top level")
    return data


def deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    """Recursively merge ``override`` onto ``base`` without mutating either.

    Mappings merge key by key.  Every other type, lists included, replaces
    wholesale.  Replacing lists is deliberate: a campaign that lists three
    configurations means those three and not those three appended to the
    default.
    """
    out: dict[str, Any] = copy.deepcopy(dict(base))
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(out.get(key), Mapping):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def _coerce(text: str) -> Any:
    """Turn a command-line override string into a typed value."""
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return text


def apply_overrides(config: Mapping[str, Any], overrides: Iterable[str]) -> dict[str, Any]:
    """Apply ``a.b.c=value`` strings onto a configuration."""
    out: dict[str, Any] = copy.deepcopy(dict(config))
    for item in overrides:
        if "=" not in item:
            raise ConfigError(f"override {item!r} is not of the form key.path=value")
        dotted, _, raw = item.partition("=")
        node: MutableMapping[str, Any] = out
        parts = dotted.strip().split(".")
        for part in parts[:-1]:
            nxt = node.get(part)
            if not isinstance(nxt, MutableMapping):
                nxt = {}
                node[part] = nxt
            node = nxt
        node[parts[-1]] = _coerce(raw.strip())
    return out


@dataclass(frozen=True)
class Campaign:
    """A resolved campaign, ready to hand to the pipeline."""

    name: str
    config: Mapping[str, Any]
    source: Path

    @property
    def tiers(self) -> Mapping[str, Any]:
        return self.config.get("tiers", {})

    @property
    def configurations(self) -> list[str]:
        return list(self.config.get("configurations", []))

    @property
    def seed(self) -> int:
        return int(self.config.get("seed", 0))

    @property
    def requirements(self) -> str:
        return str(self.config.get("requirements", "requirements/requirements"))


def load_campaign(name: str | Path, overrides: Iterable[str] = ()) -> Campaign:
    """Load, merge and validate a campaign configuration."""
    path = resolve_path(name, subdir="campaigns")
    raw = load_yaml(path)

    default_path = config_dir() / "default.yaml"
    merged: dict[str, Any] = load_yaml(default_path) if default_path.exists() else {}

    for parent in raw.get("inherits", []) or []:
        merged = deep_merge(merged, load_yaml(resolve_path(parent)))

    merged = deep_merge(merged, raw)
    merged = apply_overrides(merged, overrides)
    merged.pop("inherits", None)

    _validate(merged, path)
    return Campaign(name=merged.get("name", path.stem), config=merged, source=path)


def _validate(config: Mapping[str, Any], source: Path) -> None:
    """Fail early and legibly on a campaign that cannot run."""
    from tss.registry import LAYERS, available

    tiers = config.get("tiers")
    if not isinstance(tiers, Mapping):
        raise ConfigError(f"{source}: 'tiers' must be a mapping of layer to tier name")

    known = available()
    for layer in LAYERS:
        if layer not in tiers:
            raise ConfigError(f"{source}: no tier selected for layer {layer!r}")
        chosen = tiers[layer]
        chosen_name = chosen if isinstance(chosen, str) else chosen.get("name")
        if chosen_name not in known[layer]:
            options = ", ".join(known[layer]) or "nothing"
            raise ConfigError(
                f"{source}: {layer} tier {chosen_name!r} is not registered; "
                f"available: {options}"
            )

    if not config.get("configurations"):
        raise ConfigError(
            f"{source}: at least one configuration must be named. "
            "Section 1.3 of the proposal reports results per configuration, so a "
            "run without one cannot be scored."
        )
