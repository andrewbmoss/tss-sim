"""Tier registry.

A campaign configuration names the implementation it wants for each layer, for
example ``l0: iri2020`` or ``l1: raytrace``.  This module maps those names onto
classes so that adding a new fidelity tier means writing a class and decorating
it, and does not mean editing the pipeline.

Usage::

    from tss.registry import register, build

    @register("l0", "chapman")
    class ChapmanIonosphere:
        ...

    tier = build("l0", "chapman", scale_height_km=60.0)

Registration happens as a side effect of importing the layer packages, which
``tss/__init__.py`` does.  A name that is registered but whose implementation is
still a stub will raise :class:`NotImplementedError` when it is called, not when
it is built, so a campaign can be validated before it is run.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any, TypeVar

__all__ = ["register", "build", "tier_class", "available", "LAYERS", "UnknownTier"]

LAYERS: tuple[str, ...] = ("l0", "l1", "l2", "l3", "l4", "l5")

_REGISTRY: dict[str, dict[str, type]] = {layer: {} for layer in LAYERS}

T = TypeVar("T", bound=type)


class UnknownTier(KeyError):
    """Raised when a configuration names a tier that is not registered."""


def register(layer: str, name: str) -> Callable[[T], T]:
    """Class decorator registering an implementation under ``layer``/``name``."""
    if layer not in _REGISTRY:
        raise ValueError(f"unknown layer {layer!r}, expected one of {LAYERS}")

    def _decorate(cls: T) -> T:
        if name in _REGISTRY[layer]:
            existing = _REGISTRY[layer][name]
            if existing is not cls:
                raise ValueError(
                    f"{layer}/{name} is already registered to {existing.__name__}"
                )
        _REGISTRY[layer][name] = cls
        cls.name = name
        cls.layer = layer
        return cls

    return _decorate


def tier_class(layer: str, name: str) -> type:
    """The class registered under ``layer``/``name``, without building it."""
    try:
        return _REGISTRY[layer][name]
    except KeyError as exc:
        known = ", ".join(sorted(_REGISTRY.get(layer, {}))) or "nothing"
        raise UnknownTier(
            f"no {layer} tier named {name!r}; registered: {known}"
        ) from exc


def build(layer: str, name: str, /, **kwargs: Any) -> Any:
    """Instantiate the implementation registered under ``layer``/``name``."""
    return tier_class(layer, name)(**kwargs)


def available(layer: str | None = None) -> dict[str, Iterable[str]]:
    """Return the registered tier names, for ``tss tiers`` and for tests."""
    if layer is not None:
        return {layer: sorted(_REGISTRY[layer])}
    return {lay: sorted(names) for lay, names in _REGISTRY.items()}
