"""Import direction and private names across the package (code review F4, F11).

The layers are meant to depend only downward: L0 on nothing above it, L1 on L0,
and so on to L5, with the pipeline, the studies, the verification suite, the
figures and the command line on top, and the support modules (constants,
units, interfaces, registry, configuration, provenance, requirements, the sweep
law, the window tapers, the data fetchers, the platform) below every layer.
Until 2026-11-08 five imports ran the other way (L0 into L1, L2 into L3, L3
into L4 twice, and L3 into a study).  This test reads every module under
``src/tss`` and fails on any import, at module level or inside a function, of a
module that sits higher than the one importing it.

The second test holds the studies to public names: the helpers every study,
test and script uses to build a campaign's inputs were private
(``_sweep_from_config``, ``_grids``, ``_platform_for``, ``_base_config``) until
2026-11-08, and so were three helpers one study borrowed from another.
"""

from __future__ import annotations

import ast
from pathlib import Path

SRC = Path(__file__).resolve().parents[2] / "src" / "tss"

#: Rank of each top-level module or package; anything not listed is support (0).
RANK: dict[str, int] = {
    "l0_ionosphere": 1,
    "l1_propagation": 2,
    "l2_instrument": 3,
    "l3_signal": 4,
    "l4_retrieval": 5,
    "l5_assessment": 6,
    "pipeline": 7,
    "studies": 8,
    "verify": 8,
    "figures": 8,
    "cli": 9,
}


def _modules():
    for path in sorted(SRC.rglob("*.py")):
        rel = path.relative_to(SRC)
        if rel == Path("__init__.py"):
            continue  # the package root imports every layer to register the tiers
        yield rel, ast.parse(path.read_text(encoding="utf-8"))


def _top(rel: Path) -> str:
    return rel.parts[0][:-3] if len(rel.parts) == 1 else rel.parts[0]


def _targets(rel: Path, tree: ast.AST):
    """(line, dotted target under tss, names) for every tss import in ``tree``."""
    package = list(rel.parts[:-1])
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            for a in n.names:
                if a.name.startswith("tss."):
                    yield n.lineno, a.name[4:], ()
        elif isinstance(n, ast.ImportFrom):
            if n.level:
                base = package[: len(package) - (n.level - 1)]
                mod = ".".join(base + ([n.module] if n.module else []))
            elif n.module and (n.module == "tss" or n.module.startswith("tss.")):
                mod = n.module[4:]
            else:
                continue
            names = tuple(a.name for a in n.names)
            if mod == "":
                for name in names:  # from tss import sweep_law
                    yield n.lineno, name, ()
            else:
                yield n.lineno, mod, names


def test_no_module_imports_a_module_above_it():
    upward = []
    for rel, tree in _modules():
        mine = RANK.get(_top(rel), 0)
        for line, target, _ in _targets(rel, tree):
            theirs = RANK.get(target.split(".")[0], 0)
            if theirs > mine:
                upward.append(f"{rel.as_posix()}:{line} imports tss.{target}")
    assert not upward, "imports that run up the layer stack:\n" + "\n".join(upward)


def test_no_module_imports_a_private_name_from_another_package_or_a_study():
    private = []
    for rel, tree in _modules():
        mine = _top(rel)
        for line, target, names in _targets(rel, tree):
            theirs = target.split(".")[0]
            for name in names:
                if name.startswith("_") and not name.startswith("__") and (
                    theirs != mine or theirs == "studies"
                ):
                    private.append(f"{rel.as_posix()}:{line} imports tss.{target}.{name}")
    assert not private, "private names imported across packages or between studies:\n" + "\n".join(private)
