"""tss: a simulator for a low-SWaP FMCW topside ionospheric sounder.

The package implements the layered architecture of Section 3.4 of the
dissertation proposal, truth in at L0 and a scored density profile out at L5.

    >>> import tss
    >>> tss.available()["l0"]            # doctest: +SKIP
    ['analytic', 'chapman', 'iri2020', 'matlab_iri', 'tabulated']

Importing this package imports every layer package, which is what populates the
tier registry.  That import is kept cheap: a tier module may not import a heavy
or optional dependency at module scope.  PHaRLAP, iri2016 and the MATLAB engine
are all imported lazily inside the tier that needs them, so that ``import tss``
works on a machine that has none of them.
"""

from __future__ import annotations

__version__ = "0.1.0"

# Importing the layer packages registers their tiers.  Order matters only in
# that a layer may reference an interface from `interfaces`, never from a
# sibling layer.
from tss import (  # noqa: F401  # noqa: F401,E402
    constants,
    interfaces,
    l0_ionosphere,
    l1_propagation,
    l2_instrument,
    l3_signal,
    l4_retrieval,
    l5_assessment,
    units,
)
from tss.registry import available, build, register  # noqa: F401

__all__ = [
    "__version__",
    "available",
    "build",
    "register",
    "constants",
    "interfaces",
    "units",
]
