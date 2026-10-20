"""Layer: L0
WP: WP1
Claims: C1, C2
Status: partial
Requirements: INS-01

The ionosphere model backends, migrated from ``MinFreqStudy/iri_adapter.py``.

The study logic must not know which model produced a profile.  In the original
module every backend implemented one method,

    profiles(lat_deg, lon_deg, when, alt_grid_m, **opts) -> ne [n_samples, n_alt]

so that swapping IRI-2020's topside option, or substituting NeQuick2 later,
changed one constructor argument and nothing else.  That method is kept here
unchanged, because the band-edge study of WP0 calls it with a million sample
points at a time and a per-point interface would be the wrong shape for that
work.  What is added is :meth:`state_at`, the single-point
:class:`~tss.interfaces.L0Ionosphere` contract the layered pipeline uses, which
is implemented in terms of ``profiles`` so the two cannot diverge.

Backends
--------
``analytic``      alpha-Chapman topside with a prescribed diurnal and
                  latitudinal variation.  NOT a model of the real ionosphere.
                  It exists so that the harness can be tested end to end without
                  IRI and so a regression suite can run anywhere.  Never report
                  science results from it.
``iri2020``       IRI-2020 through PyIRI 0.1 (spherical-harmonic library), pure
                  Python and vectorised.  Needs a stated daily F10.7.
``iri_wrapper``   the ``iri2020`` / ``iri2016`` Fortran wrappers (Hirsch et al.),
                  one call per sample point.
``iri``           IRI-2020 with the implementation a configuration option,
                  ``model: pharlap`` (the default, ``matlab_iri``) or
                  ``model: pyiri`` (``iri2020``).  Campaigns should name this.
``matlab_iri``    IRI-2020 through the PHaRLAP mex, driven from MATLAB in batch
                  mode by ``matlab/tss/tss_iri_profiles.m`` and read back from
                  the CSV contract that file documents.  This is the likely path
                  on the author's machine, which has PHaRLAP.

Every backend records ``label`` and ``provenance``, both of which go into
``IonosphereState.source`` and into the run manifest.  A study whose figures do
not say which topside option produced them is not reproducible.

THE AVAILABILITY LESSON, WHICH IS THE REASON THIS MODULE IS CAREFUL
-------------------------------------------------------------------
An earlier version of this tooling decided whether PHaRLAP was usable with
``exist('iri2020', 'file')``.  That test rejected a working installation twice,
for two different reasons, and both are worth keeping in front of a reader.

The first is that ``exist`` answers a question about the MATLAB path, not a
question about the disk.  At the moment the test ran, the PHaRLAP directories
had not been added to the path, so a perfectly good install answered no.  The
second is that ``exist`` collapses a MEX file and an ``.m`` file into one return
code, and PHaRLAP ships an ``iri2020.m`` help stub in ``src/matlab`` beside the
mex, so even a test that waits for the ``addpath`` cannot distinguish the stub
from the real thing by return code alone.

The fix in both cases was the same: make a real trial call inside a try/except
and let the call itself settle it.  A trial also proves things no file test can,
notably that ``DIR_MODELS_REF_DAT`` points somewhere with coefficient files in
it and that the mex matches this platform and this MATLAB version.  That
approach now lives in :mod:`tss.data.local_models`, which is what this module
asks rather than carrying its own detection code, and it is repeated MATLAB-side
in the ``ensure_pharlap`` helper of ``tss_iri_profiles.m`` because a
``matlab -batch`` session inherits no interactive path.
"""

from __future__ import annotations

import datetime as dt
import json
import os
import subprocess
import tempfile
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from tss.constants import plasma_frequency
from tss.interfaces import IonosphereState
from tss.l0_ionosphere.base import (
    IonosphereBase,
    collision_frequency,
    default_te_k,
    iri_f107_from_r12,
    iri_r12_from_f107,
    parse_epoch,
    r12_to_f107,
    validate_alt_grid,
)
from tss.registry import register

__all__ = [
    "Backend",
    "AnalyticBackend",
    "PyIRIBackend",
    "IRIWrapperBackend",
    "MatlabIRIBackend",
    "IRIBackend",
    "IRI_MODELS",
    "IRI_DEFAULT_MODEL",
    "make_backend",
    "r12_to_f107",
    "IRI_TOPSIDE_OPTIONS",
    "EXTRA_COLUMNS",
    "PHARLAP_DEFAULT",
    "iri_index_span",
]

#: Topside models the PHaRLAP iri2020 wrapper actually exposes, spelled exactly
#: as ``iri_options.Ne_topside_model`` requires.  IRI itself also carries a
#: Gulyaeva option; this wrapper does not reach it, so the study spans three
#: models, not four.  NeQuick is the wrapper default, which is why running
#: "IRI vs NeQuick2" would have compared NeQuick largely to itself.
IRI_TOPSIDE_OPTIONS = ("NeQuick", "IRI-2001", "IRI-2001 corrected")

#: Column order of the per-point extras.  The MATLAB generator writes these with
#: the longer header names ``nmf2_m3, hmf2_km, solar_zenith_deg, dip_lat_deg,
#: modip_deg, spread_f_prob, rz12, f107_daily, kp, ap, ig12, f107_81day,
#: ap_daily, fof2_storm_ratio``, in this order.  The short names are kept
#: because the study code reads them.  Four were added on 2026-10-19 so that an
#: epoch-driven run records the indices IRI actually used, and the O+/H+
#: transition height (``h_oplus_hplus_km``, from IRI's ion composition) on
#: 2026-10-20 for SCI-02's goal.  A file written before then has ten or
#: fourteen and is read with the rest NaN.
EXTRA_COLUMNS = (
    "nmf2",
    "hmf2_km",
    "solar_zenith_deg",
    "dip_lat_deg",
    "modip_deg",
    "spreadf_prob",
    "rz12",
    "f107",
    "kp",
    "ap",
    "ig12",
    "f107_81",
    "ap_daily",
    "fof2_storm_ratio",
    "h_oplus_hplus_km",
)

#: Where PHaRLAP lives on the machine this study was developed against.  Passed
#: through to the MATLAB generator so a ``-batch`` session, which inherits
#: neither the path nor ``DIR_MODELS_REF_DAT``, can find the mex without relying
#: on its own guess list.  ``TSS_PHARLAP_PATH`` and ``PHARLAP_ROOT`` override it.
PHARLAP_DEFAULT = os.environ.get(
    "TSS_PHARLAP_PATH",
    os.environ.get("PHARLAP_ROOT", r"C:\Users\Andrew\PhD\Chartier Verification\pharlap_4.7.4"),
)


def iri_index_span(pharlap_root: str | None = None) -> tuple[dt.date, dt.date]:
    """First and last dates on which IRI's own index files can drive a run.

    An epoch-driven call (R12 = -1) reads R12 and IG12 from ``ig_rz.dat`` by
    month and F10.7 and ap from ``apf107.dat`` by day, and the PHaRLAP wrapper
    raises for a date outside either file.  The span is the overlap: from the
    later of the two first dates to the earlier of the two last.  ``ig_rz.dat``
    gives its span in its third line as ``month,year,month,year``; the last line
    of ``apf107.dat`` is its last day.
    """
    import calendar

    root = Path(pharlap_root or PHARLAP_DEFAULT) / "dat" / "iri2020"
    ig_lines = (root / "ig_rz.dat").read_text().splitlines()
    m0, y0, m1, y1 = (int(v) for v in ig_lines[2].strip().rstrip(",").split(",")[:4])
    ig_first = dt.date(y0, m0, 1)
    ig_last = dt.date(y1, m1, calendar.monthrange(y1, m1)[1])

    def _apf_date(line: str) -> dt.date:
        yy, mm, dd = int(line[0:3]), int(line[3:6]), int(line[6:9])
        return dt.date(1900 + yy if yy >= 58 else 2000 + yy, mm, dd)

    apf = [ln for ln in (root / "apf107.dat").read_text().splitlines() if ln.strip()]
    first = max(ig_first, _apf_date(apf[0]))
    last = min(ig_last, _apf_date(apf[-1]))
    return first, last


def _repo_root() -> Path:
    """Repository root, located from this file rather than from the cwd."""
    return Path(__file__).resolve().parents[3]


def _matlab_quote(text: str) -> str:
    """Escape a Python string for use inside a single-quoted MATLAB literal."""
    return str(text).replace("'", "''")


def _as_sequence(value: Any) -> list[Any]:
    """One value or many, always returned as a list."""
    if isinstance(value, (str, bytes)):
        return [value]
    if isinstance(value, Iterable) and not isinstance(value, dict):
        return list(value)
    return [value]


class Backend(IonosphereBase):
    """Base for the migrated backends.

    Carries the batch method the study calls and the single-point contract the
    pipeline calls.  ``profiles`` is abstract here exactly as it was in the
    original module; ``state_at`` is implemented once, in terms of ``profiles``,
    so that a backend author has one method to write and the two interfaces
    cannot disagree with each other.
    """

    label = "abstract"
    provenance = "none"

    def profiles(
        self,
        lat_deg: Sequence[float],
        lon_deg: Sequence[float],
        when: Sequence[dt.datetime],
        alt_grid_m: np.ndarray,
        **opts: Any,
    ) -> np.ndarray:
        """Electron density, m^-3, shaped ``[n_samples, n_alt]``."""
        raise NotImplementedError

    # -- the single-point contract ------------------------------------------

    def _extra_for_point(self, index: int) -> dict[str, float]:
        """Per-point diagnostics a backend chose to expose.  Empty by default."""
        del index
        return {}

    def state_at(
        self,
        lat_deg: float,
        lon_deg: float,
        epoch: str,
        alt_m: np.ndarray,
    ) -> IonosphereState:
        """One profile as an :class:`~tss.interfaces.IonosphereState`.

        The altitude grid is copied before use and the caller's array is never
        touched.  Electron temperature is not returned by ``profiles``, so it
        comes from :func:`tss.l0_ionosphere.base.default_te_k`, anchored on the
        peak height this profile actually has rather than on a fixed 300 km.
        The collision frequency follows from Ne and Te.  A backend that can
        supply a real Te should override this method; none of the migrated ones
        can, because the original study never asked for one.
        """
        alt = validate_alt_grid(alt_m)
        when = parse_epoch(epoch)
        ne = np.asarray(
            self.profiles([float(lat_deg)], [float(lon_deg)], [when], alt), dtype=float
        )
        ne = np.atleast_2d(ne)[0][: alt.size]
        if ne.size != alt.size:
            raise ValueError(
                f"{type(self).__name__}.profiles returned {ne.size} densities for a "
                f"grid of {alt.size} altitudes"
            )

        diagnostics = self.peak_diagnostics(alt, ne)
        diagnostics.update(self._extra_for_point(0))
        hmf2_m = float(diagnostics.get("hmf2_km", 300.0)) * 1.0e3
        te = default_te_k(alt, hmf2_m)
        nu = collision_frequency(ne, te)
        b = self.magnetic_field(lat_deg, lon_deg, alt, epoch)

        diagnostics["te_is_stand_in"] = True
        return IonosphereState(
            alt_m=alt,
            ne=ne,
            nu=np.asarray(nu, dtype=float),
            te=np.asarray(te, dtype=float),
            b_field_t=b,
            lat_deg=float(lat_deg),
            lon_deg=float(lon_deg),
            epoch=str(epoch),
            source=self.label,
            extra=diagnostics,
        )


# ---------------------------------------------------------------------------
@register("l0", "analytic")
class AnalyticBackend(Backend):
    """A stand-in, not a model.

    An alpha-Chapman topside whose peak density and height vary with solar
    zenith angle and latitude in roughly the right direction and by roughly the
    right amount.  Use it to exercise the harness; never to set a requirement.

    It differs from :class:`tss.l0_ionosphere.chapman.ChapmanIonosphere` in that
    the peak moves with location and time of day.  That makes it the right thing
    for testing an orbit-sampling harness, where a profile that does not change
    along the track would hide a bug in the sampling, and the wrong thing for a
    limiting-case check, where the closed-form tier is easier to reason about.

    ``seed`` is kept from the original signature.  Nothing in the density
    calculation draws from the generator it makes; the generator is held so that
    a perturbation layered on top has somewhere to draw from.  Passing ``rng``
    instead supplies that generator directly, which is what contract rule 4
    asks for and what a pipeline should do.
    """

    def __init__(
        self,
        f107: float = 120.0,
        plasma_scale_height: float = 100e3,
        seed: int = 0,
        *,
        rng: np.random.Generator | None = None,
        geomagnetic_model: str = "auto",
    ) -> None:
        super().__init__(geomagnetic_model=geomagnetic_model)
        self.f107 = f107
        self.hp = plasma_scale_height
        self.seed = seed
        self.rng = np.random.default_rng(seed) if rng is None else rng
        self.label = f"analytic-standin(F10.7={f107:.0f})"
        self.provenance = "AnalyticBackend -- NOT a science result"

    def _nmf2_hmf2(
        self,
        lat_deg: Sequence[float],
        lon_deg: Sequence[float],
        when: Sequence[dt.datetime],
    ) -> tuple[np.ndarray, np.ndarray]:
        """Peak density and height from solar zenith angle and latitude.

        Carried across from ``iri_adapter.py`` unchanged, including the
        constants: the 0.02 floor on the zenith cosine that keeps the night side
        from collapsing to zero, the equatorial ionisation anomaly bump centred
        on 15 degrees magnetic with a 12 degree width, and the 4e10 m^-3 floor on
        the peak density.
        """
        lat = np.radians(np.asarray(lat_deg, dtype=float))
        doy = np.array([w.timetuple().tm_yday for w in np.atleast_1d(when)])
        ut = np.array([w.hour + w.minute / 60.0 for w in np.atleast_1d(when)])
        lt = (ut + np.asarray(lon_deg, dtype=float) / 15.0) % 24.0
        decl = np.radians(23.44) * np.sin(2 * np.pi * (doy - 80) / 365.25)
        ha = np.radians(15.0 * (lt - 12.0))
        cosx = np.sin(lat) * np.sin(decl) + np.cos(lat) * np.cos(decl) * np.cos(ha)
        cosx = np.clip(cosx, 0.02, 1.0)
        solar = 0.4 + 0.6 * (self.f107 - 70.0) / 110.0
        eia = 1.0 + 0.55 * np.exp(-((np.degrees(np.abs(lat)) - 15.0) / 12.0) ** 2)
        nmf2 = 1.1e12 * solar * eia * cosx**0.55
        nmf2 = np.maximum(nmf2, 4e10)
        hmf2 = 280e3 + 60e3 * (1.0 - cosx) + 30e3 * np.cos(2 * lat)
        return nmf2, hmf2

    def profiles(
        self,
        lat_deg: Sequence[float],
        lon_deg: Sequence[float],
        when: Sequence[dt.datetime],
        alt_grid_m: np.ndarray,
        **opts: Any,
    ) -> np.ndarray:
        nmf2, hmf2 = self._nmf2_hmf2(lat_deg, lon_deg, when)
        h = np.asarray(alt_grid_m, dtype=float)[None, :]
        hchap = self.hp / 2.0
        z = (h - hmf2[:, None]) / hchap
        return nmf2[:, None] * np.exp(0.5 * (1.0 - z - np.exp(-z)))

    def _settings(self) -> dict[str, Any]:
        return {
            "f107": float(self.f107),
            "plasma_scale_height_m": float(self.hp),
            "seed": self.seed,
            "science_grade": False,
        }


# ---------------------------------------------------------------------------
@register("l0", "iri2020")
class PyIRIBackend(Backend):
    """IRI-2020 through PyIRI: pure Python, vectorised, and the default IRI route.

    PyIRI (Forsythe et al. 2023, Space Weather; major update Servan-Schreiber
    et al. 2026, Space Weather 24, e2025SW004853) is an independent Python
    implementation of the International Reference Ionosphere, not a wrapper of
    the Fortran code.  This tier calls its spherical-harmonic library,
    ``PyIRI.sh_library.IRI_density_1day``, which is the 0.1 series' main
    entry point.  The options default to IRI-2020's own defaults:

    * foF2 from the URSI coefficients (``foF2_coeff="URSI"``; ``"CCIR"`` also);
    * hmF2 from Shubin 2015 (``hmF2_model="SHU2015"``; ``"AMTB2013"`` and
      ``"BSE1979"`` also);
    * the bottomside from ABT-2009 (B0, B1), which PyIRI fixes;
    * the topside an Epstein layer with the NeQuick-2 thickness law and its
      growing scale height (g = 0.125, r = 100), which PyIRI also fixes.

    So ``topside`` can only be ``"NeQuick"``.  IRI's other topside options
    (IRI-2001, IRI-2001 corrected) are reached through ``matlab_iri``.  Being a
    reimplementation, PyIRI agrees with Fortran IRI-2020 closely but not to the
    last digit; ``docs/pyiri-validation.md`` records the comparison.

    Solar activity is the daily F10.7 the caller gives, ``f107`` in SFU.  It is
    required: PyIRI takes no epoch-driven index, and a default would put a
    solar level into every result without anyone choosing it.  A campaign
    states it, for example ``l0: {name: iri2020, f107: 120}``.

    The import is inside the constructor, so ``import tss`` does not need
    PyIRI and a campaign that does not name this tier never pays for it.
    """

    #: Topside models this tier can give.  PyIRI builds one, NeQuick-2's.
    TOPSIDES: tuple[str, ...] = ("NeQuick",)
    FOF2_COEFFS: tuple[str, ...] = ("URSI", "CCIR")
    HMF2_MODELS: tuple[str, ...] = ("SHU2015", "AMTB2013", "BSE1979")

    def __init__(
        self,
        topside: str = "NeQuick",
        f107: float | None = None,
        *,
        foF2_coeff: str = "URSI",
        hmF2_model: str = "SHU2015",
        geomagnetic_model: str = "auto",
    ) -> None:
        super().__init__(geomagnetic_model=geomagnetic_model)
        if topside not in self.TOPSIDES:
            raise ValueError(
                f"PyIRI builds only the NeQuick-2 topside; got topside={topside!r}.  "
                "Use l0: matlab_iri for IRI's IRI-2001 topside options."
            )
        if foF2_coeff not in self.FOF2_COEFFS:
            raise ValueError(f"foF2_coeff must be one of {self.FOF2_COEFFS}, got {foF2_coeff!r}")
        if hmF2_model not in self.HMF2_MODELS:
            raise ValueError(f"hmF2_model must be one of {self.HMF2_MODELS}, got {hmF2_model!r}")
        if f107 is None:
            raise ValueError(
                "the iri2020 tier needs a daily F10.7 (SFU): state it in the campaign, "
                "e.g. l0: {name: iri2020, f107: 120}.  PyIRI has no epoch-driven index."
            )
        if not 50.0 <= float(f107) <= 400.0:
            raise ValueError(f"f107 = {f107} SFU is outside 50..400")

        from importlib.metadata import version

        import PyIRI
        import PyIRI.sh_library as sh

        self.sh = sh
        self.coeff_dir = PyIRI.coeff_dir
        self.pyiri_version = version("PyIRI")
        self.topside = topside
        self.f107 = float(f107)
        self.foF2_coeff = foF2_coeff
        self.hmF2_model = hmF2_model
        self._last_extra: list[dict[str, float]] = []
        self.label = (
            f"IRI-2020 via PyIRI {self.pyiri_version}(F10.7={self.f107:g}, "
            f"foF2={foF2_coeff}, hmF2={hmF2_model}, topside={topside})"
        )
        self.provenance = (
            f"PyIRI {self.pyiri_version} sh_library.IRI_density_1day; "
            f"coefficients {self.coeff_dir}"
        )

    def profiles(
        self,
        lat_deg: Sequence[float],
        lon_deg: Sequence[float],
        when: Sequence[dt.datetime],
        alt_grid_m: np.ndarray,
        **opts: Any,
    ) -> np.ndarray:
        """Electron density, m^-3, shaped ``[n_samples, n_alt]``.

        Samples are grouped by date and universal time, because PyIRI evaluates
        a set of points at one day and a set of times; each group is one call.
        ``f107`` in ``opts`` overrides the tier's value for this call only.
        """
        lat = np.asarray(lat_deg, dtype=float).ravel()
        lon = np.asarray(lon_deg, dtype=float).ravel()
        times = list(when) if not isinstance(when, dt.datetime) else [when]
        if not lat.size == lon.size == len(times):
            raise ValueError("lat_deg, lon_deg and when must be parallel")
        alt_km = np.asarray(alt_grid_m, dtype=float).ravel() / 1.0e3
        f107 = float(opts.get("f107", self.f107))

        groups: dict[tuple[int, int, int, float], list[int]] = {}
        for i, w in enumerate(times):
            ut = w.hour + w.minute / 60.0 + (w.second + w.microsecond / 1e6) / 3600.0
            groups.setdefault((w.year, w.month, w.day, ut), []).append(i)

        out = np.empty((lat.size, alt_km.size), dtype=float)
        extra: list[dict[str, float]] = [{} for _ in range(lat.size)]
        for (year, month, day, ut), idx in groups.items():
            ix = np.asarray(idx, dtype=int)
            # Longitude to (-180, 180]: the SH basis is periodic, but the
            # local-time terms are not guaranteed to be.
            lon_g = (lon[ix] + 180.0) % 360.0 - 180.0
            F2, F1, E, _es, _sun, mag, edp = self.sh.IRI_density_1day(
                year, month, day, np.array([ut]), lon_g, lat[ix], alt_km, f107,
                coeff_dir=self.coeff_dir, foF2_coeff=self.foF2_coeff,
                hmF2_model=self.hmF2_model, coord="GEO", old_output=False,
            )
            # edp is (N_T, N_V, N_G) with N_T = 1.
            out[ix] = np.asarray(edp, dtype=float)[0].T
            for k, i in enumerate(ix):
                extra[i] = {
                    "nmf2": float(F2["Nm"][0, k]),
                    "hmf2_km": float(F2["hm"][0, k]),
                    "fof2_hz": float(F2["fo"][0, k]) * 1.0e6,
                    "b_top_km": float(F2["B_top"][0, k]),
                    "b0_km": float(F2["B0"][0, k]),
                    "b1": float(F2["B1"][0, k]),
                    "modip_deg": float(np.asarray(mag["modip"]).ravel()[k]),
                    "f107": f107,
                }
        self._last_extra = extra
        return out

    def _extra_for_point(self, index: int) -> dict[str, float]:
        """The model's own NmF2, hmF2 and thicknesses, not the grid's peak."""
        if 0 <= index < len(self._last_extra):
            return dict(self._last_extra[index])
        return {}

    def _settings(self) -> dict[str, Any]:
        return {
            "topside": self.topside,
            "f107": self.f107,
            "foF2_coeff": self.foF2_coeff,
            "hmF2_model": self.hmF2_model,
            "pyiri_version": self.pyiri_version,
        }


# ---------------------------------------------------------------------------
@register("l0", "iri_wrapper")
class IRIWrapperBackend(Backend):
    """The ``iri2020`` / ``iri2016`` Fortran wrappers.  One call per sample point.

    The module is imported by name in the constructor so that the choice between
    ``iri2016`` and ``iri2020`` is a configuration value rather than an import
    statement, and so that ``import tss`` needs neither.

    Note that :mod:`tss.data.local_models` probes ``iri2016`` with a trial call
    rather than with an import, for the reason the module docstring gives: an
    importable wrapper whose compiled extension was built against a different
    numpy imports cleanly and then fails on the first call.  A campaign that
    names this tier should call ``require("iri2016")`` at configuration time.
    """

    def __init__(
        self,
        topside: str = "NeQuick",
        module: str = "iri2016",
        *,
        geomagnetic_model: str = "auto",
    ) -> None:
        super().__init__(geomagnetic_model=geomagnetic_model)
        self.mod = __import__(module)
        self.module_name = module
        self.topside = topside
        self.label = f"{module}(topside={topside})"
        self.provenance = f"{module} Fortran wrapper"

    def profiles(
        self,
        lat_deg: Sequence[float],
        lon_deg: Sequence[float],
        when: Sequence[dt.datetime],
        alt_grid_m: np.ndarray,
        **opts: Any,
    ) -> np.ndarray:
        alt_km = np.asarray(alt_grid_m) / 1e3
        rng = (alt_km[0], alt_km[-1], alt_km[1] - alt_km[0])
        out = np.empty((len(lat_deg), len(alt_km)))
        for i, (la, lo, w) in enumerate(zip(lat_deg, lon_deg, when, strict=False)):
            ds = self.mod.IRI(w, rng, float(la), float(lo))
            out[i] = np.asarray(ds["ne"]).ravel()[: len(alt_km)]
        return out

    def _settings(self) -> dict[str, Any]:
        return {"topside": self.topside, "module": self.module_name}


# ---------------------------------------------------------------------------
@register("l0", "matlab_iri")
class MatlabIRIBackend(Backend):
    """Batch bridge to the PHaRLAP iri2020 mex, through MATLAB.

    The MATLAB side is ``matlab/tss/tss_iri_profiles.m``, which wraps

        [iono, iono_extra] = iri2020(lat, lon, R12, UT, ht_start, ht_step,
                                     num_hts, iri_options)

    with ``iono(1,:)`` the electron density in m^-3, and writes the three-CSV
    contract that file documents.  This class writes the point list and the
    altitude grid, invokes MATLAB in batch mode, and reads the contract back.

    WHY R12 RATHER THAN AN EPOCH
    ----------------------------
    Solar activity is set by a user R12 rather than by epoch.  Passing R12 = -1
    would make IRI read its indices from ``ig_rz.dat`` by date, and the wrapper
    raises if the epoch falls outside that file's span, which would abort a long
    batch partway through.  A user R12 removes that dependency and gives
    reproducible control over solar level, at the documented cost of the foF2
    storm model being forced off.  Storm sensitivity is a separate run against
    real epochs.

    EPOCH MODE (added 2026-10-19)
    -----------------------------
    ``solar_mode="epoch"`` passes R12 = -1, so IRI runs on each sample's own
    R12, IG12, daily and 81-day F10.7 and ap from its index files, and the
    storm models follow ``iri_options`` (PHaRLAP turns the foF2 storm model on
    by default in this mode; a campaign should say which it wants).  The abort
    that kept this route closed is prevented rather than risked:
    :meth:`profiles` checks every epoch against :func:`iri_index_span` before
    MATLAB starts and raises naming the first epoch outside it.  The reason it
    is needed is the solar-minimum finding of 2026-10-18: with a user R12 the
    mex derives IG12 from R12 by formula (9.1 at R12 = 15), while the index
    file gives 17.3 for the solar-minimum ISIS-1 epochs, and the densities move
    by about 16 per cent.

    AVAILABILITY
    ------------
    Both MATLAB and PHaRLAP are checked with :func:`tss.data.local_models.require`
    at construction, which makes a real trial call rather than testing for a
    file.  The module docstring explains why that distinction cost the author
    time.  Pass ``probe=False`` to skip the check, which is what a unit test
    that only exercises the file handling wants, and which is not what a
    campaign wants.

    ``extra`` holds the per-point ``iono_extra`` quantities from the last call,
    in :data:`EXTRA_COLUMNS` order.  NmF2 and hmF2 from there are more accurate
    than reading the peak off the sampled altitude grid, so
    :meth:`_extra_for_point` prefers them.
    """

    def __init__(
        self,
        topside: str = "NeQuick",
        r12: float = 10.0,
        matlab: str = "matlab",
        script_dir: str | None = None,
        workdir: str | None = None,
        iri_options: dict[str, Any] | None = None,
        timeout: float = 14400,
        pharlap_root: str | None = None,
        *,
        study: str = "iri",
        tag: str = "",
        progress: int = 2000,
        manifest: bool = False,
        probe: bool = True,
        geomagnetic_model: str = "auto",
        solar_mode: str = "r12",
    ) -> None:
        super().__init__(geomagnetic_model=geomagnetic_model)
        if topside not in IRI_TOPSIDE_OPTIONS:
            raise ValueError(
                f"topside must be one of {IRI_TOPSIDE_OPTIONS}, got {topside!r}"
            )
        if solar_mode not in ("r12", "epoch"):
            raise ValueError(f"solar_mode must be 'r12' or 'epoch', got {solar_mode!r}")
        if solar_mode == "epoch":
            r12 = -1.0
        elif not (1.0 <= r12 <= 200.0):
            raise ValueError(
                "R12 must be in 1..200 for a user-specified solar level; use "
                "solar_mode='epoch' for IRI's own indices -- see the docstring"
            )
        self.solar_mode = solar_mode
        self.topside = topside
        self.r12 = float(r12)
        self.matlab = matlab
        self.script_dir = script_dir or str(_repo_root() / "matlab" / "tss")
        self.workdir = workdir or tempfile.mkdtemp(prefix="tss_iri_")
        self.iri_options = iri_options or {}
        self.pharlap_root = pharlap_root or PHARLAP_DEFAULT
        self.timeout = timeout
        self.study = study
        self.tag = tag or (
            self.default_tag(topside, self.r12) if solar_mode == "r12"
            else f"{str(topside).replace(' ', '-').replace('.', '').lower()}_epoch"
        )
        self.progress = int(progress)
        self.manifest = bool(manifest)
        self.extra: np.ndarray | None = None
        self.last_info: dict[str, Any] = {}
        if solar_mode == "epoch":
            self.label = f"IRI-2020(topside={topside}, epoch indices)"
            self.provenance = (
                f"PHaRLAP iri2020 mex via {self.script_dir}/tss_iri_profiles.m; "
                "R12 = -1: R12, IG12, F10.7 and ap from ig_rz.dat and apf107.dat by "
                f"epoch; storm models per iri_options {json.dumps(self.iri_options)}"
            )
        else:
            self.label = f"IRI-2020(topside={topside}, R12={r12:g})"
            self.provenance = (
                f"PHaRLAP iri2020 mex via {self.script_dir}/tss_iri_profiles.m; "
                f"R12={r12:g} -> F10.7={r12_to_f107(r12):.1f}; "
                f"foF2 storm model off (forced by user R12)"
            )
        self.availability: dict[str, Any] = {}
        if probe:
            self._require_local_models()

    # -- availability --------------------------------------------------------

    def _require_local_models(self) -> None:
        """Ask :mod:`tss.data.local_models` whether this route can run at all.

        The probing lives there and is not reimplemented here.  It makes a real
        trial call for both MATLAB and PHaRLAP, which is the point: a file test
        answers "is something there" and a trial call answers "does it work",
        and only the second is the question a campaign needs answered.  The
        import is inside the method so that building this class never drags the
        probe machinery in when ``probe=False``.
        """
        from tss.data.local_models import require

        matlab_status = require("matlab")
        pharlap_status = require("pharlap")
        self.availability = {
            "matlab": matlab_status.describe(),
            "pharlap": pharlap_status.describe(),
        }

    # -- file names ----------------------------------------------------------

    @staticmethod
    def default_tag(topside: str, r12: float) -> str:
        """Reproduce the default tag ``tss_iri_profiles.m`` builds.

        MATLAB side:
        ``sprintf('%s_r12-%g', lower(strrep(strrep(topside, ' ', '-'), '.', '')), R12)``
        so ``NeQuick`` at R12 = 105 gives ``nequick_r12-105`` and
        ``IRI-2001 corrected`` gives ``iri-2001-corrected_r12-105``.
        """
        stem = str(topside).replace(" ", "-").replace(".", "").lower()
        return f"{stem}_r12-{_matlab_g(r12)}"

    def _expected(self, what: str, ext: str = "csv") -> Path:
        return Path(self.workdir) / f"{self.study}_{what}_{self.tag}.{ext}"

    def _find_output(self, what: str) -> Path:
        """Locate one output file, falling back to a glob on the stem.

        The exact name is built from the same rule the MATLAB side uses, so it
        normally matches.  The glob exists because that rule involves MATLAB's
        ``%g`` formatting of R12 and a lower-cased, punctuation-stripped topside
        name, and a mismatch there would otherwise present as a missing file
        with no hint about why.  When the glob finds exactly one candidate it is
        used and the discrepancy is recorded in ``last_info``.
        """
        exact = self._expected(what)
        if exact.is_file():
            return exact
        candidates = sorted(Path(self.workdir).glob(f"{self.study}_{what}_*.csv"))
        if len(candidates) == 1:
            self.last_info.setdefault("name_mismatch", {})[what] = {
                "expected": str(exact),
                "found": str(candidates[0]),
            }
            return candidates[0]
        raise FileNotFoundError(
            f"tss_iri_profiles did not leave a readable {what} file in {self.workdir}.  "
            f"Expected {exact.name}; found {[c.name for c in candidates] or 'nothing'}."
        )

    # -- the batch method ----------------------------------------------------

    def profiles(
        self,
        lat_deg: Sequence[float],
        lon_deg: Sequence[float],
        when: Sequence[dt.datetime],
        alt_grid_m: np.ndarray,
        **opts: Any,
    ) -> np.ndarray:
        """Electron density for every sample point, m^-3, ``[n_samples, n_alt]``.

        The grid checks are the ones the original carried and the ones the mex
        actually imposes: a uniformly spaced grid and between 2 and 1000
        samples.  They are made here as well as MATLAB-side so that a bad grid
        fails in a second rather than after a MATLAB start-up.
        """
        alt_km = np.asarray(alt_grid_m, dtype=float) / 1e3
        step = np.diff(alt_km)
        if step.size == 0 or not np.allclose(step, step[0]):
            raise ValueError("iri2020 needs a uniformly spaced altitude grid")
        if not (2 <= len(alt_km) <= 1000):
            raise ValueError(f"num_hts must be 2..1000, got {len(alt_km)}")

        lats = _as_sequence(lat_deg)
        lons = _as_sequence(lon_deg)
        epochs = _as_sequence(when)
        if not (len(lats) == len(lons) == len(epochs)):
            raise ValueError(
                "lat_deg, lon_deg and when must be the same length; got "
                f"{len(lats)}, {len(lons)}, {len(epochs)}"
            )
        if self.solar_mode == "epoch":
            self._check_epochs_in_span(epochs)

        os.makedirs(self.workdir, exist_ok=True)
        rows = [
            [w.year, w.month, w.day, w.hour, w.minute, float(la), float(lo)]
            for la, lo, w in zip(lats, lons, epochs, strict=False)
        ]
        points_csv = os.path.join(self.workdir, "points.csv")
        alt_csv = os.path.join(self.workdir, "alt_km.csv")
        np.savetxt(points_csv, np.array(rows, dtype=float), delimiter=",", fmt="%.10g")
        np.savetxt(alt_csv, alt_km, delimiter=",", fmt="%.10g")

        cmd = self._matlab_command(points_csv, alt_csv)

        # Capture rather than inherit.  A MATLAB error under -batch goes to the
        # child's stderr and is otherwise lost behind CalledProcessError, which
        # reports only the exit code -- the single most common way to waste an
        # hour on this binding.
        proc = subprocess.run(
            [self.matlab, "-batch", cmd],
            timeout=self.timeout,
            cwd=self.workdir,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                "MATLAB exited {} running tss_iri_profiles.\n"
                "--- command ---\n{}\n"
                "--- stdout ---\n{}\n"
                "--- stderr ---\n{}".format(
                    proc.returncode,
                    cmd,
                    (proc.stdout or "").strip(),
                    (proc.stderr or "").strip(),
                )
            )

        self.last_info = {
            "command": cmd,
            "stdout_tail": (proc.stdout or "").strip()[-2000:],
            "n_points": len(rows),
            "n_alt": int(alt_km.size),
        }
        ne = self._read_profiles(len(rows), int(alt_km.size))
        self.extra = self._read_extra(len(rows))
        return ne

    def _matlab_command(self, points_csv: str, alt_csv: str) -> str:
        """The single ``-batch`` command string.

        ``addpath`` and nothing else.  ``restoredefaultpath`` is never used
        anywhere in this package: it would throw away whatever the operator has
        deliberately put on the path, PHaRLAP included, and the failure it causes
        surfaces far from its cause.
        """
        # The options go to MATLAB as a JSON file, not as a JSON string: the
        # double quotes of a string on the -batch command line are removed by
        # Windows argument parsing, and jsondecode then fails (2026-10-19).
        opts_json = ""
        if self.iri_options:
            opts_path = os.path.join(self.workdir, "iri_options.json")
            with open(opts_path, "w", encoding="utf-8") as fh:
                json.dump(self.iri_options, fh)
            opts_json = opts_path
        parts = [
            f"addpath('{_matlab_quote(self.script_dir)}');",
            "tss_iri_profiles(",
            f"'{_matlab_quote(self.workdir)}',",
            f"'{_matlab_quote(points_csv)}',",
            f"'{_matlab_quote(alt_csv)}',",
            f"'Topside', '{_matlab_quote(self.topside)}',",
            f"'R12', {_matlab_g(self.r12)},",
            f"'IriOptions', '{_matlab_quote(opts_json)}',",
            f"'PharlapRoot', '{_matlab_quote(self.pharlap_root)}',",
            f"'Study', '{_matlab_quote(self.study)}',",
            f"'Tag', '{_matlab_quote(self.tag)}',",
            f"'Progress', {self.progress:d},",
            f"'Manifest', {'true' if self.manifest else 'false'}",
            "); exit",
        ]
        return " ".join(parts)

    def _read_profiles(self, n_points: int, n_alt: int) -> np.ndarray:
        """Read the density block out of the profiles CSV.

        The contract is a header row then ``point_id`` and seven epoch and
        location columns, then ``n_alt`` density columns in grid order.  The
        eight leading columns are dropped here rather than parsed, because the
        caller already knows the points it asked for and re-reading them would
        create a second copy of the point list that can drift out of step with
        the first.
        """
        path = self._find_output("profiles")
        table = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
        if table.shape[0] != n_points:
            raise ValueError(
                f"{path.name} has {table.shape[0]} rows for {n_points} points asked for"
            )
        if table.shape[1] < 8 + n_alt:
            raise ValueError(
                f"{path.name} has {table.shape[1]} columns, expected at least {8 + n_alt}"
            )
        return np.asarray(table[:, 8 : 8 + n_alt], dtype=float)

    def _read_extra(self, n_points: int) -> np.ndarray | None:
        """Read the per-point extras, or ``None`` if the generator wrote none."""
        try:
            path = self._find_output("extra")
        except FileNotFoundError:
            return None
        table = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
        if table.shape[0] != n_points or table.shape[1] < 8 + 10:
            return np.asarray(table, dtype=float)
        # A file from before 2026-10-19 carries ten extras, and one from before
        # 2026-10-20 fourteen; pad the rest with NaN
        # so that a column is never read from the wrong position.
        out = np.full((table.shape[0], len(EXTRA_COLUMNS)), np.nan)
        k = min(len(EXTRA_COLUMNS), table.shape[1] - 8)
        out[:, :k] = table[:, 8 : 8 + k]
        return out

    def _extra_for_point(self, index: int) -> dict[str, float]:
        """IRI's own NmF2 and hmF2 and the rest, when the extras were written.

        These override the values read off the sampled grid, because IRI knows
        where its own peak is and a 5 km grid does not.  ``-1`` has already been
        turned into NaN MATLAB-side, so a NaN here means IRI declined to answer
        and not that the value is zero.
        """
        if self.extra is None:
            return {}
        arr = np.atleast_2d(self.extra)
        if index >= arr.shape[0] or arr.shape[1] < len(EXTRA_COLUMNS):
            return {}
        row = arr[index]
        out = {name: float(row[i]) for i, name in enumerate(EXTRA_COLUMNS)}
        if np.isfinite(out.get("nmf2", np.nan)):
            out["fof2_hz"] = float(plasma_frequency(out["nmf2"]))
        return out

    def _check_epochs_in_span(self, epochs: Sequence[dt.datetime]) -> None:
        """Raise before MATLAB starts if any epoch is outside the index files."""
        first, last = iri_index_span(self.pharlap_root)
        days = [w.date() if isinstance(w, dt.datetime) else w for w in epochs]
        bad = [d for d in days if not (first <= d <= last)]
        if bad:
            raise ValueError(
                f"{len(bad)} of {len(days)} epochs fall outside IRI's index files "
                f"({first} to {last}), for example {bad[0]}; an epoch-driven run "
                "would abort there"
            )

    def _settings(self) -> dict[str, Any]:
        epoch = self.solar_mode == "epoch"
        return {
            "topside": self.topside,
            "solar_mode": self.solar_mode,
            "r12": None if epoch else self.r12,
            "f107_from_r12": None if epoch else float(r12_to_f107(self.r12)),
            "iri_options": dict(self.iri_options),
            "pharlap_root": self.pharlap_root,
            "script_dir": self.script_dir,
            "workdir": self.workdir,
            "study": self.study,
            "tag": self.tag,
            "storm_model": ("per iri_options (PHaRLAP default: foF2 storm on)" if epoch
                            else "off, forced by user R12"),
            "availability": dict(self.availability),
        }


# ---------------------------------------------------------------------------
#: The IRI-2020 implementations the ``iri`` tier chooses between, by the name a
#: campaign uses, and the registered tier that answers for each.
IRI_MODELS: dict[str, str] = {"pharlap": "matlab_iri", "pyiri": "iri2020"}

#: The implementation the ``iri`` tier uses when a campaign does not name one.
IRI_DEFAULT_MODEL = "pharlap"


@register("l0", "iri")
class IRIBackend(Backend):
    """IRI-2020, with the implementation chosen by configuration.

    ``model`` names the code that answers::

        l0: {name: iri, f107: 120}                  # PHaRLAP, the default
        l0: {name: iri, model: pyiri, f107: 120}    # PyIRI

    ``pharlap`` (the default)
        The Fortran IRI-2020 inside PHaRLAP's ``iri2020`` mex, run through
        MATLAB by :class:`MatlabIRIBackend` (the ``matlab_iri`` tier).
    ``pyiri``
        PyIRI's Python reimplementation, :class:`PyIRIBackend` (the
        ``iri2020`` tier).

    PHaRLAP is the default because it is the reference code and because the
    two do not agree in the topside.  PyIRI 0.1.7 sets the NeQuick topside
    thickness with the NeQuick G law, H0 = k B2bot / v, where IRI-2020 uses
    H0 = k' B2bot with k' = (k e^{2(k-1)} + 1) / (e^{2(k-1)} + 1).  Its H0 is
    0.63 to 0.66 of IRI-2020's, and at 900 km and F10.7 = 120 its density is
    0.43 of IRI-2020's for that reason alone; with IRI-2020's law in its place
    the two agree to 0.6 to 2.5 per cent.  Against 5,118 ISIS-1 topside
    profiles IRI-2020's law is the closer of the two (median fitted H0 0.93 of
    it; PyIRI's 2.3 to 3.0 dB low).  ``docs/topside-model-campaign.md`` has the
    trace.  PyIRI stays available because it needs no MATLAB and is
    vectorised, which suits regression runs and quick looks; a topside or
    band-edge result should not be made on it.

    Solar activity is stated once, as ``f107`` (daily F10.7, SFU) or as
    ``r12``, not both.  Whichever the chosen code does not take is converted
    with IRI's own relation (:data:`~tss.l0_ionosphere.base.IRI_F107_COEFFS`),
    the one the PHaRLAP mex applies to a user R12, so ``f107: 120`` runs both
    implementations at F10.7 = 120 (R12 = 71.1 for PHaRLAP).  There is no
    default, for the reason :class:`PyIRIBackend` gives.

    ``epoch_indices=True`` instead runs PHaRLAP on each epoch's own indices
    (R12 = -1; :class:`MatlabIRIBackend` epoch mode), with neither ``f107`` nor
    ``r12`` stated.  PyIRI has no such mode and refuses it.

    ``topside`` is passed through: PHaRLAP gives ``NeQuick``, ``IRI-2001`` and
    ``IRI-2001 corrected``; PyIRI only ``NeQuick``.  Every other keyword goes to
    the chosen implementation unchanged, so ``pharlap_root``, ``timeout``,
    ``iri_options`` or ``probe`` reach :class:`MatlabIRIBackend` and
    ``foF2_coeff`` or ``hmF2_model`` reach :class:`PyIRIBackend`.  An option
    the chosen implementation does not take is an error, not silently dropped.

    The choice is recorded: ``label`` (and so ``IonosphereState.source``) and
    ``provenance`` name it, and :meth:`describe` carries ``model``, the
    implementing tier and that tier's own settings.  Attributes this class
    does not have, such as ``extra`` or ``last_info`` on the PHaRLAP route,
    are read from the implementation.
    """

    MODELS: dict[str, str] = IRI_MODELS
    DEFAULT_MODEL: str = IRI_DEFAULT_MODEL

    def __init__(
        self,
        model: str = IRI_DEFAULT_MODEL,
        topside: str = "NeQuick",
        f107: float | None = None,
        r12: float | None = None,
        *,
        epoch_indices: bool = False,
        geomagnetic_model: str = "auto",
        **options: Any,
    ) -> None:
        super().__init__(geomagnetic_model=geomagnetic_model)
        key = str(model).strip().lower()
        if key not in IRI_MODELS:
            raise ValueError(
                f"model must be one of {tuple(IRI_MODELS)}, got {model!r}"
            )
        if epoch_indices:
            if f107 is not None or r12 is not None:
                raise ValueError("epoch_indices=True takes the indices from the epoch; "
                                 "do not also state f107 or r12")
            if key != "pharlap":
                raise ValueError("epoch_indices needs model: pharlap; PyIRI takes only a stated F10.7")
            try:
                impl = MatlabIRIBackend(topside=topside, geomagnetic_model=geomagnetic_model,
                                        solar_mode="epoch", **options)
            except TypeError as exc:
                raise TypeError(
                    f"l0/iri model={key} passes its options to MatlabIRIBackend, "
                    f"which refused them: {exc}"
                ) from exc
            self.impl = impl
            self.model = key
            self.topside = topside
            self.f107 = None
            self.r12 = None
            self.solar_index_stated = "epoch"
            self.label = f"{impl.label} [l0/iri model={key}]"
            self.provenance = f"l0/iri model={key} -> l0/{IRI_MODELS[key]}; {impl.provenance}"
            return
        if (f107 is None) == (r12 is None):
            raise ValueError(
                "the iri tier needs one solar index, f107 (daily F10.7, SFU) or "
                "r12, and not both: e.g. l0: {name: iri, f107: 120}"
            )
        stated = "f107" if f107 is not None else "r12"
        if f107 is not None:
            f107_v = float(f107)
            r12_v = float(iri_r12_from_f107(f107_v))
        else:
            r12_v = float(r12)
            f107_v = float(iri_f107_from_r12(r12_v))

        if key == "pharlap":
            if not 1.0 <= r12_v <= 200.0:
                lo, hi = iri_f107_from_r12(1.0), iri_f107_from_r12(200.0)
                raise ValueError(
                    f"{stated} = {f107 if f107 is not None else r12} gives R12 = "
                    f"{r12_v:.1f}; the PHaRLAP route takes R12 in 1..200, which "
                    f"is F10.7 {lo:.1f}..{hi:.1f} SFU"
                )
            impl_cls: type[Backend] = MatlabIRIBackend
            solar: dict[str, float] = {"r12": r12_v}
        else:
            impl_cls = PyIRIBackend
            solar = {"f107": f107_v}

        try:
            impl = impl_cls(
                topside=topside, geomagnetic_model=geomagnetic_model, **solar, **options
            )
        except TypeError as exc:
            raise TypeError(
                f"l0/iri model={key} passes its options to {impl_cls.__name__}, "
                f"which refused them: {exc}"
            ) from exc

        self.impl = impl
        self.model = key
        self.topside = topside
        self.f107 = f107_v
        self.r12 = r12_v
        self.solar_index_stated = stated
        self.label = f"{impl.label} [l0/iri model={key}]"
        self.provenance = (
            f"l0/iri model={key} -> l0/{IRI_MODELS[key]}; stated {stated}, "
            f"F10.7={f107_v:.2f} R12={r12_v:.2f} by IRI's relation; {impl.provenance}"
        )

    def __getattr__(self, name: str) -> Any:
        # Only reached when ordinary lookup fails.  ``impl`` itself is excluded
        # so that a half-built or unpickled object raises instead of recursing.
        if name == "impl":
            raise AttributeError(name)
        return getattr(self.impl, name)

    def profiles(
        self,
        lat_deg: Sequence[float],
        lon_deg: Sequence[float],
        when: Sequence[dt.datetime],
        alt_grid_m: np.ndarray,
        **opts: Any,
    ) -> np.ndarray:
        """Electron density, m^-3, ``[n_samples, n_alt]``, from the chosen code."""
        return self.impl.profiles(lat_deg, lon_deg, when, alt_grid_m, **opts)

    def _extra_for_point(self, index: int) -> dict[str, float]:
        return self.impl._extra_for_point(index)

    def _settings(self) -> dict[str, Any]:
        return {
            "model": self.model,
            "implementation": f"l0/{IRI_MODELS[self.model]}",
            "implementation_class": type(self.impl).__name__,
            "topside": self.topside,
            "solar_index_stated": self.solar_index_stated,
            "f107": self.f107,
            "r12": self.r12,
            "implementation_settings": dict(self.impl._settings()),
        }


def _matlab_g(value: float) -> str:
    """Format a number the way MATLAB's ``%g`` would, for name and argument use."""
    text = f"{float(value):g}"
    return text


def make_backend(name: str, **kw: Any) -> Backend:
    """Build a backend by the short name the study used.

    Kept so that scripts written against ``iri_adapter.make_backend`` keep
    working.  New code should go through :func:`tss.registry.build`, which knows
    the layer as well as the name, and whose names are the registered tier names
    ``analytic``, ``iri``, ``iri2020``, ``iri_wrapper`` and ``matlab_iri``.
    """
    return {
        "analytic": AnalyticBackend,
        "pyiri": PyIRIBackend,
        "iriwrapper": IRIWrapperBackend,
        "matlab": MatlabIRIBackend,
    }[name](**kw)
