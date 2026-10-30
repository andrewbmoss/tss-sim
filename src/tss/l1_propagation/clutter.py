"""Layer: L1
WP: WP4
Claims: C1, C2, C4
Status: implemented
Requirements: INS-10, SCI-05, SCI-07

The clutter generator (WP4 Phase 1).

Until Phase 1 every ionogram tss-sim made held one clean O trace and one clean X
trace and nothing else, so a trace extractor tested on it was tested against an
ionosphere that does not exist.  Real topside ionograms (Alouette, ISIS) carry
other returns, and this module adds the ones that dominate them:

* **the ground echo** above the penetration frequency, where the wave goes
  through the F layer to the Earth and back (``clutter:ground``);
* **plasma resonances** at the spacecraft, the long ringing of the local plasma
  at f_N, f_UH and the harmonics of f_H that the transmitter excites, which
  appear as spikes across the whole range scale (``clutter:resonance:<name>``);
* **spread F**, oblique returns from field-aligned irregularities at ranges
  beyond the vertical echo, including a little above the penetration frequency
  (``clutter:spread_f``);
* **broadcast interference**, ground transmitters in the HF broadcast bands,
  which reach the spacecraft only above foF2 (rendered by L2, see below);
* **the front-end compression image** of WP3, which L2 already produces when
  the leakage drives the front end into compression.  It is listed so the
  extractor is tested against it; nothing is generated here.

The second-hop echo of the WP4 plan is left out on purpose (decision of
2026-10-30): a path of twice the group range needs a reflector at the
spacecraft, which is negligible, so it is not what clutters a topside ionogram.

HOW IT ENTERS THE PIPELINE
--------------------------
Echo-like terms (ground echo, spread F, resonances) are returned as extra rows of
the :class:`~tss.interfaces.EchoSet`, flagged ``clutter:<kind>``, on the
sounding-frequency grid the echo set already uses.  L2 synthesises them exactly
as it does an echo (stepped-sweep phase, step transition, Doppler), except that
a clutter row never joins a mode's dispersion trace.  Truth scoring
(:func:`tss.l5_assessment.scoring.truth_echo_cells`, and so INS-10 and every
trace measure) excludes every flagged row.  Carriers are not echoes: the
generator lists them on ``EchoSet.interference`` and L2 adds them at the
receiver input (:mod:`tss.l2_instrument.interference`).

The generator draws from its own random generator, which
:class:`tss.pipeline.Pipeline` derives from the seed and the configuration
separately from L2's, so turning clutter on does not by itself redraw the
receiver noise.

OCCURRENCE AND STRENGTH
-----------------------
By default every term is present at its declared strength in every run
(``climatology: false``), so the extractor is always tested.  The declared
values are representative, not measured; they are stated here and in
``config/default.yaml`` so a reader can change them:

* ground echo: surface reflection -3 dB (between sea, about -1 dB, and land,
  about -8 dB); non-deviative D-region loss A / (f + f_L)^2 two-way, with
  A = 100 dB MHz^2 and f_L = 1.2 MHz (George 1971's form at a modest level);
  the F-region loss is integrated through the truth profile
  (:func:`tss.l1_propagation.absorption.absorption_db`);
* spread F: total power 6 dB below the echo it accompanies, in 12 rows at
  offsets beyond the echo drawn from an exponential of mean 30 km truncated at
  100 km, and extending 5 per cent above each mode's penetration frequency;
* resonances: -70 dB (the EchoSet's path-gain units, per 3 km of range) at the
  start, decaying in power with time constant 0.5 ms (f_N), 1 ms (f_UH) and
  3 ms (n f_H, n = 1 to 3), to 30 dB down or 1000 km;
* interference: three carriers in each ITU HF broadcast band, on a 5 kHz
  raster, at -70 dBm (isotropic antenna) with a 6 dB log-normal spread, only
  above foF2.

``climatology: true`` makes the occurrence and two strengths depend on where and
when the run is (geographic latitude and local time, a stated simplification):

* spread F occurs with probability 0.7 at |lat| < 20 deg between 19 and 03 LT,
  0.5 poleward of 60 deg, 0.1 at night elsewhere, 0.02 otherwise;
* the D-region loss scales with cos(chi)^0.75 of the solar zenith angle at the
  spacecraft's sub-point (0.05 of the day value at night);
* the carriers are 10 dB weaker between 06 and 18 LT (D-region absorption of
  the ground signal).

RESONANCE ROWS
--------------
A resonance rings at its own frequency for milliseconds after the sweep passes
it.  It is modelled as a diffuse return confined to the sounding step that
contains the resonance frequency: rows every ``row_spacing_km`` of virtual range
with power falling as ``exp(-tau / tau_d)``.  The spectral fine structure inside
the step is not modelled.  ``level_db`` is the amplitude of the first row, so
the power per kilometre of range scales with the row spacing.
"""

from __future__ import annotations

import copy
import math
from collections.abc import Mapping
from dataclasses import replace
from typing import Any

import numpy as np

from tss.interfaces import EchoSet, IonosphereState, PlatformState

__all__ = [
    "CLUTTER_PREFIX",
    "DEFAULTS",
    "BROADCAST_BANDS_MHZ",
    "ClutterGenerator",
    "is_clutter_flag",
    "local_time_h",
    "d_region_loss_db",
]

C_LIGHT = 299_792_458.0
CLUTTER_PREFIX = "clutter:"
_NO_ECHO = ("evanescent", "penetrating")
_DB_PER_NEPER_POWER = 10.0 / math.log(10.0)  # 4.343 dB per e-fold of power

#: ITU HF broadcast bands (Radio Regulations Art. 5), MHz.
BROADCAST_BANDS_MHZ: tuple[tuple[float, float], ...] = (
    (3.90, 4.00), (4.75, 5.06), (5.90, 6.20), (7.20, 7.45), (9.40, 9.90),
    (11.60, 12.10), (13.57, 13.87), (15.10, 15.80), (17.48, 17.90), (18.90, 19.02),
)

DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "climatology": False,
    "ground_echo": {
        "enabled": True,
        "surface_reflection_db": -3.0,
        "d_region_a_db_mhz2": 100.0,
        "d_region_f_l_mhz": 1.2,
    },
    "spread_f": {
        "enabled": True,
        "level_db": -6.0,
        "rows_per_echo": 12,
        "mean_km": 30.0,
        "extent_km": 100.0,
        "beyond_penetration_fraction": 0.05,
    },
    "plasma_resonances": {
        "enabled": True,
        "row_spacing_km": 3.0,
        "extent_km": 1000.0,
        "floor_db": 30.0,
        "fN": {"enabled": True, "level_db": -70.0, "decay_ms": 0.5},
        "fUH": {"enabled": True, "level_db": -70.0, "decay_ms": 1.0},
        "nfH": {"enabled": True, "level_db": -70.0, "decay_ms": 3.0, "harmonics": 3},
    },
    "interference": {
        "enabled": True,
        "bands_mhz": [list(b) for b in BROADCAST_BANDS_MHZ],
        "carriers_per_band": 3,
        "power_dbm_iso": -70.0,
        "sigma_db": 6.0,
        "raster_hz": 5000.0,
        "only_above_fof2": True,
    },
    "compression_image": {"enabled": True},
}


def is_clutter_flag(flag: Any) -> bool:
    """Whether an EchoSet row flag marks a clutter row."""
    return str(flag).startswith(CLUTTER_PREFIX)


def _merge(base: Mapping[str, Any], over: Mapping[str, Any] | None) -> dict[str, Any]:
    out = copy.deepcopy(dict(base))
    for k, v in (over or {}).items():
        if isinstance(v, Mapping) and isinstance(out.get(k), Mapping):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def local_time_h(lon_deg: float, epoch: str) -> float:
    """Local solar time (h) at a longitude, from the epoch's UT."""
    ut = float(epoch[11:13]) + float(epoch[14:16]) / 60.0 if len(epoch) >= 16 else 12.0
    return float((ut + lon_deg / 15.0) % 24.0)


def _day_of_year(epoch: str) -> int:
    from datetime import date

    try:
        d = date(int(epoch[0:4]), int(epoch[5:7]), int(epoch[8:10]))
    except (ValueError, IndexError):
        return 1
    return d.timetuple().tm_yday


def cos_solar_zenith(lat_deg: float, lon_deg: float, epoch: str) -> float:
    """cos(chi) at the sub-point, from declination and hour angle (no equation of time)."""
    decl = math.radians(-23.44 * math.cos(2.0 * math.pi * (_day_of_year(epoch) + 10) / 365.0))
    hour = math.radians(15.0 * (local_time_h(lon_deg, epoch) - 12.0))
    lat = math.radians(lat_deg)
    return math.sin(lat) * math.sin(decl) + math.cos(lat) * math.cos(decl) * math.cos(hour)


def d_region_loss_db(f_hz: float, a_db_mhz2: float, f_l_mhz: float, scale: float = 1.0) -> float:
    """Two-way non-deviative D-region loss, ``scale * A / (f + f_L)^2`` dB (f in MHz)."""
    return float(scale * a_db_mhz2 / (f_hz / 1.0e6 + f_l_mhz) ** 2)


def _spread_f_probability(lat_deg: float, lt_h: float) -> float:
    night = lt_h >= 19.0 or lt_h < 3.0
    if abs(lat_deg) < 20.0:
        return 0.7 if night else 0.02
    if abs(lat_deg) >= 60.0:
        return 0.5
    return 0.1 if (lt_h >= 18.0 or lt_h < 6.0) else 0.02


class ClutterGenerator:
    """Adds the topside clutter terms to an echo set (module docstring)."""

    def __init__(self, settings: Mapping[str, Any] | None = None) -> None:
        self.settings = _merge(DEFAULTS, settings)
        #: What the last call to :meth:`augment` added, for the run record.
        self.last_report: dict[str, Any] = {}

    @classmethod
    def from_config(cls, block: Mapping[str, Any] | None) -> ClutterGenerator | None:
        """A generator from a ``clutter`` config block, or ``None`` when it is off.

        A missing or empty block is off: a campaign gets clutter from
        ``config/default.yaml``, which declares it on, and a mapping assembled
        by hand without it gets none.
        """
        if not block or not bool(block.get("enabled", False)):
            return None
        return cls(block)

    # -- contract ------------------------------------------------------------

    def describe(self) -> dict[str, Any]:
        s = self.settings
        return {
            "generator": "tss.l1_propagation.clutter",
            "climatology": bool(s["climatology"]),
            "terms": {k: (dict(s[k]) if isinstance(s[k], Mapping) else s[k]) for k in
                      ("ground_echo", "spread_f", "plasma_resonances", "interference", "compression_image")},
            "compression_image_source": (
                "L2 front end (tss.l2_instrument.receive_chain): present when the leakage drives it into "
                "compression, e.g. with the canceller off; not generated by the clutter generator"),
        }

    def augment(self, echoes: EchoSet, state: IonosphereState, platform: PlatformState,
                freq_hz: np.ndarray, rng: np.random.Generator) -> EchoSet:
        """The echo set with clutter rows appended and carriers listed.  Inputs are not modified."""
        from tss.l1_propagation.base import critical_frequency_hz, plasma_frequency_at

        s = self.settings
        grid = np.unique(np.asarray(freq_hz, dtype=float))
        h_sat = float(platform.alt_m)
        epoch = str(platform.epoch or state.epoch or "")
        lat, lon = float(platform.lat_deg), float(platform.lon_deg)
        lt = local_time_h(lon, epoch) if epoch else 12.0
        fo = float(critical_frequency_hz(state))
        f_n = float(plasma_frequency_at(state, h_sat))
        report: dict[str, Any] = {"fof2_hz": fo, "f_n_sat_hz": f_n, "local_time_h": lt,
                                  "climatology": bool(s["climatology"])}
        parts: list[dict[str, np.ndarray]] = []

        if s["ground_echo"].get("enabled", True):
            scale = 1.0
            if s["climatology"]:
                c = cos_solar_zenith(lat, lon, epoch)
                scale = max(c, 0.0) ** 0.75 if c > 0.0 else 0.05
                report["d_region_scale"] = scale
            p = self._ground_echo(echoes, state, h_sat, scale)
            report["ground_echo_rows"] = int(p["freq"].size)
            parts.append(p)

        if s["spread_f"].get("enabled", True):
            occurs = True
            if s["climatology"]:
                prob = _spread_f_probability(lat, lt)
                occurs = bool(rng.random() < prob)
                report["spread_f_probability"] = prob
            report["spread_f_occurs"] = occurs
            p = self._spread_f(echoes, rng) if occurs else _empty()
            report["spread_f_rows"] = int(p["freq"].size)
            parts.append(p)

        if s["plasma_resonances"].get("enabled", True):
            p, res = self._resonances(state, h_sat, f_n, grid)
            report["resonances"] = res
            report["resonance_rows"] = int(p["freq"].size)
            parts.append(p)

        carriers = None
        if s["interference"].get("enabled", True) and grid.size:
            day = 6.0 <= lt < 18.0
            carriers = self._carriers(fo, grid, rng, -10.0 if (s["climatology"] and day) else 0.0)
            report["carriers"] = int(carriers["freq_hz"].size)

        report["compression_image"] = "listed; produced by L2 when the front end compresses"
        self.last_report = report
        return _append(echoes, parts, carriers)

    # -- terms ---------------------------------------------------------------

    def _ground_echo(self, echoes: EchoSet, state: IonosphereState, h_sat: float, scale: float
                     ) -> dict[str, np.ndarray]:
        from tss.l1_propagation import magnetoionic
        from tss.l1_propagation.absorption import absorption_db
        from tss.l1_propagation.analytic import group_range_m as ff_group_range_m

        g = self.settings["ground_echo"]
        flags = _flags(echoes)
        mode = _modes(echoes)
        freq = np.asarray(echoes.freq_hz, dtype=float)
        dopp = _col(echoes.doppler_hz, freq.size)
        h_b = float(np.min(state.alt_m))
        b = state.b_field_t
        field_aware = b is not None and np.ndim(b) == 2 and np.shape(b)[0] == 3
        out = _rows()
        for i in np.nonzero(np.array([f == "penetrating" for f in flags], dtype=bool))[0]:
            f, m = float(freq[i]), str(mode[i])
            try:
                if field_aware:
                    gr = magnetoionic.group_range_m(state, f, h_b, h_sat, m)
                elif m == "O":
                    gr = ff_group_range_m(state, f, h_b, h_sat)
                else:
                    continue
            except (ValueError, FloatingPointError):
                continue
            if not np.isfinite(gr):
                continue
            p_group = float(gr) + h_b  # vacuum below the bottom of the grid
            a_iono = float(absorption_db(state, f, h_sat, h_refl_m=h_b))
            a_d = d_region_loss_db(f, float(g["d_region_a_db_mhz2"]), float(g["d_region_f_l_mhz"]), scale)
            amp = (20.0 * math.log10((C_LIGHT / f) / (4.0 * math.pi * 2.0 * p_group))
                   - a_iono - a_d + float(g["surface_reflection_db"]))
            _add(out, f, m, 2.0 * p_group / C_LIGHT, amp, a_iono + a_d, dopp[i], 0.0, "clutter:ground")
        return _arr(out)

    def _spread_f(self, echoes: EchoSet, rng: np.random.Generator) -> dict[str, np.ndarray]:
        sp = self.settings["spread_f"]
        n = max(int(sp["rows_per_echo"]), 1)
        mean_m = float(sp["mean_km"]) * 1e3
        ext_m = float(sp["extent_km"]) * 1e3
        each_db = float(sp["level_db"]) - 10.0 * math.log10(n)
        freq = np.asarray(echoes.freq_hz, dtype=float)
        delay = np.asarray(echoes.group_delay_s, dtype=float)
        amp = np.asarray(echoes.amplitude_db, dtype=float)
        absn = _col(echoes.absorption_db, freq.size)
        dopp = _col(echoes.doppler_hz, freq.size)
        flags = _flags(echoes)
        mode = _modes(echoes)
        parent = (np.isfinite(freq) & np.isfinite(delay) & (delay > 0) & np.isfinite(amp)
                  & np.array([f not in _NO_ECHO and not is_clutter_flag(f) for f in flags], dtype=bool))
        out = _rows()

        def offsets() -> np.ndarray:
            # Exponential of mean ``mean_m`` truncated at ``ext_m``, by inverse CDF.
            u = rng.random(n)
            return -mean_m * np.log1p(-u * (1.0 - math.exp(-ext_m / mean_m)))

        for i in np.nonzero(parent)[0]:
            for d in offsets():
                _add(out, freq[i], mode[i], delay[i] + 2.0 * d / C_LIGHT, amp[i] + each_db, absn[i],
                     dopp[i], np.nan, "clutter:spread_f")
        frac = float(sp.get("beyond_penetration_fraction", 0.0) or 0.0)
        if frac > 0.0:
            for m in sorted(set(mode[parent].tolist())):
                mine = parent & (mode == m)
                top = int(np.nonzero(mine)[0][np.argmax(freq[mine])])
                pen = np.nonzero(np.array([fl == "penetrating" for fl in flags], dtype=bool) & (mode == m)
                                 & (freq > freq[top]) & (freq <= freq[top] * (1.0 + frac)))[0]
                for i in pen:
                    for d in offsets():
                        _add(out, freq[i], m, delay[top] + 2.0 * d / C_LIGHT, amp[top] + each_db, absn[top],
                             dopp[i], np.nan, "clutter:spread_f")
        return _arr(out)

    def _resonances(self, state: IonosphereState, h_sat: float, f_n: float, grid: np.ndarray
                    ) -> tuple[dict[str, np.ndarray], list[dict[str, Any]]]:
        from tss.l1_propagation import magnetoionic

        r = self.settings["plasma_resonances"]
        try:
            f_h, _theta = magnetoionic.field_at(state, h_sat)
            f_h = float(np.asarray(f_h).ravel()[0])
        except (ValueError, TypeError):
            f_h = float("nan")
        lines: list[tuple[str, float, Mapping[str, Any]]] = []
        if r["fN"].get("enabled", True):
            lines.append(("fN", f_n, r["fN"]))
        if r["fUH"].get("enabled", True) and np.isfinite(f_h):
            lines.append(("fUH", math.hypot(f_n, f_h), r["fUH"]))
        if r["nfH"].get("enabled", True) and np.isfinite(f_h):
            for k in range(1, int(r["nfH"].get("harmonics", 3)) + 1):
                lines.append((f"{k}fH", k * f_h, r["nfH"]))
        dr = float(r["row_spacing_km"]) * 1e3
        out = _rows()
        record: list[dict[str, Any]] = []
        for name, f_r, spec in lines:
            j = _step_for(f_r, grid)
            tau_d = float(spec["decay_ms"]) * 1e-3
            n_floor = float(r["floor_db"]) * tau_d / _DB_PER_NEPER_POWER  # delay at which it is floor_db down
            n_rows = int(min(float(r["extent_km"]) * 1e3 / dr, n_floor * C_LIGHT / (2.0 * dr)))
            record.append({"name": name, "f_hz": f_r, "grid_hz": None if j is None else float(grid[j]),
                           "rows": 0 if j is None else n_rows})
            if j is None or n_rows < 1:
                continue
            for q in range(1, n_rows + 1):
                tau = 2.0 * q * dr / C_LIGHT
                amp = float(spec["level_db"]) - _DB_PER_NEPER_POWER * tau / tau_d
                _add(out, grid[j], "-", tau, amp, 0.0, np.nan, np.nan, f"clutter:resonance:{name}")
        return _arr(out), record

    def _carriers(self, fo: float, grid: np.ndarray, rng: np.random.Generator, shift_db: float
                  ) -> dict[str, np.ndarray]:
        it = self.settings["interference"]
        lo, hi = float(grid.min()), float(grid.max())
        f_list, p_list = [], []
        raster = float(it["raster_hz"])
        for b_lo, b_hi in it["bands_mhz"]:
            n_slots = int((float(b_hi) - float(b_lo)) * 1e6 // raster)
            if n_slots < 1:
                continue
            k = int(it["carriers_per_band"])
            slots = rng.choice(n_slots + 1, size=min(k, n_slots + 1), replace=False)
            powers = float(it["power_dbm_iso"]) + shift_db + float(it["sigma_db"]) * rng.standard_normal(slots.size)
            for sl, pw in zip(slots, powers, strict=True):
                f = float(b_lo) * 1e6 + float(sl) * raster
                if not lo <= f <= hi:
                    continue
                if bool(it.get("only_above_fof2", True)) and not f > fo:
                    continue
                f_list.append(f)
                p_list.append(float(pw))
        n = len(f_list)
        return {"freq_hz": np.asarray(f_list, dtype=float), "power_dbm_iso": np.asarray(p_list, dtype=float),
                "phase_rad": rng.uniform(0.0, 2.0 * math.pi, n)}


# -- helpers -------------------------------------------------------------------


def _step_for(f: float, grid: np.ndarray) -> int | None:
    """Index of the grid frequency (step centre) whose step contains ``f``, or ``None``."""
    if not np.isfinite(f) or grid.size == 0:
        return None
    j = int(np.argmin(np.abs(grid - f)))
    if grid.size == 1:
        return j
    lo = grid[j] - 0.5 * (grid[j] - grid[j - 1]) if j > 0 else grid[0] - 0.5 * (grid[1] - grid[0])
    hi = grid[j] + 0.5 * (grid[j + 1] - grid[j]) if j < grid.size - 1 else grid[-1] + 0.5 * (grid[-1] - grid[-2])
    return j if lo <= f <= hi else None


def _flags(e: EchoSet) -> np.ndarray:
    n = np.size(e.freq_hz)
    if e.flags is None:
        return np.array([""] * n, dtype=object)
    return np.array([str(x) for x in np.asarray(e.flags, dtype=object).ravel()[:n]], dtype=object)


def _modes(e: EchoSet) -> np.ndarray:
    n = np.size(e.freq_hz)
    if e.mode is None or np.size(e.mode) != n:
        return np.array([""] * n, dtype=object)
    return np.array([str(x) for x in np.asarray(e.mode, dtype=object).ravel()], dtype=object)


def _col(x: Any, n: int) -> np.ndarray:
    if x is None or np.size(x) != n:
        return np.full(n, np.nan)
    return np.asarray(x, dtype=float)


_KEYS = ("freq", "mode", "delay", "amp", "abs", "dopp", "alt", "flag")


def _rows() -> dict[str, list[Any]]:
    return {k: [] for k in _KEYS}


def _add(out: dict[str, list[Any]], f: float, m: str, tau: float, amp: float, absn: float, dopp: float,
         alt: float, flag: str) -> None:
    for k, v in zip(_KEYS, (f, m, tau, amp, absn, dopp, alt, flag), strict=True):
        out[k].append(v)


def _arr(out: dict[str, list[Any]]) -> dict[str, np.ndarray]:
    return {k: np.asarray(v, dtype=object if k in ("mode", "flag") else float) for k, v in out.items()}


def _empty() -> dict[str, np.ndarray]:
    return _arr(_rows())


def _append(e: EchoSet, parts: list[dict[str, np.ndarray]], carriers: dict[str, np.ndarray] | None) -> EchoSet:
    n = np.size(e.freq_hz)
    new = {k: np.concatenate([p[k] for p in parts]) if parts else np.array([], dtype=float) for k in _KEYS}
    add = int(np.size(new["freq"]))
    src = e.source + (f" + clutter ({add} rows)" if add else " + clutter")
    if add == 0:
        return replace(e, interference=carriers, source=src)
    alt = _col(e.reflection_alt_m, n) if e.reflection_alt_m is not None else None
    return replace(
        e,
        freq_hz=np.concatenate([np.asarray(e.freq_hz, dtype=float), new["freq"].astype(float)]),
        mode=np.concatenate([_modes(e), new["mode"].astype(object)]),
        group_delay_s=np.concatenate([np.asarray(e.group_delay_s, dtype=float), new["delay"].astype(float)]),
        amplitude_db=np.concatenate([np.asarray(e.amplitude_db, dtype=float), new["amp"].astype(float)]),
        absorption_db=np.concatenate([_col(e.absorption_db, n), new["abs"].astype(float)]),
        doppler_hz=(None if e.doppler_hz is None
                    else np.concatenate([_col(e.doppler_hz, n), new["dopp"].astype(float)])),
        aoa_deg=(None if e.aoa_deg is None
                 else np.concatenate([_col(e.aoa_deg, n), np.full(add, np.nan)])),
        reflection_alt_m=(None if alt is None else np.concatenate([alt, new["alt"].astype(float)])),
        flags=np.concatenate([_flags(e), new["flag"].astype(object)]),
        interference=carriers,
        source=src,
    )
