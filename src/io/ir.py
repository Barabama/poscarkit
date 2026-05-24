"""Intermediate Representation — unified data model for imported SOF data."""

import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

# Fallback phase name → sublattice site ratios (must sum to 1.0).
# When a TDB is provided, ratios are extracted from the TDB's PHASE records
# and take precedence over this map via get_site_ratios().
PHASE_SITE_RATIOS: dict[str, list[float]] = {
    "FCC": [0.25, 0.75],
    "BCC": [0.5, 0.5],
    "HCP": [0.25, 0.75],
}

PHASE_SUBLATTICE_MAP: dict[str, dict[int, str]] = {
    "FCC": {1: "1a", 2: "3c"},
    "BCC": {1: "1a", 2: "1b"},
    "HCP": {1: "2a", 2: "6c"},
}

_RE_MAP_SPEC = re.compile(r"(\d+):(\w+)")


def parse_sublattice_map(spec: str) -> dict[int, str]:
    """Parse a sublattice-map string like '1:1a,2:3c' into {1: '1a', 2: '3c'}.

    Tolerates common input variations:
    Chinese punctuation (：，) to ASCII (: ,), extra whitespace.
    """
    spec = spec.replace("：", ":").replace("，", ",")
    spec = spec.replace("．", ".").replace("、", ",")
    spec = re.sub(r"\s+", "", spec)

    result: dict[int, str] = {}
    for part in spec.split(","):
        if not part:
            continue
        m = _RE_MAP_SPEC.search(part)
        if not m:
            raise ValueError(f"Invalid sublattice mapping '{part}'. Expected 'N:name'")
        result[int(m.group(1))] = m.group(2)
    return result


def get_site_ratios(
    phase: str,
    tdb_phase_ratios: dict[str, list[float]] | None = None,
) -> list[float]:
    """Resolve site ratios for a phase, preferring TDB data over fallback map."""
    phase_up = phase.upper()
    if tdb_phase_ratios and phase_up in tdb_phase_ratios:
        return tdb_phase_ratios[phase_up]
    if phase_up in PHASE_SITE_RATIOS:
        return PHASE_SITE_RATIOS[phase_up]
    raise KeyError(
        f"Unknown phase '{phase}'. Add it to PHASE_SITE_RATIOS "
        f"or ensure the TDB file contains a PHASE definition for it."
    )


@dataclass
class SOFData:
    """Unified intermediate representation emitted by all format readers."""

    source_path: str
    phase: str
    site_ratios: list[float]
    T: np.ndarray
    Y_subl: list[dict[str, np.ndarray]]
    composition: dict[str, float]
    elements: list[str]
    G_real: np.ndarray | None = None

    def __post_init__(self):
        if not self.site_ratios:
            raise ValueError("site_ratios must not be empty")
        if abs(sum(self.site_ratios) - 1.0) > 1e-6:
            raise ValueError(f"site_ratios must sum to 1.0, got {self.site_ratios}")
        if len(self.Y_subl) != len(self.site_ratios):
            raise ValueError(
                f"Y_subl has {len(self.Y_subl)} sublattices but "
                f"site_ratios has {len(self.site_ratios)}"
            )

    @property
    def n_subl(self) -> int:
        return len(self.site_ratios)

    @property
    def n_points(self) -> int:
        return len(self.T)

    def dump_csv(self, path: str | Path) -> Path:
        """Export SOFData as a flat CSV for inspection or downstream use."""
        path = Path(path)
        rows = []
        for i in range(self.n_points):
            row = {"T": self.T[i]}
            for k, y_dict in enumerate(self.Y_subl):
                for elem in self.elements:
                    row[f"Y_subl{k+1}_{elem}"] = y_dict.get(elem, np.zeros(1))[i]
            if self.G_real is not None:
                row["G_real"] = self.G_real[i]
            rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(path, index=False)
        return path


def get_sofs_at(
    ir: SOFData, T: float, sublattice_map: dict[int, str] | None = None
) -> dict[str, dict[str, float]]:
    """Extract SOFs at a single temperature, mapped to Wyckoff site names.

    Args:
        ir: SOFData from a reader.
        T: Target temperature.
        sublattice_map: CALPHAD sub-lattice # to Wyckoff name. Defaults to
            PHASE_SUBLATTICE_MAP[phase].

    Returns:
        {wyckoff_site: {element: fraction, ...}, ...}
    """
    if sublattice_map is None:
        sublattice_map = PHASE_SUBLATTICE_MAP.get(ir.phase.upper(), {})

    available = ir.T
    if T in available:
        idx = np.where(available == T)[0][0]
    else:
        nearest = min(available, key=lambda x: abs(x - T))
        idx = np.where(available == nearest)[0][0]
        logging.warning(
            f"Temperature {T} not found in data. Using nearest: {nearest}"
        )

    result: dict[str, dict[str, float]] = {}
    for calphad_idx, y_dict in enumerate(ir.Y_subl, start=1):
        site_name = sublattice_map.get(calphad_idx, f"subl{calphad_idx}")
        sofs: dict[str, float] = {}
        for elem in ir.elements:
            if elem in y_dict:
                val = float(y_dict[elem][idx])
                if val > 1e-12:
                    sofs[elem] = val
        total = sum(sofs.values())
        if total > 0 and abs(total - 1.0) > 1e-6:
            logging.warning(
                f"SOFs for site {site_name} sum to {total:.6f}, normalizing."
            )
            sofs = {e: f / total for e, f in sofs.items()}
        elif total == 0:
            logging.warning(f"No SOFs found for site {site_name} at T={T}")
            continue
        result[site_name] = sofs

    return result


def build_structure_info(
    config_cfg: dict, sofs_by_site: dict[str, dict[str, float]]
) -> dict:
    """Merge config geometry with SOF data into structure_info dict.

    Args:
        config_cfg: Full config dict (from config.toml).
        sofs_by_site: {wyckoff_site: {element: fraction, ...}} from get_sofs_at.

    Returns:
        structure_info dict ready for run_modeling().
    """
    result = {"cell": config_cfg.get("cell", [])}
    for site_name, sofs in sofs_by_site.items():
        site_data = config_cfg.get(site_name, {})
        result[site_name] = {
            "atoms": site_data.get("atoms", []),
            "sofs": sofs,
        }
    return result
