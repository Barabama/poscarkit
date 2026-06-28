"""Thermo pipeline — read data + TDB → compute Sconf + DeltaG → save results + plot."""

import logging
import os
from pathlib import Path

import numpy as np

from src.io.readers import read as read_ir
from src.io.ir import SOFData, get_site_ratios
from src.thermo.tdb import parse_tdb
from src.thermo.sconf import calc_sconf, calc_sconf_random
from src.thermo.deltaG import calc_deltaG_random
from src.thermo.plot import save_plot


def run_thermo(
    data_path: str,
    tdb_path: str,
    name: str = "thermo",
    outdir: str = "output",
) -> tuple[Path, Path]:
    """Read SOF data and TDB, compute Sconf + DeltaG, save results.

    Args:
        data_path: Path to CSV/XLSX data file.
        tdb_path: Path to TDB thermodynamic database.
        name: Work name for output subdirectory.
        outdir: Output root directory.

    Returns:
        (csv_path, png_path) — output file paths.
    """
    if not os.path.isfile(data_path):
        raise FileNotFoundError(f"Data file not found: {data_path}")
    if not os.path.isfile(tdb_path):
        raise FileNotFoundError(f"TDB file not found: {tdb_path}")

    logging.info(f"Parsing TDB: {tdb_path}")
    tdb_data = parse_tdb(tdb_path)
    logging.info(
        f"  Functions: {len(tdb_data['functions'])}, " f"Parameters: {len(tdb_data['params'])}"
    )
    for p, r in tdb_data.get("phase_ratios", {}).items():
        logging.info(f"  Phase {p}: {r} ({len(r)} sublattices)")

    logging.info(f"Reading data: {data_path}")
    ir = read_ir(data_path)

    tdb_ratios = tdb_data.get("phase_ratios", {})
    if tdb_ratios:
        ir.site_ratios = get_site_ratios(ir.phase, tdb_ratios)

    logging.info(f"  Phase: {ir.phase}, site ratios: {ir.site_ratios}")
    logging.info(f"  Elements: {ir.elements}")
    logging.info(f"  T range: {ir.T[0]:.0f} ~ {ir.T[-1]:.0f} K ({ir.n_points} points)")

    # Sconf_real(T) — configurational entropy from actual site fractions
    Sconf_real = np.zeros(ir.n_points)
    for i in range(ir.n_points):
        sofs_at_T = [
            {e: float(ir.Y_subl[k][e][i]) for e in ir.elements if e in ir.Y_subl[k]}
            for k in range(ir.n_subl)
        ]
        for sofs in sofs_at_T:
            _normalize(sofs)
        Sconf_real[i] = calc_sconf(sofs_at_T, ir.site_ratios)

    logging.info(f"  Sconf_real range: {Sconf_real.min():.4f} ~ {Sconf_real.max():.4f} J/(mol*K)")

    # Sconf_random — configurational entropy for random mixing (scalar)
    Sconf_random = calc_sconf_random(ir.composition, ir.site_ratios)
    logging.info(f"  Sconf_random = {Sconf_random:.4f} J/(mol*K)")

    # DeltaG_random(T)
    G_random = calc_deltaG_random(
        ir.composition,
        tdb_data["functions"],
        tdb_data["params"],
        ir.T,
        Sconf_random,
        ir.phase,
        ir.site_ratios,
    )
    logging.info(f"  DeltaG_random range: {G_random[0]:.2f} ~ {G_random[-1]:.2f} J/mol")

    # S_random(T), H_random(T) via numerical differentiation of G_random
    dG_dT = np.gradient(G_random, ir.T)
    S_random = -dG_dT  # total entropy (includes configurational + vibrational + magnetic)
    H_random = G_random + ir.T * S_random

    outdir_p = Path(outdir) / name
    outdir_p.mkdir(parents=True, exist_ok=True)

    csv_path = _save_csv(ir, G_random, Sconf_real, Sconf_random, H_random, S_random, outdir_p)
    logging.info(f"  -> data: {csv_path}")

    png_path = save_plot(ir, G_random, Sconf_real, Sconf_random, H_random, S_random, outdir_p)
    logging.info(f"  -> plot: {png_path}")

    return csv_path, png_path


def _normalize(sofs: dict[str, float]) -> None:
    """Normalize site fractions in place to sum to 1."""
    total = sum(sofs.values())
    if total > 1e-15:
        for e in sofs:
            sofs[e] /= total


def _save_csv(
    ir: SOFData,
    G_random: np.ndarray,
    Sconf_real: np.ndarray,
    Sconf_random: float,
    H_random: np.ndarray,
    S_random: np.ndarray,
    outdir: Path,
) -> Path:
    """Save results as CSV."""
    import pandas as pd

    nan_arr = np.full(ir.n_points, np.nan)
    comp_str = _composition_str(ir.composition)
    df = pd.DataFrame(
        {
            "T (K)": ir.T,
            "DeltaG_f_real (J/mol-atom)": ir.G_real if ir.G_real is not None else nan_arr,
            "DeltaG_f_random (J/mol-atom)": G_random,
            "DeltaH_f_real (J/mol-atom)": ir.H_real if ir.H_real is not None else nan_arr,
            "DeltaH_f_random (J/mol-atom)": H_random,
            "DeltaS_f_real (J/mol-atom-K)": ir.S_real if ir.S_real is not None else nan_arr,
            "DeltaS_f_random (J/mol-atom-K)": S_random,
            "S_conf_real (J/mol-atom-K)": Sconf_real,
            "S_conf_random (J/mol-atom-K)": np.full(ir.n_points, Sconf_random),
        }
    )
    fpath = outdir / f"Thermo_{comp_str}.csv"
    df.to_csv(fpath, index=False)
    return fpath


def _composition_str(comp: dict[str, float]) -> str:
    parts = []
    for elem in sorted(comp.keys()):
        v = comp[elem]
        parts.append(f"{elem}{v:.3g}".replace(".", "_"))
    return "-".join(parts)
