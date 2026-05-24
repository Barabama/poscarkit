"""Import SOF data from CSV/XLSX and run modeling — shared by CLI, interactive, GUI."""

import logging
import tomllib
import json
from pathlib import Path

from src.io.readers import detect_format, read as read_sof_data
from src.io.ir import get_sofs_at, build_structure_info, parse_sublattice_map
from src.workflow.modeling import run_modeling
from src.config import normalize_config_keys


def run_import_to_model(
    csv_path: str,
    config_path: str,
    phase: str,
    temperatures: list[float],
    sublattice_map: str | None = None,
    outdir: str = "output",
    name: str = "modeling",
    factors: tuple[int, int, int] = (3, 3, 3),
    seeds: list[int | None] | None = None,
    batch_size: int = 1,
    enable_sqs: bool = False,
    iterations: int = int(1e7),
    output_mode: str = "run",
) -> list[Path]:
    """Read CSV/XLSX site-fraction data and run modeling for each temperature.

    Args:
        csv_path: Path to CSV or XLSX file.
        config_path: Path to config.toml (geometry source).
        phase: Phase name (FCC, BCC, HCP).
        temperatures: List of temperatures to process.
        sublattice_map: Optional CALPHAD to Wyckoff mapping override.
        outdir: Output root directory.
        name: Work name prefix.
        factors: Supercell factors.
        seeds: Shuffle seeds (None = no seed).
        batch_size: Modeling batch size.
        enable_sqs: Enable SQS generation.
        iterations: SQS iterations.
        output_mode: "run" | "save-config" | "print".

    Returns:
        List of output paths (empty for print mode).
    """
    fmt = detect_format(csv_path)
    logging.info(f"Detected format: {fmt}")
    ir = read_sof_data(csv_path, phase_hint=phase)

    smap = parse_sublattice_map(sublattice_map) if sublattice_map else None

    with open(config_path, "rb") as f:
        cfg = normalize_config_keys(tomllib.load(f))
    phase_cfg = cfg.get(phase.upper(), {})

    results: list[Path] = []
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    for T in temperatures:
        sofs_by_site = get_sofs_at(ir, T, smap)
        logging.info(f"T={T} K: {sofs_by_site}")

        structure_info = build_structure_info(phase_cfg, sofs_by_site)

        if output_mode == "print":
            print(json.dumps({"T": T, "sofs": sofs_by_site}, indent=2))
            continue

        if output_mode == "save-config":
            out_config = outdir_p / f"config_{int(T)}K.toml"
            _write_config_snippet(config_path, phase, sofs_by_site, out_config)
            logging.info(f"Config saved to {out_config}")
            results.append(out_config)
            continue

        # output_mode == "run"
        name_t = f"{name}_{int(T)}K"
        poscar_path = outdir_p / name_t
        poscar_path.mkdir(parents=True, exist_ok=True)
        # Generate unitcell from structure_info for run_modeling
        from src.modeling.supercell import unitcell2file
        uc = unitcell2file(structure_info=structure_info, outdir=poscar_path)

        run_results = run_modeling(
            name=name_t,
            poscar=uc,
            outdir=outdir_p,
            supercell_factors=factors,
            structure_info=structure_info,
            shuffle_seeds=seeds or [None],
            batch_size=batch_size,
            enable_sqs=enable_sqs,
            iterations=iterations,
        )
        results.extend(run_results)

    return results


def _write_config_snippet(
    config_path: str,
    phase: str,
    sofs_by_site: dict[str, dict[str, float]],
    out_path: Path,
):
    """Write a modified config.toml with updated SOF sections."""
    import re
    phase_upper = phase.upper()
    with open(config_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    wanted = {f"{phase_upper}.{site}.sofs": sofs
              for site, sofs in sofs_by_site.items()}
    new_lines = []
    skip_section = False
    seen: set[str] = set()
    for line in lines:
        m = re.match(r"^\s*\[([^\]]+)\]\s*$", line)
        if m:
            section = m.group(1).strip()
            skip_section = section in wanted
            if skip_section:
                new_lines.append(line)
                for elem, val in sorted(wanted[section].items()):
                    new_lines.append(f"{elem} = {val:.6g}\n")
                new_lines.append("\n")
                seen.add(section)
                continue
        if skip_section and line.strip() and not line.strip().startswith("#"):
            continue
        new_lines.append(line)
    for section, sofs in wanted.items():
        if section not in seen:
            new_lines.append(f"\n[{section}]\n")
            for elem, val in sorted(sofs.items()):
                new_lines.append(f"{elem} = {val:.6g}\n")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.writelines(new_lines)
