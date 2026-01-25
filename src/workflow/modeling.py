# src/workflow/modeling.py

# import os
import logging
from typing import Any
from pathlib import Path

from sqsgenerator import write
from sqsgenerator.core import Structure

from src.modeling.base import SimplePoscar, Struct
from src.modeling.supercell import unitcell2file
from src.modeling.model import ModelStruct


def run_modeling(
    name: str,
    poscar: Path,
    outdir: Path,
    supercell_factors: tuple[int, int, int],
    structure_info: dict[str, Any] = {},
    shuffle_seeds: list[int | None] = [None],
    batch_size: int = 1,
    enable_sqs: bool = False,
    iterations: int = 1e7,
) -> list[Path]:
    """
    Run modeling workflow, first generate supercell, then shuffle/sqs structures.

    Args:
        name: Name of the modeling
        poscar: Path to POSCAR file
        outdir: Output directory
        supercell_factors: Factors for supercell (nx, ny, nz)
        structure_info: Structure info
        shuffle_seeds: Seeds just for shuffle
        batch_size: Batch size for modeling
        iterations: Number of iterations for sqsgenerator (default: 1e7)
    Returns:
        List: List of output files
    """
    # Output directory
    outdir = outdir.joinpath(Path(name))
    outdir.mkdir(parents=True, exist_ok=True)  # Ensure the directory is created

    # Get unitcell
    if poscar.is_file():
        logging.info(f"Using poscar {poscar}")
    elif structure_info:
        logging.info(f"Using structure info {structure_info}")
        poscar = unitcell2file(structure_info=structure_info, outdir=outdir)
    else:
        raise ValueError("Either POSCAR file path or structure_info should be provided")

    modeler = ModelStruct(
        name=name,
        poscar=poscar,
        factors=supercell_factors,
        struct_info=structure_info,
        batch_size=batch_size,
    )

    # Choose model engine
    generator = (
        modeler.model_by_sqsgen(iterations=iterations)
        if enable_sqs
        else modeler.model_by_shuffle(seeds=shuffle_seeds)
    )

    # Save structures
    filepaths = []
    for filename, struct in generator:
        poscar = outdir.joinpath(filename)
        SimplePoscar.write_poscar(poscar=poscar, struct=struct, comment=str(filename))
        logging.info(f"Structure saved to {poscar}")
        filepaths.append(poscar)

    logging.info(f"Structures generation completed, files saved to {outdir}")
    return filepaths
