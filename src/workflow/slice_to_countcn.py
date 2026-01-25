# src/workflow/slice_to_countcn.py

import logging
import shutil
from pathlib import Path

from src.modeling.countcn import CNCounter
from src.modeling.slice import Slicer


def slice2files_with_countcn(
    name: str,
    poscar: Path,
    outdir: Path,
    miller_index: tuple[int, int, int],
) -> list[Path]:
    """
    Slice a structure into layers and count the CNs in each layer.

    Args:
        name: Name of the structure.
        poscar: Path to the POSCAR file.
        outdir: Path to the output directory.
        miller_index: Miller index of the slice.

    Returns:
        A list of paths to the output dirs.
    """
    slicer = Slicer(name, poscar, miller_index)
    layer_files = slicer.slice2files(outdir)
    results = []
    for layer_file in layer_files:
        counter = CNCounter(name, layer_file)
        results.append(counter.countCN2files(outdir=layer_file.parent))

    return results
