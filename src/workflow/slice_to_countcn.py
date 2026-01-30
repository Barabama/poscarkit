# src/workflow/slice_to_countcn.py

import logging
import shutil
from pathlib import Path

from src.modeling.base import SimplePoscar
from src.modeling.countcn import CNCounter
from src.modeling.slice import Slicer
from src.utils.progress import progress


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
    basis = slicer.basis
    transfd = slicer.transformed

    # Get basis
    logging.info(f"Basis: {basis}")

    # Output directory
    miller_index_str = "".join(str(d) for d in miller_index)
    outdir = Path(outdir) if isinstance(outdir, str) else outdir
    outdir = outdir.joinpath(f"{name}-sliced({miller_index_str})")
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Transform
    output = outdir.joinpath(f"Transformed({miller_index_str}).vasp")
    SimplePoscar.write_poscar(poscar=output, struct=transfd, comment=str(output.stem))

    # Group by normal
    layers = [ls for ls in slicer.group_by_normal()]
    num_layers = len(layers)
    logging.info(f"Number of layers: {num_layers}")
    ll = len(str(num_layers))
    results = []
    for i, (proj, layer) in progress(
        enumerate(layers, start=1),
        total=num_layers,
        desc="Slicing layers",
    ):
        logging.info(f"Layer {i}: {proj:.2f}")
        logging.info(f"Layer {layer}")
        # Save layer
        layer_name = f"Transformed({miller_index_str})-layer{i:0{ll}d}"
        output = outdir.joinpath(f"{layer_name}.vasp")
        SimplePoscar.write_poscar(poscar=output, struct=layer, comment=str(output.stem))

        # Count CN
        counter = CNCounter(layer_name, poscar=output)
        cn_count_dir = counter.countCN2files(outdir=output.parent)
        results.append(cn_count_dir)

        # Export to Excel
        xlspath = output.with_suffix(".xlsx")
        slicer.export_layer_xls(layer, xlspath)

        # Plot layer
        imgpath = output.with_suffix(".png")
        title = f"Transformed({' '.join(str(d) for d in miller_index)})-layer{i:0{ll}d}"
        slicer.plot_layer(
            imgpath=imgpath,
            title=title,
            layer=layer,
            pair_counts=counter.pair_counts,
        )

    logging.info(f"Sliced Saved to {outdir}")
    return results
