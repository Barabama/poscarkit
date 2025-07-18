# AtomSlice.py

import logging
import os
from collections import defaultdict
from itertools import groupby

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

from .SimplePoscar import Atoms, read_poscar, write_poscar, to_ase_atoms, from_ase_atoms
from .Utils import color_map

basis_map = {
    (0, 0, 1): [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    (1, 1, 0): [(1, -1, 0), (0, 0, -1), (1, 1, 0)],
    (1, 1, 1): [(1, -1, 0), (1, 1, -2), (1, 1, 1)], }


def _normalize(vector: np.ndarray) -> np.ndarray:
    return vector / np.linalg.norm(vector)


def _get_basis(miller_index: tuple[int, int, int]) -> np.ndarray:
    """Find the base vectors by the miller_index of the plane.

    Args:
        miller_index (tuple[int, int, int]): The miller_index of the plane.
    Returns:
        ndarray: 3 base vectors.
    """
    if miller_index in basis_map:
        basis = np.array(basis_map[miller_index])
    else:
        n = np.array(miller_index)
        # Find two base vectors to the miller index
        t0 = np.array([1, 0, 0]) if abs(n[0]) < abs(n[1]) else np.array([0, 1, 0])
        b1 = np.cross(n, t0)
        b2 = np.cross(n, b1)
        basis = np.column_stack([_normalize(v) for v in [b1, b2, n]])
    return basis


def _convert(atoms: Atoms, basis: np.ndarray) -> Atoms:
    from ase.atoms import Atoms as ASEAtoms
    from ase.build.tools import cut
    ase_atoms: ASEAtoms = to_ase_atoms(atoms)
    a, b, c = basis
    converted = cut(ase_atoms, a, b, c, maxatoms=len(atoms))
    # from ase.io.vasp import write_vasp
    # from ase.visualize import view
    # view(ase_atoms)
    # view(converted)
    # write_vasp("original.vasp", ase_atoms, direct=True, sort=True, vasp5=True)
    # write_vasp("converted.vasp", converted, direct=True, sort=True, vasp5=True)
    return from_ase_atoms(converted)


def group_by_normal(atoms: Atoms, precision: int = 2):
    """Group atoms by projection distance along the normal of base vectors.

    Args:
        atoms (Atoms): Atoms object.
        precision (int, optional): Number of decimal places to round to. Defaults to 6.
    Yields:
        tuple[float, Atoms]: Projection, layer.
    """
    basis = np.array([(1, 0, 0), (0, 1, 0), (0, 0, 1)])
    # Calculate and Round projections
    # coords = atoms.cartesian_coords
    coords = atoms.direct_coords
    projs = np.dot(coords, basis[-1])  # Projections onto the normal
    projs = np.round(projs, precision)

    # Sort atoms based on rounded projections
    sorted_indices = np.argsort(projs)
    for proj, group in groupby(sorted_indices, key=lambda x: projs[x]):
        layer = atoms.rebuild([atoms[i] for i in group])
        yield proj, layer


def plot_layer(layer: Atoms, basis: np.ndarray, title: str, filepath: str):
    """Plot layer by base vectors.

    Args:
        layer (list[Atom]): list of atoms in layer.
        basis (ndarray): Base vectors.
        title (str): Title of plot.
        filepath (str): File path to save plot.
    """
    # Calculate projections onto the normal to get projected coordinates
    b1, b2, n = basis
    coords = layer.cartesian_coords
    n_projs = np.dot(coords, n)  # Projections onto the normal
    p_projs = coords - np.outer(n_projs, n)  # Projections onto plane
    # xs = np.dot(p_projs, b1)  # Components of on b1
    # ys = np.dot(p_projs, b2)  # Components of on b2
    proj_coords = np.column_stack((np.dot(p_projs, b1), np.dot(p_projs, b2)))

    # Group projected coordinates by symbol
    symbol_coords = defaultdict(list)
    for atom, coord in zip(layer, proj_coords):
        symbol_coords[atom.symbol].append(coord)

    # Plot layer with projected coordinates
    plt.figure(figsize=(6, 6))
    for symbol, coords in symbol_coords.items():
        color = color_map.get(symbol, "magenta")
        x, y = zip(*coords)
        plt.scatter(x, y, marker="o", s=10, color=color, alpha=1.0, label=symbol)

    plt.title(title)
    plt.xlabel(f"[{' '.join(str(v) for v in basis[0])}] Coordinate (Å)")
    plt.ylabel(f"[{' '.join(str(v) for v in basis[1])}] Coordinate (Å)")
    plt.axis("equal")
    plt.grid()
    plt.legend(title="Symbols", bbox_to_anchor=(1, 1), loc="upper left")
    # plt.tight_layout(rect=[0, 0, 1, 0])
    plt.savefig(filepath, bbox_inches="tight")
    plt.close()


def slice2file(filepath: str, miller_index: tuple[int, int, int]) -> str:
    """Slice POSCAR by the miller index."""
    miller_index_str = "".join(str(d) for d in miller_index)
    output = f"{os.path.splitext(os.path.abspath(filepath))[0]}-({miller_index_str})-sliced"
    os.makedirs(output, exist_ok=True)  # Force directory creation

    # Read POSCAR to get atoms
    atoms = read_poscar(filepath)
    symbols_str = "".join(s for s, c in atoms.symbol_count)
    logging.debug(atoms)

    # Get_basis, Regarding miller index as the normal
    basis = _get_basis(miller_index)  # ndarray([b1, b2, n])
    logging.debug(f"Basis: {basis}")

    # Convert atoms alone with basis
    converted = _convert(atoms, basis)
    comment = f"{symbols_str}-({miller_index_str})-converted"
    filename = os.path.join(output, f"POSCAR-{comment}.vasp")
    write_poscar(filename, converted, comment)

    # Group atoms by the normal
    layers = [ls for ls in group_by_normal(converted)]
    num_layers = len(layers)
    logging.info(f"Found {num_layers} layers")

    # Save layers as POSCAR and plot layers
    l = len(str(num_layers))
    for i, (proj, layer) in enumerate(tqdm(layers, desc="Processing layers",
                                           total=num_layers, ncols=80), start=1):
        logging.debug(f"Layer {i:0{l}d} proj={proj:.4f}")
        logging.debug(f"layer: {layer}")

        # Save layer to POSCAR file
        comment = f"{symbols_str}-({miller_index_str})-Layer{i:0{l}d}"
        filename = os.path.join(output, f"POSCAR-{comment}.vasp")
        write_poscar(filename, layer, comment)

        # Plot layer by base vectors
        imgname = os.path.join(output, f"{comment}.png")
        plot_layer(layer, basis, comment, imgname)
        # break  # for test

    logging.info(f"Results saved in {output}")
    return output
