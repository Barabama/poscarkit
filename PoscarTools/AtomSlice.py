# AtomSlice.py

import logging
import os
from collections import defaultdict
from itertools import groupby

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

from .SimplePoscar import Atoms, SimplePoscar
from .Utils import color_map

basis_map = {
    (0, 0, 1): [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    (1, 1, 0): [(-1, 1, 0), (0, 0, 1), (1, 1, 0)],
    (1, 1, 1): [(-1, 1, 0), (-1, -1, 2), (1, 1, 1)], }


def _normalize(vector: np.ndarray) -> np.ndarray:
    return vector / np.linalg.norm(vector)


def _get_basis(direction: tuple[int, int, int]) -> np.ndarray:
    """Find the base vectors by the direction of the plane.

    Args:
        direction (tuple[int, int, int]): Direction of the plane.
    Returns:
        ndarray: 3 base vectors.
    """
    if direction in basis_map:
        basis = np.array(basis_map[direction])
    else:
        n = np.array(direction)

        # Find two base vectors to the direction
        t0 = np.array([1, 0, 0]) if abs(n[0]) < abs(n[1]) else np.array([0, 1, 0])
        b1 = np.cross(n, t0)
        b2 = np.cross(n, b1)
        basis = np.column_stack([_normalize(v) for v in [b1, b2, n]])

    return basis


def group_by_direction(atoms: Atoms, basis: np.ndarray, precision: int = 6):
    """Group atoms by projection distance along a direction.

    Args:
        atoms (Atoms): Atoms object.
        basis (ndarray): Base vectors.
        precision (int, optional): Number of decimal places to round to. Defaults to 6.
    Yields:
        tuple[float, Atoms]: Projection, layer.
    """
    # Calculate and Round projections
    coords = atoms.cartesian_coords
    projs = np.dot(coords, basis[-1])  # Projections onto direction
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
    # Calculate projections onto direction to get projected coordinates
    b1, b2, n = basis
    coords = layer.cartesian_coords
    n_projs = np.dot(coords, n)  # Projections onto direction
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


def slice2file(filepath: str, direction: tuple[int, int, int]):
    """Slice POSCAR by direction."""
    output = f"{os.path.splitext(os.path.abspath(filepath))[0]}-sliced"
    os.makedirs(output, exist_ok=True)  # Force directory creation
    direct_str = "".join(str(d) for d in direction)

    # Read POSCAR to get atoms
    poscar = SimplePoscar()
    atoms = poscar.read_poscar(filepath)
    symbols_str = "".join(s for s, c in atoms.symbol_count)
    logging.debug(atoms)

    # Get_basis, Regarding direction as z_axis
    basis = _get_basis(direction)  # ndarray([b1, b2, n])
    logging.debug(f"Basis: {basis}")

    # Group atoms by direction
    layers = [ls for ls in group_by_direction(atoms, basis)]
    num_layers = len(layers)
    logging.info(f"Found {num_layers} layers")

    # Save layers as POSCAR and plot layers
    l = len(str(num_layers))
    for i, (proj, layer) in enumerate(tqdm(layers, desc="Processing layers",
                                           total=num_layers, ncols=80), start=1):
        logging.debug(f"Layer {i:0{l}d} proj={proj:.4f}")
        logging.debug(f"layer: {layer}")

        # Save layer to POSCAR file
        filename = os.path.join(output, f"POSCAR-({direct_str})-Layer{i:0{l}d}.vasp")
        comment = f"{symbols_str}-({direct_str})-Layer{i:0{l}d}"
        poscar.write_poscar(filename, layer, comment)

        # Plot layer by base vectors
        imgname = os.path.join(output, f"({direct_str})-Layer{i:0{l}d}.png")
        plot_layer(layer, basis, comment, imgname)
        # break  # for test

    logging.info(f"Results saved in {output}")
