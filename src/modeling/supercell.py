# src/modeling/supercell.py

import logging
from pathlib import Path

import numpy as np
from ase.build import make_supercell as ase_make_supercell

from src.modeling.base import Atom, Struct, SimplePoscar
from src.utils.progress import progress

def _clean_matrix(matrix: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Clean small values in a matrix by setting them to zero.

    Args:
        matrix: Input matrix
        eps: Threshold for small values. Values smaller than `eps` are set to zero

    Returns:
        np.ndarray: Cleaned matrix
    """
    matrix = np.array(matrix)  # Ensure input is a NumPy array
    matrix[np.abs(matrix) < eps] = 0  # Set small values to zero
    return matrix


def make_supercell(struct: Struct, factors: tuple[int, int, int]) -> Struct:
    """Make supercell from struct.

    Args:
        struct: Struct object
        factors: Supercell factors

    Returns:
        struct: Supercell struct
    """
    struct = struct.copy()
    # Original coordinates
    coords = struct.get_coords(direct=True)  # Ensure in direct coordinates

    # Supercell indices
    n, m, p = factors
    i, j, k = np.mgrid[0:n, 0:m, 0:p]  # shape(n, m, p)
    indices = np.stack([i, j, k], axis=-1)  # shape(n, m, p, 3)
    indices = indices.reshape(-1, 3)  # shape(n*m*p, 3)

    # Broadcast
    super_coords = (coords[:, np.newaxis, :] + indices[np.newaxis, :, :]) / [n, m, p]
    super_coords = super_coords.reshape(-1, 3)  # shape(N * n * m * p, 3)
    super_coords = super_coords % 1.0

    # New struct
    matrix = np.array([[n, 0, 0], [0, m, 0], [0, 0, p]])
    new_cell = _clean_matrix(np.dot(struct.cell, matrix))
    new_struct = Struct(cell=new_cell, is_direct=True)
    for idx, coord in progress(enumerate(super_coords), total=len(super_coords)):
        atom = struct[idx // (n * m * p)]
        note = atom.note if atom.note else atom.symbol
        new_atom = Atom(
            index=idx,
            symbol=atom.symbol,
            coord=coord,
            constr=atom.constr,
            note=note,
            meta=atom.meta,
        )
        new_struct.append(new_atom)

    return new_struct


def unitcell2file(structure_info: dict[str, dict], outdir: Path) -> Path:
    """Save structure info to file."""
    # Get cell info
    if "cell" not in structure_info:
        raise ValueError(f"Cell not found in structure_info {structure_info}.")
    cell = np.array(structure_info["cell"])
    if cell.size == 3:
        diagonal_matrix = np.zeros((3, 3))
        np.fill_diagonal(diagonal_matrix, cell)
        cell = diagonal_matrix
    if cell.shape != (3, 3):
        raise ValueError(f"Cell {cell} not a 3x3 matrix.")

    # Get atoms info
    atom_list = []
    for idx, (site, data) in enumerate(structure_info.items()):
        if site == "cell":
            continue
        if "atoms" not in data:
            raise ValueError(f"Atoms_info not found for site {site}.")
        symbol, coords = data["atoms"]
        coords = np.array(coords)
        note = f"{site}-{symbol}"
        substruct = [
            Atom(index=idx, symbol=symbol, coord=coord, note=note) for coord in coords
        ]
        atom_list.extend(substruct)

    struct = Struct(cell=cell, is_direct=True, atom_list=atom_list)
    symbol_count_str = "".join(f"{s}{c}" for s, c in struct.symbol_count)
    output = outdir.joinpath(f"{symbol_count_str}-unitcell.vasp")

    # Ensure the directory exists
    output.parent.mkdir(parents=True, exist_ok=True)

    comment = "Unitcell"
    SimplePoscar.write_poscar(poscar=output, struct=struct, comment=comment)

    logging.info(f"Unitcell saved to {output}")
    return output


def supercell2file(
    poscar: Path, outdir: Path, factors: tuple[int, int, int], by_ase: bool = False
) -> Path:
    """Make supercell and save to file."""
    # Read original POSCAR
    struct = SimplePoscar.read_poscar(poscar)
    logging.debug(f"struct: {struct}")

    if by_ase:
        # Make supercell by ase.make_supercell()
        atoms = SimplePoscar.struct2atoms(struct)
        a, b, c = factors
        matrix = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])
        supercell = ase_make_supercell(atoms, matrix, order="atom-major")
        new_struct = SimplePoscar.atoms2struct(supercell)
    else:
        # Make supercell
        new_struct = make_supercell(struct, factors)

    logging.debug(f"supercell: {new_struct}")

    # Save to file
    factors_str = "x".join(str(f) for f in factors)
    symbol_count_str = "".join(f"{s}{c}" for s, c in struct.symbol_count)
    output = outdir.joinpath(f"{symbol_count_str}-supercell-{factors_str}.vasp")

    # Ensure the directory exists
    output.parent.mkdir(parents=True, exist_ok=True)

    comment = f"Supercell-{factors_str}"
    SimplePoscar.write_poscar(poscar=output, struct=new_struct, comment=comment)

    logging.info(f"Supercell saved to {output}")
    return output
