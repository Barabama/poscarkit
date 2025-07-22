# AtomSupercell.py

import logging
import os

import numpy as np
from tqdm import tqdm

from .SimplePoscar import Atom, Atoms, SimplePoscar


def _clean_matrix(matrix: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Clean small values in a matrix by setting them to zero.

    Args:
        matrix (np.ndarray): Input matrix.
        eps (float): Threshold for small values. Values smaller than `eps` are set to zero.

    Returns:
        np.ndarray: Cleaned matrix.
    """
    matrix = np.array(matrix)  # Ensure input is a NumPy array
    matrix[np.abs(matrix) < eps] = 0  # Set small values to zero
    return matrix


def make_supercell(atoms: Atoms, factors: tuple[int, int, int]) -> Atoms:
    """Make supercell from atoms.

    Args:
        atoms (Atoms): Atoms object.
        factors (tuple[int, int, int]): Supercell factors.

    Returns:
        Atoms: Supercell atoms.
    """
    atoms = atoms.copy()
    # Original coordinates
    atoms.switch_coords(direct=True)  # Ensure in direct coordinates
    coords = atoms.direct_coords

    # Supercell indices
    n, m, p = factors
    i, j, k = np.mgrid[0:n, 0:m, 0:p]  # shape(n, m, p)
    indices = np.stack([i, j, k], axis=-1)  # shape(n, m, p, 3)
    indices = indices.reshape(-1, 3)  # shape(n*m*p, 3)

    # Broadcast
    super_coords = (coords[:, np.newaxis, :] + indices[np.newaxis, :, :]) / [n, m, p]
    super_coords = super_coords.reshape(-1, 3)  # shape(N * n * m * p, 3)
    super_coords = super_coords % 1.0

    # New atoms
    matrix = np.array([[n, 0, 0], [0, m, 0], [0, 0, p]])
    new_cell = _clean_matrix(np.dot(atoms.cell, matrix))
    new_atoms = Atoms(cell=new_cell, is_direct=True)
    for idx, coord in enumerate(tqdm(super_coords, ncols=80, desc="Generating Supercell")):
        atom = atoms[idx // (n * m * p)]
        note = atom.note if atom.note else atom.symbol
        new_atoms.append(Atom(index=idx, symbol=atom.symbol, coord=coord,
                              constr=atom.constr, note=note, meta=atom.meta))

    return new_atoms


def unitcell2file(struct_info: dict[str, dict], outdir: str) -> str:
    """Save structure info to file."""
    # Get cell info
    if "cell" not in struct_info:
        raise ValueError("No cell information provided.")
    cell = np.array(struct_info.pop("cell"))
    if cell.size == 3:
        diagonal_matrix = np.zeros((3, 3))
        np.fill_diagonal(diagonal_matrix, cell)
        cell = diagonal_matrix
    if cell.shape != (3, 3):
        raise ValueError("Cell must be a 3-element vector or a 3x3 matrix.")

    # Get atoms info
    atom_list = []
    for idx, (site, data) in enumerate(struct_info.items()):
        if "atoms" not in data:
            raise ValueError(f"No atoms information provided for site({site})")
        symbol, coords = data["atoms"]
        coords = np.array(coords)
        note = f"{site}-{symbol}"
        atom_list.extend([Atom(index=idx, symbol=symbol, coord=coord, note=note)
                          for coord in coords])

    atoms = Atoms(cell=cell, is_direct=True, atom_list=atom_list)
    output = os.path.join(outdir, "POSCAR-unitcell.vasp")
    SimplePoscar.write_poscar(filepath=output, atoms=atoms, comment="Unitcell")

    logging.info(f"Unitcell saved to {output}")
    return output


def supercell2file(filepath: str, outdir: str, factors: tuple[int, int, int]) -> str:
    """Make supercell and save to file."""
    # Read original POSCAR
    atoms = SimplePoscar.read_poscar(filepath)
    logging.debug(f"atoms: {atoms}")

    # Make supercell
    new_atoms = make_supercell(atoms, factors)
    logging.debug(f"supercell: {new_atoms}")

    # ### Make supercell by ASE.make_supercell()
    # from ase.build import make_supercell as ASEmake_supercell
    # ase_atoms = SimplePoscar.to_ase_atoms(atoms)
    # a, b, c = factors
    # matrix = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])
    # supercell = ASEmake_supercell(ase_atoms, matrix, order="atom-major")
    # new_atoms = SimplePoscar.from_ase_atoms(supercell)

    # Save to file
    factors_str = "x".join(str(f) for f in factors)
    output = os.path.join(outdir, f"POSCAR-supercell-{factors_str}.vasp")
    comment = f"Supercell-{factors_str}"
    SimplePoscar.write_poscar(filepath=output, atoms=new_atoms, comment=comment)

    logging.info(f"Supercell saved to {output}")
    return output
