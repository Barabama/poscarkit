# Supercell.py

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

    l = len(str(len(super_coords)))
    for idx, coord in enumerate(tqdm(super_coords, ncols=80, desc="Generating Supercell")):
        atom = atoms[idx // (n * m * p)]
        new_atom = Atom(index=idx,
                        symbol=atom.symbol,
                        coord=coord,
                        constr=atom.constr,
                        comment=f"{atom.comment}")
                        # comment=f"{atom.symbol}-{idx+1:0{l}d}")
        new_atoms.append(new_atom)

    return new_atoms


def supercell2file(filepath: str, factors: tuple[int, int, int]) -> str:
    """Make supercell and save to file"""
    # Read original POSCAR
    poscar = SimplePoscar()
    atoms = poscar.read_poscar(filepath)
    logging.debug(f"atoms: {atoms}")

    # Make supercell
    logging.info(f"Expanding {filepath}...")
    new_atoms = make_supercell(atoms, factors)
    logging.debug(f"supercell: {new_atoms}")

    # Save to file
    factors_str = "".join(str(f) for f in factors)
    symbols_str = "".join(f"{s}" for s, c in new_atoms.symbol_count)
    output = f"{os.path.splitext(filepath)[0]}-{factors_str}.vasp"
    comment = f"Supercell-{symbols_str}-{factors_str}"
    poscar.write_poscar(output, new_atoms, comment)

    # # ===============================================================================
    # """Make supercell by ASE.make_supercell()"""
    # from ase.atom import Atom as ASEAtom
    # from ase.atoms import Atoms as ASEAtoms
    # from ase.io.vasp import read_vasp, write_vasp
    # from ase.build import make_supercell as ASEmake_supercell
    # atoms: ASEAtoms = read_vasp(filepath)
    # a, b, c = factors
    # matrix = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])
    # supercell: ASEAtoms = ASEmake_supercell(atoms, matrix, order="atom-major")
    # output = f"{os.path.splitext(filepath)[0]}-{''.join(str(f) for f in factors)}.vasp"
    # write_vasp(output, supercell, direct=True, sort=True, vasp5=True)

    logging.info(f"Supercell saved to {output}")
    return output
