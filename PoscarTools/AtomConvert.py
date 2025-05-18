# AtomConvert.py

import logging
import os

import numpy as np

from .SimplePoscar import Atoms, SimplePoscar
from .Utils import color_map

basis_map = {
    (0, 0, 1): [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    (1, 1, 0): [(-1, 1, 0), (0, 0, 1), (1, 1, 0)],
    (1, 1, 1): [(-1, 1, 0), (-1, -1, 2), (1, 1, 1)],
}


def _normalize(vector: np.ndarray) -> np.ndarray:
    return vector / np.linalg.norm(vector)


def _get_basis(direction: tuple[int, int, int]) -> np.ndarray:
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


def convert_atoms(atoms: Atoms, basis: np.ndarray):
    # coords = atoms.cartesian_coords
    # basis = np.array([_normalize(v) for v in basis])
    # new_cell = np.dot(atoms.cell, basis.T)
    # new_coords = np.dot(coords, np.linalg.inv(basis))
    # new_atoms = Atoms(cell=new_cell, is_direct=True, atom_list=atoms.atom_list)
    # for atom, coord in zip(new_atoms, new_coords):
    #     atom.coord = coord
    # return new_atoms
    old_cell = atoms.cell
    old_coords = atoms.cartesian_coords
    basis_inv = np.linalg.inv(basis)
    new_cell = (basis_inv @ old_cell.T).T
    new_coords = (basis_inv @ old_coords.T).T
    new_atoms = Atoms(cell=new_cell, is_direct=False, atom_list=atoms.atom_list)
    for atom, coord in zip(new_atoms, new_coords):
        atom.coord = coord
    return new_atoms


def convert2file(filepath: str, direction: tuple[int, int, int]):
    poscar = SimplePoscar()
    atoms = poscar.read_poscar(filepath)
    logging.debug(atoms)

    basis = _get_basis(direction)
    logging.debug(f"Basis: {basis}")

    new_atoms = convert_atoms(atoms, basis)
    logging.debug(new_atoms)

    direct_str = "".join(str(d) for d in direction)
    output = f"{os.path.splitext(filepath)[0]}-{direct_str}.vasp"
    poscar.write_poscar(output, new_atoms)
    logging.info(f"Converted POSCAR saved to {output}")

    return output
