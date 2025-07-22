# SimplePoscar.py

from itertools import groupby
import logging
import os
import re
from collections import Counter
from collections.abc import Iterable
from copy import copy, deepcopy
from dataclasses import dataclass
from typing import Any

import numpy as np
from ase import Atoms as ASEAtoms


key_funcs = {"symbol": lambda atom: atom.symbol,
             "coord": lambda atom: tuple(atom.coord),
             "x": lambda atom: atom.coord[0],
             "y": lambda atom: atom.coord[1],
             "z": lambda atom: atom.coord[2],
             "note": lambda atom: atom.note}


@dataclass
class Atom:
    index: int
    symbol: str
    coord: np.ndarray
    constr: list[str] | None = None
    note: str = ""
    meta: Any = None

    def __eq__(self, other):
        if isinstance(other, Atom):
            return (self.symbol == other.symbol and np.allclose(self.coord, other.coord))
        return False

    def __hash__(self):
        return hash((self.symbol, tuple(self.coord)))


class Atoms:
    def __init__(self, cell: np.ndarray, is_direct: bool = True, atom_list: list[Atom] | None = None):
        self.cell = cell.copy()
        self.is_direct = is_direct
        self.atom_list = deepcopy(atom_list) if atom_list is not None else []

    def __len__(self) -> int:
        return len(self.atom_list)

    def __str__(self) -> str:
        symbol_count = "".join(f"{s}{c}" for s, c in self.symbol_count)
        return f"Atoms(cell={self.cell}, is_direct={self.is_direct}, atoms={symbol_count})"

    def __iter__(self):
        yield from self.atom_list

    def __getitem__(self, idx: int) -> Atom:
        return self.atom_list[idx]

    def __setitem__(self, idx: int, atom: Atom):
        self.atom_list[idx] = atom

    def __delitem__(self, idx: int):
        del self.atom_list[idx]

    def append(self, atom: Atom):
        self.atom_list.append(atom)

    def extend(self, atoms: "list[Atom] | Atoms"):
        if isinstance(atoms, Atoms):
            atom_list = atoms.atom_list
        elif isinstance(atoms, Iterable):
            atom_list = [atom for atom in atoms if isinstance(atom, Atom)]
        self.atom_list.extend(atom_list)

    def insert(self, idx: int, atom: Atom):
        self.atom_list.insert(idx, atom)

    def remove(self, atom: Atom):
        self.atom_list.remove(atom)

    def pop(self, idx: int) -> Atom:
        return self.atom_list.pop(idx)

    def clear(self):
        self.atom_list.clear()

    def copy(self, clean: bool = False, atom_list: list[Atom] | None = None) -> "Atoms":
        if clean:
            atom_list = []
        elif atom_list is None:
            atom_list = self.atom_list
        return Atoms(cell=self.cell, is_direct=self.is_direct, atom_list=deepcopy(atom_list))

    def sort(self, key: str = "symbol", reverse: bool = False) -> "Atoms":
        if key not in key_funcs:
            raise ValueError(f"Invalid key: {key}. Valid keys are {list(key_funcs.keys())}")
        return self.copy(atom_list=sorted(self.atom_list, key=key_funcs[key], reverse=reverse))

    def group_atoms(self, key: str = "symbol", reverse: bool = False) -> list[tuple[str, "Atoms"]]:
        """Returns: list: [(key, Atoms), ...]."""
        sort_atoms = self.sort(key=key, reverse=reverse)
        return [(str(k), self.copy(atom_list=list(v)))
                for k, v in groupby(sort_atoms, key=key_funcs[key])]

    @property
    def symbol_count(self) -> list[tuple[str, int]]:
        """Returns list: [(symbol, count), ...]."""
        symbol_count = list(Counter(atom.symbol for atom in self.atom_list).items())
        return sorted(symbol_count, key=lambda x: x[0])

    @property
    def direct_coords(self) -> np.ndarray:
        if not self.is_direct:
            self.switch_coords(direct=True)
        return np.array([atom.coord for atom in self.atom_list])

    @property
    def cartesian_coords(self) -> np.ndarray:
        return np.dot(self.direct_coords, self.cell)

    def switch_coords(self, direct: bool = True):
        """Switch coordinates to wanted format."""
        if direct == self.is_direct:
            return  # Already in correct format
        inv_cell = np.linalg.inv(self.cell)
        for atom in self.atom_list:
            atom.coord = np.dot(atom.coord, inv_cell) % 1.0 if direct \
                else np.dot(atom.coord, self.cell)
        self.is_direct = direct

    @property
    def duplicates(self) -> str:
        """Return a list of atoms with duplicate coordinates."""
        results = []
        for coord, subatoms in self.group_atoms(key="coord"):
            if len(subatoms) <= 1:
                continue
            symbol_count = "".join(f"{s}{c}" for s, c in subatoms.symbol_count)
            results.append(f"in {coord} {len(subatoms)} atoms {symbol_count}")
        return "\n".join(results)

    def remove_duplicates(self, keep_old: bool = False):
        """Remove duplicate atoms."""
        atom_list = []
        for coord, subatoms in self.group_atoms(key="coord"):
            keep_idx = 0 if keep_old else -1
            atom_list.append(subatoms[keep_idx])
        self.atom_list = atom_list

    def compare(self, atoms2: "Atoms") -> tuple[bool, str]:
        """Compare two Atoms objects. Return tuple(flag, msg)."""
        atoms1 = self

        # Check cell
        cell1 = atoms1.cell
        cell2 = atoms2.cell
        if not np.array_equal(cell1, cell2):
            return False, f"{cell1} != {cell2}"

        # Check symbols counts
        sc1 = atoms1.symbol_count
        sc2 = atoms2.symbol_count
        if sc1 != sc2:
            return False, f"{sc1} != {sc2}"

        # Check atoms
        atoms1.switch_coords(direct=True)
        atoms2.switch_coords(direct=True)
        total_atoms = atoms1.copy()
        total_atoms.extend(atoms2)
        msg = []
        for coord, subatoms in total_atoms.group_atoms(key="coord"):
            if len(subatoms) <= 1 and subatoms[0].symbol != subatoms[1].symbol:
                continue
            msg.append(f"{coord} has atoms {''.join(a.symbol for a in subatoms)}\n")

        return (False, "".join(msg)) if msg \
            else (True, f"{atoms1} equals {atoms2}")


class SimplePoscar:

    @staticmethod
    def _parse_comment(line: str) -> tuple[str, int, Any]:
        """Try to return note, index, meta of a line."""
        result = ("", -1, None)
        if not line or "#" not in line:
            return result

        comment_str = line.split("#", 1)[1].strip()
        if not comment_str:
            return result

        match = re.search(r"(\d+[a-z]-[A-Za-z]+)-#(\d+)(.)?", comment_str)
        if match:
            try:
                return (match.group(1), int(match.group(2)) - 1, match.group(3))
            except (ValueError, AttributeError):
                return result
        return result

    @staticmethod
    def read_poscar(filepath: str) -> Atoms:
        """Read POSCAR file.

        Args:
            filepath (str): Path to POSCAR file.

        Returns:
            Atoms: Atoms from POSCAR.
        """
        with open(filepath, "r") as f:
            lines = f.readlines()

        # Read comment line
        comment = lines[0].strip()

        # Read scale factor
        scale = np.array(list(map(float, lines[1].strip().split())))

        # Read cell vectors
        cell = np.array([list(map(float, line.split())) for line in lines[2:5]])

        # Apply scale factor to cell vectors
        scale = scale if scale[0] >= 0.0 else \
            np.cbrt(-1.0 * scale / np.linalg.det(cell))
        cell *= scale  # Working for both one and three scale factors

        # Read symbols and counts
        symbols = lines[5].split()
        counts = list(map(int, lines[6].split()))

        # Check if selective dynamics is present
        constrainted = "selective" in lines[7].lower()

        # Check coordinate type (Direct or Cartesian)
        coord_type = lines[7 + constrainted].strip().lower()[0]
        is_direct = coord_type == "d"

        # Read atoms (coordinates, constraints, note)
        atoms = Atoms(cell=cell, is_direct=is_direct)
        start_idx = 8 + constrainted
        for symbol, count in zip(symbols, counts):
            for i, line in enumerate(lines[start_idx:start_idx + count]):
                parts = line.split()
                note, idx, meta = SimplePoscar._parse_comment(line)
                idx = idx if idx != -1 else i
                coord = np.array(list(map(float, parts[:3])))
                constr = parts[3:6] if constrainted else []
                note = note if note else symbol
                # f"{symbol}-#{idx + 1:0{len(str(count))}d}"
                # Apply scale factor to Cartesian coordinates
                if not is_direct:
                    coord *= scale
                atoms.append(Atom(index=idx, symbol=symbol, coord=coord,
                                  constr=constr, note=note, meta=meta))
            start_idx += count

        # Check for duplicates
        if duplicates := atoms.duplicates:
            logging.warning(f"Duplicate atoms found: {duplicates}")
            atoms.remove_duplicates()

        # Switch to direct coordinates
        is_direct = True
        atoms.switch_coords(is_direct)

        return atoms

    @staticmethod
    def write_poscar(filepath: str, atoms: Atoms, comment: str = "",
                     is_direct: bool = True, constrainted: bool = True):
        """Write POSCAR file.

        Args:
            filepath (str): POSCAR file path to write to.
            atoms (Atoms): Atoms to write to POSCAR file.
            comment (str, optional): Comment line to write to POSCAR file.
            is_direct (bool, optional): Whether to write in direct coordinates. Defaults to True.
            constrainted (bool, optional): Whether to write constrainted POSCAR file. Defaults to True.
        """
        # Check for duplicates
        if duplicates := atoms.duplicates:
            logging.warning(f"Duplicate atoms found: {duplicates}")
            atoms.remove_duplicates()

        lines = []

        # Write comment line
        lines.append(comment.split("\n")[0])

        # Write scale factor as 1.0
        lines.append(f"{1.0:19.16f}")

        # Write cell vectors
        for vec in atoms.cell:
            lines.append(" " + " ".join(f"{v:21.16f}" for v in vec))

        # Write symbol count
        symbols, counts = zip(*atoms.symbol_count) if len(atoms) > 0 else ([], [])
        lines.append(" " + " ".join(f"{s:>3s}" for s in symbols))
        lines.append(" " + " ".join(f"{c:>3d}" for c in counts))

        # Write if selective dynamics are present
        if constrainted and any(a.constr is not None for a in atoms):
            lines.append("Selective dynamics")

        # Write direct or cartesian coordinates
        lines.append("Direct" if is_direct else "Cartesian")
        atoms.switch_coords(is_direct)

        # Write atoms (coordinates, constraint, note)
        for atom in atoms.sort():
            coord_str = " " + " ".join(f"{c:19.16f}" for c in atom.coord)
            constr_str = " " + " ".join(c for c in atom.constr) \
                if constrainted and atom.constr is not None else ""
            meta = atom.meta if atom.meta is not None else ""
            comment_str = f" # {atom.note}-#{atom.index + 1:0{len(str(len(atoms)))}d} {meta}"
            lines.append(f" {coord_str}{constr_str}{comment_str}")

        with open(filepath, "w") as f:
            f.write("\n".join(lines) + "\n")

    @staticmethod
    def to_ase_atoms(atoms: Atoms) -> ASEAtoms:
        """Convert Atoms to ASEAtoms."""
        symbols = [atom.symbol for atom in atoms]
        cell = atoms.cell
        positions = atoms.cartesian_coords
        return ASEAtoms(symbols=symbols, cell=cell, positions=positions, pbc=True)

    @staticmethod
    def from_ase_atoms(ase_atoms: ASEAtoms) -> Atoms:
        """Convert ASEAtoms to Atoms."""
        cell = np.array(ase_atoms.get_cell().copy())
        atom_list = []
        for idx, (symbol, position) in enumerate(
                zip(ase_atoms.get_chemical_symbols(), ase_atoms.get_positions())):
            atom_list.append(Atom(index=idx, symbol=symbol, coord=copy(position)))
        atoms = Atoms(cell=cell, is_direct=False, atom_list=atom_list)
        return atoms

    @staticmethod
    def compare_poscar(filepath1: str, filepath2: str):
        atoms1 = SimplePoscar.read_poscar(filepath1)
        atoms2 = SimplePoscar.read_poscar(filepath2)

        flag, msg = atoms1.compare(atoms2)
        logging.info(f"{flag}, {msg}")

    @staticmethod
    def merge_poscar(filepath1: str, filepath2: str, outdir: str):
        atoms1 = SimplePoscar.read_poscar(filepath1)
        atoms2 = SimplePoscar.read_poscar(filepath2)
        atoms = atoms1.copy(atom_list=atoms1.atom_list.extend(atoms2.atom_list))
        output = os.path.join(outdir, f"POSCAR-merged.vasp")
        SimplePoscar.write_poscar(filepath=output, atoms=atoms, comment="Merged")

    @staticmethod
    def separate_poscar(filepath: str, outdir: str, key: str = "note"):
        atoms = SimplePoscar.read_poscar(filepath)
        for key, subatoms in atoms.group_atoms(key=key):
            output = os.path.join(outdir, f"POSCAR-group-{key}.vasp")
            SimplePoscar.write_poscar(filepath=output, atoms=subatoms, comment=f"Group-{key}")
