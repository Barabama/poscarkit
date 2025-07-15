# SimplePoscar.py

import logging
import re
from collections import defaultdict, Counter
from collections.abc import Generator, Iterable
from copy import deepcopy
from dataclasses import dataclass

import numpy as np


@dataclass
class Atom:
    index: int
    symbol: str
    coord: np.ndarray
    constr: list[str] | None = None
    comment: str = ""
    # meta: Any

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
            self.atom_list.extend(atoms.atom_list)
        elif isinstance(atoms, list):
            self.atom_list.extend(atoms)
        elif isinstance(atoms, Iterable):
            for atom in atoms:
                self.atom_list.append(atom)

    def insert(self, idx: int, atom: Atom):
        self.atom_list.insert(idx, atom)

    def remove(self, atom: Atom):
        self.atom_list.remove(atom)

    def pop(self, idx: int) -> Atom:
        return self.atom_list.pop(idx)

    def clear(self):
        self.atom_list.clear()

    def copy(self, clean: bool = False) -> 'Atoms':
        atom_list = [] if clean else self.atom_list
        return Atoms(cell=self.cell,
                     is_direct=self.is_direct,
                     atom_list=atom_list)

    def rebuild(self, atom_list: list[Atom] = []) -> 'Atoms':
        return Atoms(cell=self.cell,
                     is_direct=self.is_direct,
                     atom_list=atom_list)

    def sort(self, key="symbol", reverse: bool = False):
        key_funcs = {"symbol": lambda atom: atom.symbol,
                     "coord": lambda atom: tuple(atom.coord),
                     "x": lambda atom: atom.coord[0],
                     "y": lambda atom: atom.coord[1],
                     "z": lambda atom: atom.coord[2],
                     "comment": lambda atom: atom.comment}
        key_func = key_funcs[key]if key in key_funcs else lambda key: key

        self.atom_list.sort(key=key_func, reverse=reverse)

    @property
    def symbol_count(self) -> list[tuple[str, int]]:
        """Returns list: [(symbol, count), ...]."""
        # self.sort(key="symbol")
        symbol_count = list(Counter(atom.symbol for atom in self.atom_list).items())
        return sorted(symbol_count, key=lambda x: x[0])

    @property
    def group_atoms(self) -> Generator[tuple[str, list[Atom]], None, None]:
        """Yields tuple: (symbol, [atoms])."""
        grouped_atoms = defaultdict(list)
        for atom in self.atom_list:
            grouped_atoms[atom.symbol].append(atom)
        yield from ((symbol, atoms) for symbol, atoms in grouped_atoms.items())

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
    def coord_map(self) -> dict[tuple, list[Atom]]:
        """Return a map of coordinates to atoms."""
        coord_map = defaultdict(list)
        for atom in self.atom_list:
            coord_map[tuple(atom.coord.tolist())].append(atom)
        return coord_map

    @property
    def duplicates(self) -> str:
        """Return a list of atoms with duplicate coordinates."""
        results = []
        for coord, atom_list in self.coord_map.items():
            if len(atom_list) <= 1:
                continue
            symbol_count = "".join(
                f"{s}{c}" for s, c in Atoms(self.cell, atom_list=atom_list).symbol_count)
            results.append(f"in {coord} {len(atom_list)} atoms {symbol_count}")
        return "\n".join(results)

    def remove_duplicates(self, reserve_old: bool = False):
        """Remove duplicate atoms."""
        atom_list = []
        for al in self.coord_map.values():
            reserve_idx = 0 if reserve_old else -1
            atom_list.append(al[reserve_idx])
        self.atom_list = atom_list

    def compare(self, atoms2: 'Atoms') -> tuple[bool, str]:
        """Compare two Atoms objects."""
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
        total_atoms = atoms1.copy(clean=False)
        total_atoms.extend(atoms2)
        msg = []
        for coord, atom_list in total_atoms.coord_map.items():
            if len(atom_list) != 2 or atom_list[0].symbol != atom_list[1].symbol:
                msg.append(f"{coord} has atoms {''.join(a.symbol for a in atom_list)}\n")
            else:
                continue

        return (False, "".join(msg)) if msg \
            else (True, f"{atoms1} equals {atoms2}")


class SimplePoscar:
    def __init__(self):
        self.comment_line = ""
        self.selective_dynamics = False
        self.direct_coordinates = True

    def _parse_comment(self, line: str) -> tuple[int, str]:
        result = (-1, "")
        if not line or "#" not in line:
            return result

        comment_part = line.split("#", 1)[1].strip()
        if not comment_part:
            return result

        match = re.search(r"([^#]+)-#(\d+)", comment_part)
        if match:
            try:
                comment = match.group(1) + "-#" + match.group(2)
                return (int(match.group(2)) - 1, comment)
            except (ValueError, AttributeError):
                return result
        return result

    def read_poscar(self, filepath: str) -> Atoms:
        """Read POSCAR file.

        Args:
            filepath (str): Path to POSCAR file.

        Returns:
            Atoms: Atoms from POSCAR.
        """
        with open(filepath, "r") as f:
            lines = f.readlines()

        # Read comment line
        self.comment_line = lines[0].strip()

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
        self.selective_dynamics = "selective" in lines[7].lower()

        # Check coordinate type (Direct or Cartesian)
        coord_type = lines[7 + self.selective_dynamics].strip().lower()[0]
        self.direct_coordinates = coord_type == "d"

        # Read atoms (coordinates, constraints, comment)
        atoms = Atoms(cell=cell, is_direct=self.direct_coordinates)
        start_idx = 8 + self.selective_dynamics
        for symbol, count in zip(symbols, counts):
            for i, line in enumerate(lines[start_idx:start_idx + count]):
                idx, comment = self._parse_comment(line)
                idx = idx if idx != -1 else i
                parts = line.split()
                coord = np.array(list(map(float, parts[:3])))
                constr = parts[3:6] if self.selective_dynamics else []
                comment = comment if comment else \
                    f"{symbol}-#{idx + 1:0{len(str(count))}d}"
                # Apply scale factor to Cartesian coordinates
                if not self.direct_coordinates:
                    coord *= scale
                atoms.append(Atom(index=idx,
                                  symbol=symbol,
                                  coord=coord,
                                  constr=constr,
                                  comment=comment))
            start_idx += count

        # Check for duplicates
        if duplicates := atoms.duplicates:
            logging.warning(f"Duplicate atoms found: {duplicates}")
            atoms.remove_duplicates()

        # Switch to direct coordinates
        self.direct_coordinates = True
        atoms.switch_coords(self.direct_coordinates)

        return atoms

    def write_poscar(self, filepath: str, atoms: Atoms, comment_line: str = ""):
        """Write POSCAR file.

        Args:
            filepath (str): POSCAR file path to write to.
            atoms (Atoms): Atoms to write to POSCAR file.
            comment_line (str, optional): Comment line to write to POSCAR file.
        """
        # Check for duplicates
        if duplicates := atoms.duplicates:
            logging.warning(f"Duplicate atoms found: {duplicates}")
            atoms.remove_duplicates()

        self.comment_line = comment_line if comment_line else self.comment_line
        lines = []

        # Write comment line
        lines.append(self.comment_line)

        # Write scale factor as 1.0
        lines.append(f"{1.0:19.16f}")

        # Write cell vectors
        for vec in atoms.cell:
            lines.append(" " + " ".join(f"{v:21.16f}" for v in vec))

        # Write symbol count
        if len(atoms) == 0:
            symbols = []
            counts = []
        else:
            symbols, counts = zip(*atoms.symbol_count)
        
        lines.append(" " + " ".join(f"{s:>3s}" for s in symbols))
        lines.append(" " + " ".join(f"{c:>3d}" for c in counts))

        # Write if selective dynamics are present
        if self.selective_dynamics:
            lines.append("Selective dynamics")

        # Write direct or cartesian coordinates
        lines.append("Direct" if self.direct_coordinates else "Cartesian")
        atoms.switch_coords(self.direct_coordinates)

        # Write atoms (coordinates, constraint, comment)
        atoms.sort()
        for atom in atoms:
            coord_str = " " + " ".join(f"{c:19.16f}" for c in atom.coord)
            constr_str = " " + " ".join(c for c in atom.constr) \
                if self.selective_dynamics and atom.constr is not None else ""
            comment_str = " # " + atom.comment
            lines.append(f" {coord_str}{constr_str}{comment_str}")

        with open(filepath, "w") as f:
            f.write("\n".join(lines) + "\n")


def compare_poscar(filepath1: str, filepath2: str):
    poscar = SimplePoscar()

    atoms1 = poscar.read_poscar(filepath1)
    atoms2 = poscar.read_poscar(filepath2)

    flag, msg = atoms1.compare(atoms2)
    logging.info(f"{flag}, {msg}")
