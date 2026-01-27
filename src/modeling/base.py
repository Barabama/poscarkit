# src/modeling/base.py

import os
import re
import logging
from copy import deepcopy
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from itertools import groupby
from pathlib import Path
from typing import Any

import numpy as np
from ase import Atoms

from src.utils.progress import progress


@dataclass
class Atom:
    __slots__ = ["index", "symbol", "coord", "constr", "note", "meta"]

    def __init__(
        self,
        index: int,
        symbol: str,
        coord: np.ndarray,
        constr: list[str] = [],
        note: str = "",
        meta: Any = None,
    ):
        self.index = index
        self.symbol = symbol
        self.coord = coord.copy()
        self.constr = constr
        self.note = note
        self.meta = meta

    def __eq__(self, other):
        if isinstance(other, Atom):
            return self.symbol == other.symbol and np.allclose(self.coord, other.coord)
        else:
            return False

    def __hash__(self):
        return hash((self.symbol, tuple(self.coord)))


class Struct:
    """Structure class."""

    # Key functions to sort atoms
    key_funcs = {
        "symbol": lambda atom: atom.symbol,
        "coord": lambda atom: tuple(atom.coord),
        "x": lambda atom: atom.coord[0],
        "y": lambda atom: atom.coord[1],
        "z": lambda atom: atom.coord[2],
        "note": lambda atom: atom.note,
    }

    def __init__(
        self,
        cell: np.ndarray,
        is_direct: bool = True,
        atom_list: list[Atom] = [],
    ):
        self.cell = cell.copy()
        self.is_direct = is_direct
        atom_list = atom_list or []
        self.atom_list = deepcopy(atom_list)

    def __len__(self) -> int:
        return len(self.atom_list)

    def __str__(self) -> str:
        symbol_count = "".join(f"{s}{c}" for s, c in self.symbol_count)
        cell_str = "".join(f"{c:.5f} " for c in self.cell.flatten())
        return f"Struct(cell={cell_str}, is_direct={self.is_direct}, atoms={symbol_count})"

    def __iter__(self):
        yield from self.atom_list

    def __getitem__(self, idx: int) -> Atom:
        return self.atom_list[idx]

    def __setitem__(self, idx: int, atom: Atom):
        self.atom_list[idx] = atom

    def __delitem__(self, idx: int):
        del self.atom_list[idx]

    def append(self, atom: Atom):
        self.atom_list.append(deepcopy(atom))

    def extend(self, atoms: "list[Atom] | Struct"):
        if isinstance(atoms, Struct):
            atom_list = atoms.atom_list
        elif isinstance(atoms, Iterable):
            atom_list = [atom for atom in atoms if isinstance(atom, Atom)]
        self.atom_list.extend(deepcopy(atom_list))

    def insert(self, idx: int, atom: Atom):
        self.atom_list.insert(idx, atom)

    def remove(self, atom: Atom):
        self.atom_list.remove(atom)

    def pop(self, idx: int) -> Atom:
        return self.atom_list.pop(idx)

    def clear(self):
        self.atom_list.clear()

    def classify(self, symbol: str) -> "Struct":
        return self.copy(atom_list=[a for a in self.atom_list if a.symbol == symbol])

    def copy(self, clean: bool = False, atom_list: list[Atom] = []) -> "Struct":
        if clean:
            atom_list = []
        elif not atom_list:
            atom_list = self.atom_list
        return Struct(cell=self.cell, is_direct=self.is_direct, atom_list=atom_list)

    def sort(self, key: str = "symbol", reverse: bool = False) -> "Struct":
        if key not in self.key_funcs:
            raise ValueError(f"Invaild key: {key}. Valids are {list(self.key_funcs.keys())}.")
        # sorted_atoms = sorted(self.atom_list, key=self.key_funcs[key], reverse=reverse)
        # return self.copy(atom_list=sorted_atoms)
        self.atom_list.sort(key=self.key_funcs[key], reverse=reverse)
        return self

    def group_structs(
        self, key: str = "symbol", reverse: bool = False
    ) -> list[tuple[str, "Struct"]]:
        """Returns list[(key, Struct), ...]."""
        sort_struct = self.sort(key=key, reverse=reverse)
        return [
            (str(k), self.copy(atom_list=list(v)))
            for k, v in groupby(sort_struct, key=self.key_funcs[key])
        ]

    @property
    def symbols(self) -> list[str]:
        return [atom.symbol for atom in self.atom_list]

    @property
    def symbol_count(self) -> list[tuple[str, int]]:
        """Returns list[(symbol, count), ...]."""
        symbol_count = list(Counter(self.symbols).items())
        symbol_count.sort(key=lambda x: x[0])
        return symbol_count

    # @property
    # def direct_coords(self) -> np.ndarray:
    #     if not self.is_direct:
    #         self.switch_coords(direct=True)
    #     return np.array([atom.coord for atom in self.atom_list])

    # @property
    # def cartesian_coords(self) -> np.ndarray:
    #     return np.dot(self.direct_coords, self.cell)

    def get_coords(self, direct: bool = True) -> np.ndarray:
        """Switch coordinates to wanted format."""
        coords = np.array([atom.coord for atom in self.atom_list])
        if direct == self.is_direct:
            return coords  # Already in wanted format
        self.is_direct = direct
        inv_cell = np.linalg.inv(self.cell)
        coords = np.dot(coords, inv_cell) if direct else np.dot(coords, self.cell)
        for i, coord in progress(enumerate(coords), len(coords), desc="Switching Coordinates"):
            self.atom_list[i].coord = coord
        return coords

    def remove_duplicates(self, keep_new: bool = False):
        """Remove duplicate atoms."""
        coords = self.get_coords(direct=False)
        # Round coordinates to avoid floating point precision issues
        coords_rounded = np.round(coords, decimals=10)
        # Find unique coordinates (with indices)
        coords_unique, indices = np.unique(coords_rounded, axis=0, return_index=True)
        # Keep last occurrences instead of firse
        if keep_new:
            for i, coord in progress(
                enumerate(coords_unique),
                len(coords_unique),
                desc="Removing Duplicates",
            ):
                matches = np.all(coords_rounded == coord, axis=1)
                indices[i] = np.where(matches)[0][-1]
        self.atom_list = [self.atom_list[i] for i in indices]

    def compare(self, struct2: "Struct") -> tuple[bool, str]:
        """Compare two struct objects. Return (flag, msg)."""
        struct1 = self

        # Compare cell
        cell1, cell2 = struct1.cell, struct2.cell
        if not np.array_equal(cell1, cell2):
            return False, f"cells are different: {cell1} != {cell2}"

        # Compare symbols
        symbols1, symbols2 = struct1.symbols, struct2.symbols
        if symbols1 != symbols2:
            return False, f"symbols are different: {symbols1} != {symbols2}"

        # Compare atoms
        coords1 = struct1.get_coords(direct=True)
        coords2 = struct2.get_coords(direct=True)

        def coord2tuple(coord: np.ndarray) -> tuple[float, float, float]:
            return tuple(round(x, 6) for x in coord)

        atoms1 = set((s, coord2tuple(c)) for s, c in zip(symbols1, coords1))
        atoms2 = set((s, coord2tuple(c)) for s, c in zip(symbols2, coords2))
        missing1 = atoms1 - atoms2
        if missing1:
            return False, f"atoms are different, missing in struct2: {missing1}"
        missing2 = atoms2 - atoms1
        if missing2:
            return False, f"atoms are different, missing in struct1: {missing2}"
        return True, f"{struct1} equals {struct2}"


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
    def read_poscar(poscar: Path | str) -> Struct:
        """Read POSCAR file.

        Args:
            Poscar: Path to POSCAR file or file-like object.

        Returns:
            Struct: Struct from POSCAR
        """
        if isinstance(poscar, Path):
            with open(poscar, "r") as f:
                logging.info(f"Reading POSCAR: {poscar}")
                lines = f.readlines()
        elif os.path.exists(poscar):
            with open(poscar, "r") as f:
                logging.info(f"Reading POSCAR: {poscar}")
                lines = f.readlines()
        else:
            lines = poscar.splitlines()

        # Read comment line
        comment = lines[0].strip()
        logging.info(f"Comment: {comment}")

        # Read scale factor
        scale = np.array(list(map(float, lines[1].strip().split())))

        # Read cell vectors
        cell = np.array([list(map(float, line.split())) for line in lines[2:5]])

        # Apply scale factor to cell vectors
        scale = scale if scale[0] >= 0.0 else np.cbrt(-1.0 * scale / np.linalg.det(cell))

        cell *= scale  # Working for both one and three scale factors

        # Read symbols and counts
        symbols = lines[5].split()
        counts = list(map(int, lines[6].split()))
        symbol_counts = []
        for symbol, count in zip(symbols, counts):
            symbol_counts.extend([symbol] * count)

        counts_total = sum(counts)

        # Check if selective dynamics is present
        constrainted = "selective" in lines[7].lower()

        # Check coordinate type (Direct or Cartesian)
        coord_type = lines[7 + constrainted].strip().lower()[0]
        is_direct = coord_type == "d"

        # Read atoms (coordinates, constraints, note)
        atom_list = []
        start_idx = 8 + constrainted
        for i, line in progress(
            enumerate(lines[start_idx : start_idx + counts_total]),
            counts_total,
            desc="Reading atoms",
        ):
            symbol = symbol_counts[i]
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
            atom_list.append(
                Atom(
                    index=idx,
                    symbol=symbol,
                    coord=coord,
                    constr=constr,
                    note=note,
                    meta=meta,
                )
            )
        struct = Struct(cell=cell, is_direct=is_direct, atom_list=atom_list)

        # Check for duplicates
        logging.info(f"Formatting struct: {struct}")
        struct.remove_duplicates()

        # Switch to direct coordinates
        logging.info(f"Switching to direct coordinates: {struct}")
        is_direct = True
        struct.get_coords(is_direct)

        return struct

    @staticmethod
    def write_poscar(
        poscar: Path,
        struct: Struct,
        comment: str = "",
        is_direct: bool = True,
        constrainted: bool = True,
    ):
        """Write POSCAR file.

        Args:
            Poscar: POSCAR file path to write to
            struct: Struct to write to POSCAR file
            comment: Comment line to write to POSCAR file
            is_direct: Whether to write in direct coordinates. Defaults to True
            constrainted: Whether to write constrainted POSCAR file. Defaults to True
        """
        # Check for duplicates
        logging.info(f"Formatting struct: {struct}")
        struct.remove_duplicates()

        lines = []

        # Write comment line
        lines.append(comment.split("\n")[0])

        # Write scale factor as 1.0
        lines.append(f"{1.0:19.16f}")

        # Write cell vectors
        for vec in struct.cell:
            lines.append(" " + " ".join(f"{v:21.16f}" for v in vec))

        # Write symbol count
        symbols, counts = zip(*struct.symbol_count) if len(struct) > 0 else ([], [])
        lines.append(" " + " ".join(f"{s:>3s}" for s in symbols))
        lines.append(" " + " ".join(f"{c:>3d}" for c in counts))

        # Write if selective dynamics are present
        if constrainted and any(a.constr for a in struct):
            lines.append("Selective dynamics")

        # Write direct or cartesian coordinates
        s = "Direct" if is_direct else "Cartesian"
        logging.info(f"Switching to {s.lower()} coordinates")
        lines.append(s)
        struct.get_coords(is_direct)

        # Write atoms (coordinates, constraint, note)
        for atom in progress(struct.sort(), len(struct), desc="Writing atoms"):
            coord_str = " " + " ".join(f"{c:19.16f}" for c in atom.coord)
            constr_str = (
                " " + " ".join(c for c in atom.constr)
                if constrainted and atom.constr is not None
                else ""
            )
            meta = atom.meta if atom.meta is not None else ""
            note_str = f"# {atom.note}-" if atom.note is not None else ""
            comment_str = f" {note_str}#{atom.index + 1:0{len(str(len(struct)))}d} {meta}"
            lines.append(f" {coord_str}{constr_str}{comment_str}")

        # Write POSCAR file
        logging.info(f"Writing POSCAR: {poscar}")
        with open(poscar, "w") as f:
            f.write("\n".join(lines) + "\n")

    @staticmethod
    def struct2atoms(struct: Struct) -> Atoms:
        """Convert Strcut to Atoms."""
        symbols = struct.symbols
        cell = struct.cell
        positions = struct.get_coords(False)
        return Atoms(symbols=symbols, cell=cell, positions=positions, pbc=True)

    @staticmethod
    def atoms2struct(atoms: Atoms) -> Struct:
        """Convert Atoms to Struct."""
        cell = np.array(atoms.get_cell().copy())
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        atom_list = []
        for idx, (symbol, position) in enumerate(zip(symbols, positions)):
            atom_list.append(Atom(index=idx, symbol=symbol, coord=position))
        struct = Struct(cell=cell, is_direct=False, atom_list=atom_list)
        return struct

    @staticmethod
    def compare_poscar(poscar1: Path | str, poscar2: Path | str):
        """Compare two POSCAR files. Return (flag, msg).

        Args:
            poscar1: First POSCAR file path
            poscar2: Second POSCAR file path
        """
        logging.info(f"Comparing {poscar1} and {poscar2}")
        struct1 = SimplePoscar.read_poscar(poscar1)
        struct2 = SimplePoscar.read_poscar(poscar2)
        flag, msg = struct1.compare(struct2)
        if flag:
            logging.info(f"The two POSCAR files are the same.")
        else:
            logging.info(f"The two POSCAR files are different.")
            logging.info(msg)

    @staticmethod
    def merge_poscar(poscar1: Path | str, poscar2: Path | str, outdir: Path | str) -> Path:
        """Merge two POSCAR files.

        Args:
            poscar1: First POSCAR file path
            poscar2: Second POSCAR file path
            outdir: Output directory path
        Returns:
            Path: Merged POSCAR file path
        """
        poscar1 = Path(poscar1) if isinstance(poscar1, str) else poscar1
        poscar2 = Path(poscar2) if isinstance(poscar2, str) else poscar2
        outdir = Path(outdir) if isinstance(outdir, str) else outdir
        name1 = poscar1.stem
        name2 = poscar2.stem
        logging.info(f"Merging {poscar1} and {poscar2}")
        struct1 = SimplePoscar.read_poscar(poscar1)
        struct2 = SimplePoscar.read_poscar(poscar2)
        struct = struct1.copy(atom_list=struct1.atom_list + struct2.atom_list)
        output = outdir.joinpath(f"POSCAR-merged-{name1}-{name2}.vasp")
        comment = f"Merged {name1} and {name2}"
        SimplePoscar.write_poscar(poscar=output, struct=struct, comment=comment)
        return output

    @staticmethod
    def separate_poscar(
        poscar: Path | str, outdir: Path | str, key: str = "note"
    ) -> list[Path]:
        """Separate a POSCAR file by a key.

        Args:
            poscar: POSCAR file path
            outdir: Output directory path
            key: Key to separate by (default: "note")
        Returns:
            List: List of output files
        """
        poscar = Path(poscar) if isinstance(poscar, str) else poscar
        outdir = Path(outdir) if isinstance(outdir, str) else outdir
        logging.info(f"Separating {poscar} by {key}")
        struct = SimplePoscar.read_poscar(poscar)
        outputs = []
        for key, substruct in struct.group_structs(key=key):
            output = outdir.joinpath(f"POSCAR-group-{key}.vasp")
            comment = f"Group-{key}"
            SimplePoscar.write_poscar(poscar=output, struct=substruct, comment=comment)
            outputs.append(output)
        return outputs
