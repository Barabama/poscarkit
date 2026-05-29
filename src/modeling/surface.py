# src/modeling/surface.py

import csv
import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase.build.tools import cut

from src.modeling.base import Atom, Struct, SimplePoscar


@dataclass
class Layer:
    """An atomic layer identified in a transformed bulk cell.

    Attributes:
        index: Layer index (0 = lowest z, ascending).
        z_centroid: Mean fractional z coordinate of atoms in this layer.
        atoms: List of Atom objects in this layer.
        composition: Element composition string (e.g. "Co3Ni3V1").
    """
    index: int
    z_centroid: float
    atoms: list
    composition: str


class SurfaceBuilder:
    """Build asymmetric surface slabs from a bulk POSCAR.

    Parameters:
        poscar: Path to bulk POSCAR file.
        miller: Miller indices (h, k, l).
        layers: Number of slab layers (N). N and N+1 slabs are produced.
        vacuum: Total vacuum thickness in Angstrom (default 15.0).
        fix_layers: Number of bottom layers to fix (None = auto).
        fix_z_only: If True, fix only z-direction (TTF) instead of all (FFF).
        outdir: Output directory.
        precision: Decimal precision for z-coordinate grouping.
    """

    def __init__(
        self,
        poscar: Path,
        miller: tuple[int, int, int] = (0, 0, 1),
        layers: int = 3,
        vacuum: float = 15.0,
        fix_layers: int | None = None,
        fix_z_only: bool = False,
        outdir: Path | None = None,
        precision: int = 2,
    ):
        self.poscar = Path(poscar)
        self.miller = miller
        self.n_layers = layers
        self.vacuum = vacuum
        self.fix_layers = fix_layers
        self.fix_z_only = fix_z_only
        self.outdir = Path(outdir) if outdir else Path("output")
        self.precision = precision

        self.vacuum_bottom = 2.0
        self.vacuum_top = vacuum - self.vacuum_bottom

        if self.vacuum < 2.0:
            raise ValueError(
                f"Total vacuum must be >= 2.0 A (minimum 2 A bottom gap). "
                f"Got {vacuum} A."
            )

        self._bulk_struct = SimplePoscar.read_poscar(self.poscar)
        self._transformed: Struct | None = None
        self._layers: list[Layer] = []
        self._gaps: list[float] = []

    def _transform_cell(self) -> Struct:
        """Transform bulk cell so z-axis aligns with surface normal.

        Uses ASE's cut() which preserves per-atom arrays (note, meta, constr_mask).
        Returns the transformed Struct.
        """
        h, k, l = self.miller
        hkl_str = f"{h}{k}{l}"

        atoms = SimplePoscar.struct2atoms(self._bulk_struct)

        if hkl_str == "001":
            a_dir = np.array([1, 0, 0])
            b_dir = np.array([0, 1, 0])
            c_dir = np.array([0, 0, 1])
        elif hkl_str == "110":
            a_dir = np.array([0, 0, -1])
            b_dir = np.array([-1, 1, 0])
            c_dir = np.array([1, 1, 0])
        elif hkl_str == "111":
            a_dir = np.array([1, 1, -2])
            b_dir = np.array([-1, 1, 0])
            c_dir = np.array([1, 1, 1])
        else:
            # General Miller index: build basis vectors from cross products
            n = np.array([h, k, l], dtype=float)
            t0 = np.array([1, 0, 0]) if abs(n[0]) < abs(n[1]) else np.array([0, 1, 0])
            a_dir = np.cross(n, t0)
            b_dir = np.cross(n, a_dir)
            c_dir = n

        transformed_atoms = cut(atoms, a=a_dir, b=b_dir, c=c_dir)
        self._transformed = SimplePoscar.atoms2struct(transformed_atoms)
        return self._transformed

    def _identify_layers(self) -> list[Layer]:
        """Identify atomic layers in the transformed bulk using circular gap detection.

        Returns list of Layer objects ordered by fractional z ascending.
        """
        if self._transformed is None:
            raise RuntimeError("Must call _transform_cell() before _identify_layers()")

        struct = self._transformed
        coords = struct.get_coords(direct=True)
        z_vals = coords[:, 2]

        # Sort atoms by z
        sorted_idx = np.argsort(z_vals)
        sorted_z = z_vals[sorted_idx]

        n_atoms = len(sorted_z)
        if n_atoms == 0:
            return []

        # gaps between consecutive atoms
        gaps = np.diff(sorted_z)
        # wrap-around gap
        wrap_gap = np.array([1.0 - sorted_z[-1] + sorted_z[0]])

        all_gaps = np.concatenate([gaps, wrap_gap])
        threshold = np.median(all_gaps) * 3.0

        # Build layers from large-gap-delimited groups
        groups = []
        current_group = [sorted_idx[0]]

        for i in range(1, n_atoms):
            if gaps[i - 1] > threshold:
                groups.append(current_group)
                current_group = []
            current_group.append(sorted_idx[i])

        # Handle wrap: if wrap_gap is small, merge first and last
        if wrap_gap <= threshold and groups:
            groups[0] = current_group + groups[0]
        else:
            groups.append(current_group)

        # Convert to Layer objects
        result = []
        for li, indices in enumerate(groups):
            atoms_in_layer = [struct[idx] for idx in indices]
            z_centroid = float(np.mean([struct[idx].coord[2] for idx in indices]))
            z_centroid = z_centroid % 1.0

            sym_counts = Counter(a.symbol for a in atoms_in_layer)
            comp = "".join(f"{s}{c}" for s, c in sorted(sym_counts.items()))

            result.append(Layer(
                index=li,
                z_centroid=z_centroid,
                atoms=atoms_in_layer,
                composition=comp,
            ))

        # Sort by z_centroid
        result.sort(key=lambda layer: layer.z_centroid)

        # Re-index after sorting
        for i, layer in enumerate(result):
            layer.index = i

        self._layers = result
        return result
