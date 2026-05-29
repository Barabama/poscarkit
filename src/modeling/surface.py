# src/modeling/surface.py

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
        composition: Element composition string (e.g. "Co3Ni3V").
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
        """Identify atomic layers in the transformed bulk.

        Groups atoms by their rounded fractional z-coordinate, which
        reliably separates crystallographic layers regardless of cell size.
        Returns list of Layer objects ordered by fractional z ascending.
        """
        if self._transformed is None:
            raise RuntimeError("Must call _transform_cell() before _identify_layers()")

        struct = self._transformed
        coords = struct.get_coords(direct=True)
        z_vals = coords[:, 2]

        n_atoms = len(z_vals)
        if n_atoms == 0:
            return []

        # Round z to configured precision
        z_rounded = np.round(z_vals, decimals=self.precision)

        # Group atoms by rounded z
        unique_z = np.sort(np.unique(z_rounded))

        # Merge top and bottom if they represent the same layer
        # (periodic boundary: z ~ 0.0 and z ~ 1.0 round to 0.00 and 1.00)
        if len(unique_z) >= 2:
            gap = 1.0 - unique_z[-1] + unique_z[0]
            # If the wrap gap is <= the precision unit, merge
            wrap_threshold = 1.5 * 10 ** (-self.precision)
            if gap <= wrap_threshold:
                # Merge: atoms at z ~ 1.0 become part of z ~ 0.0 layer
                z_rounded[z_rounded >= unique_z[-1] - 1e-10] = unique_z[0]
                unique_z = unique_z[:-1]  # Remove the top group

        # Build Layer objects from groups
        result = []
        for li, zu in enumerate(unique_z):
            mask = np.isclose(z_rounded, zu, atol=0.5 * 10 ** (-self.precision))
            indices = np.where(mask)[0]
            atoms_in_layer = [struct[i] for i in indices]
            # Handle PBC merge: atoms near z=1.0 wrap to near z=0.0 for correct centroid
            z_layer = z_vals[indices].copy()
            if np.max(z_layer) - np.min(z_layer) > 0.5:
                z_layer[z_layer > 0.5] -= 1.0
            z_centroid = float(np.mean(z_layer)) % 1.0

            sym_counts = Counter(a.symbol for a in atoms_in_layer)
            comp = "".join(
                f"{s}{c}" if c > 1 else s
                for s, c in sorted(sym_counts.items())
            )

            result.append(Layer(
                index=li,
                z_centroid=z_centroid,
                atoms=atoms_in_layer,
                composition=comp,
            ))

        # Sort by z_centroid ascending
        result.sort(key=lambda layer: layer.z_centroid)
        for i, layer in enumerate(result):
            layer.index = i

        self._layers = result
        return result

    def _find_gaps(self) -> list[float]:
        """Find midpoints (cutting planes) between consecutive layers.

        Returns list of midpoint z-values in fractional coordinates.
        gap[i] is the midpoint between layer[i] and layer[(i+1) % n].
        Called "gaps" because each midpoint defines where to cut between layers.
        """
        if not self._layers:
            raise RuntimeError("Must call _identify_layers() before _find_gaps()")

        n = len(self._layers)
        gaps = []
        for i in range(n):
            layer_i = self._layers[i]
            layer_next = self._layers[(i + 1) % n]
            if i < n - 1:
                gap = (layer_i.z_centroid + layer_next.z_centroid) / 2.0
            else:
                # Wrap-around gap: midpoint with wrap adjustment
                gap = (layer_i.z_centroid + layer_next.z_centroid + 1.0) / 2.0
                gap = gap % 1.0
            gaps.append(gap)

        self._gaps = gaps
        return gaps

    def _validate_layer_count(self) -> None:
        """Validate that requested layers < available bulk layers."""
        if not self._layers:
            raise RuntimeError("Must call _identify_layers() before validation")

        total_layers = len(self._layers)
        if self.n_layers >= total_layers:
            raise ValueError(
                f"Requested {self.n_layers} layers but bulk only has "
                f"{total_layers} atomic layers after "
                f"({' '.join(str(d) for d in self.miller)}) transformation. "
                f"Reduce --layers or use a larger supercell."
            )

    def build_all(self, outdir=None):
        outdir = Path(outdir) if outdir else self.outdir
        self._validate_layer_count()
        return []
