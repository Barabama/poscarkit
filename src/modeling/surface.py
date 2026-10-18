# src/modeling/surface.py

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import csv
import logging

import numpy as np

from src.modeling.base import Atom, Struct, SimplePoscar
from src.modeling.layers import (
    DEFAULT_LAYER_GAP,
    cut_commensurate,
    layer_indices,
)


def _rotation_to_z(normal: np.ndarray) -> np.ndarray:
    """Return a rotation matrix that maps ``normal`` onto the +z axis."""
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    z_axis = np.array([0.0, 0.0, 1.0])

    cross = np.cross(normal, z_axis)
    sin = float(np.linalg.norm(cross))
    cos = float(np.dot(normal, z_axis))

    if sin < 1e-12:
        # Already aligned (identity) or anti-parallel (180 degree flip).
        return np.eye(3) if cos > 0 else np.diag([1.0, -1.0, -1.0])

    skew = np.array([
        [0.0, -cross[2], cross[1]],
        [cross[2], 0.0, -cross[0]],
        [-cross[1], cross[0], 0.0],
    ])
    return np.eye(3) + skew + skew @ skew * ((1.0 - cos) / (sin ** 2))


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
        precision: Deprecated and unused. Layers are identified with ``gap``.
        gap: Minimum inter-layer gap in Angstrom used to identify layers.
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
        gap: float = DEFAULT_LAYER_GAP,
        name: str = "surface",
    ):
        self.poscar = Path(poscar)
        self.name = name
        self.miller = miller
        self.n_layers = layers
        self.vacuum = vacuum
        self.fix_layers = fix_layers
        self.fix_z_only = fix_z_only
        self.outdir = Path(outdir) if outdir else Path("output")
        self.precision = precision
        self.gap = gap

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

        transformed_atoms, tolerance = cut_commensurate(
            atoms, np.array([a_dir, b_dir, c_dir], dtype=float)
        )
        logging.info(
            f"Cut with basis {[a_dir, b_dir, c_dir]} "
            f"(tolerance={tolerance}, atoms={len(transformed_atoms)})"
        )
        self._transformed = SimplePoscar.atoms2struct(transformed_atoms)
        return self._transformed

    def _identify_layers(self) -> list[Layer]:
        """Identify atomic layers in the transformed bulk.

        Layers are found by clustering atoms along the surface normal and
        splitting at real gaps (see :func:`src.modeling.layers.layer_indices`).
        This keeps a thermally broadened layer (~0.07 A wide) together, unlike
        rounding to a fixed number of decimal places, which shreds such a layer
        into several fragments.

        Returns list of Layer objects ordered by fractional z ascending.
        """
        if self._transformed is None:
            raise RuntimeError("Must call _transform_cell() before _identify_layers()")

        struct = self._transformed
        # Compute Cartesian positions without mutating the struct: _build_slab
        # relies on the stored coordinates being direct.
        coords_direct = struct.get_coords(direct=True)
        coords_cart = np.dot(coords_direct, struct.cell)

        if len(coords_direct) == 0:
            return []

        groups = layer_indices(coords_cart, struct.cell, gap=self.gap)

        # Build Layer objects from groups
        result = []
        for li, (z_centroid, indices) in enumerate(groups):
            atoms_in_layer = [struct[i] for i in indices]

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

    def _unwrap_zs(self, z_frac: np.ndarray, gap_z: float) -> np.ndarray:
        """Unwrap fractional z-coordinates for continuity below gap_z.

        Starting from gap_z going downward, adjust z by +/- 1 whenever
        the sequence would cross the periodic boundary at z=0 or z=1.
        Returns unwrapped fractional z values (monotonically decreasing
        from gap_z downward).
        """
        n = len(z_frac)

        # Initial offset: atoms with z > gap_z are "above" in direct space
        # and need -1 to move below the gap
        offsets = np.where(z_frac > gap_z, -1.0, 0.0)
        adjusted_z = z_frac + offsets

        # Sort descending from gap (highest adjusted z = closest to gap)
        sort_idx = np.argsort(adjusted_z)[::-1]

        unwrapped = np.zeros(n)
        for rank, idx in enumerate(sort_idx):
            if rank == 0:
                unwrapped[idx] = adjusted_z[idx]
            else:
                prev_idx = sort_idx[rank - 1]
                extra = 0.0
                if adjusted_z[idx] > adjusted_z[prev_idx]:
                    extra = -1.0
                unwrapped[idx] = adjusted_z[idx] + extra

        return unwrapped

    def _build_slab(self, gap_z: float, n_layers: int) -> Struct:
        """Build a single slab from a gap with n_layers below it."""
        layers = self._layers
        total_layers = len(layers)

        # Find the layer just below the gap
        gap_below_idx = None
        min_dist = float('inf')
        for i, layer in enumerate(layers):
            dist = (gap_z - layer.z_centroid) % 1.0
            if dist < min_dist:
                min_dist = dist
                gap_below_idx = i

        # Select n_layers going downward from gap_below_idx
        selected_layers = []
        for offset in range(n_layers):
            layer_idx = (gap_below_idx - offset) % total_layers
            selected_layers.append(layers[layer_idx])

        # Collect atoms with fractional coords from transformed cell
        all_atoms = []
        for layer in selected_layers:
            for atom in layer.atoms:
                all_atoms.append(Atom(
                    index=len(all_atoms),
                    symbol=atom.symbol,
                    coord=atom.coord.copy(),
                    constr=atom.constr.copy() if atom.constr else [],
                    note=atom.note,
                    meta=atom.meta,
                ))

        if self._transformed is None:
            raise RuntimeError("Must call _transform_cell() before _build_slab()")
        xform_cell = self._transformed.cell

        # 1. Unwrap fractional z-coordinates
        z_frac = np.array([a.coord[2] for a in all_atoms])
        z_unwrapped = self._unwrap_zs(z_frac, gap_z)

        # 2. Replace fractional z with unwrapped value, convert to Cartesian
        for atom, zw in zip(all_atoms, z_unwrapped):
            frac = atom.coord.copy()
            frac[2] = zw
            atom.coord = np.dot(frac, xform_cell)

        # 3. Rotate so the stacking direction points along Cartesian z, then
        #    build the new cell: a, b unchanged in-plane, c = thickness + vacuum.
        #    For (110) the transformed c axis has no z component at all, so
        #    without this rotation the slab would not lie along its own c axis.
        c_dir = xform_cell[2]
        c_dir_norm = c_dir / np.linalg.norm(c_dir)
        rotation = _rotation_to_z(c_dir_norm)

        proj = np.array([float(np.dot(a.coord, c_dir_norm)) for a in all_atoms])
        slab_thickness = float(np.max(proj) - np.min(proj))
        c_length = slab_thickness + self.vacuum_bottom + self.vacuum_top

        a_vec = rotation @ xform_cell[0]
        b_vec = rotation @ xform_cell[1]
        c_vec = np.array([0.0, 0.0, c_length])
        new_cell = np.array([a_vec, b_vec, c_vec])

        # 4. Rotate atoms and translate so the bottom of the slab sits at
        #    z = vacuum_bottom.
        z_shift = self.vacuum_bottom - float(np.min(proj))
        for atom in all_atoms:
            rotated = rotation @ atom.coord
            rotated[2] += z_shift
            atom.coord = rotated

        slab_struct = Struct(cell=new_cell, is_direct=False, atom_list=all_atoms)
        slab_struct.get_coords(direct=True)

        # 5. Apply constraints
        slab_struct = self._add_constraints(slab_struct, n_layers)

        return slab_struct

    def _add_constraints(self, slab: Struct, n_layers: int | None = None) -> Struct:
        """Apply Selective Dynamics to the slab.

        Bottom N layers are fixed; top layers are free.
        N determined by auto-fix rule or user override.
        Raises ValueError if fix_layers >= total slab layers.

        Parameters:
            slab: The slab to apply constraints to.
            n_layers: Number of atomic layers in the slab. If None, detected
                      from the gap-based layer clustering.
        """
        if self.fix_layers is not None:
            n_fix = self.fix_layers
        else:
            n_fix = (self.n_layers + 1) // 2

        coords = slab.get_coords(direct=False)
        groups = layer_indices(coords, slab.cell, gap=self.gap)

        slab_layers = n_layers if n_layers is not None else max(1, len(groups))

        if n_fix > slab_layers:
            raise ValueError(
                f"Cannot fix {n_fix} layers in a {slab_layers}-layer slab. "
                f"Reduce --fix-layers."
            )

        if slab_layers <= 1:
            for atom in slab:
                atom.constr = ["T", "T", "T"]
            return slab

        fix_constr = ["T", "T", "F"] if self.fix_z_only else ["F", "F", "F"]
        free_constr = ["T", "T", "T"]

        for atom in slab:
            atom.constr = free_constr.copy()

        # Fix the lowest layers as whole layers: a relaxed layer is a
        # continuous distribution, so a z threshold would cut through it.
        n_fix_groups = min(n_fix, len(groups) - 1) + 1
        for _, indices in groups[:n_fix_groups]:
            for idx in indices:
                slab[int(idx)].constr = fix_constr.copy()

        return slab

    def build_all(self, outdir: Path | None = None) -> list[Path]:
        """Build all possible slabs from all gaps, with N and N+1 layers.

        Returns list of output file paths.
        """
        outdir = Path(outdir) if outdir else self.outdir
        self._validate_layer_count()

        all_slabs = []

        name = self.name
        miller_str = "".join(str(d) for d in self.miller)

        subdir = outdir / f"{name}-slab-{miller_str}"
        subdir.mkdir(parents=True, exist_ok=True)

        for gap_idx, gap_z in enumerate(self._gaps):
            # Find the layer just below this gap
            gap_below_idx = None
            min_dist = float('inf')
            for i, layer in enumerate(self._layers):
                dist = (gap_z - layer.z_centroid) % 1.0
                if dist < min_dist:
                    min_dist = dist
                    gap_below_idx = i

            # Filter: gap must have enough distinct layers below without wrapping
            if gap_below_idx is not None and gap_below_idx < self.n_layers:
                continue

            term_id = f"term{gap_idx + 1:02d}"

            for nl in (self.n_layers, self.n_layers + 1):
                slab = self._build_slab(gap_z, nl)
                filename = (
                    f"{name}-slab-{miller_str}-{term_id}-layers{nl}.vasp"
                )
                slab_path = subdir / filename
                SimplePoscar.write_poscar(
                    poscar=slab_path, struct=slab,
                    comment=f"Slab {miller_str} {term_id} {nl} layers"
                )
                all_slabs.append((slab_path, slab, gap_z, nl, term_id))

        # Write summary CSV (placeholder, full impl in Task 6)
        self._write_summary(all_slabs, subdir)

        return [s[0] for s in all_slabs]

    def _compute_metrics(self, slab: Struct, n_layers: int, term_id: str) -> dict:
        """Compute all metrics for a slab."""
        coords = slab.get_coords(direct=False)
        z_vals = coords[:, 2]

        from collections import Counter
        sym_counts = Counter(a.symbol for a in slab)
        total = len(slab)
        composition = "".join(
            f"{s}{c}" if c > 1 else s
            for s, c in sorted(sym_counts.items())
        )

        # Bulk composition for deviation
        bulk_counts = Counter(a.symbol for a in self._bulk_struct)
        bulk_total = len(self._bulk_struct)

        all_elements = set(sym_counts.keys()) | set(bulk_counts.keys())
        deviations = []
        for elem in sorted(all_elements):
            surf_frac = sym_counts.get(elem, 0) / total
            bulk_frac = bulk_counts.get(elem, 0) / bulk_total
            deviations.append(abs(surf_frac - bulk_frac))
        composition_deviation = sum(deviations) / len(deviations) if deviations else 0.0

        # Identify layers by clustering, not by rounding: a relaxed layer is a
        # continuous distribution and would be split by a decimal grid.
        groups = layer_indices(coords, slab.cell, gap=self.gap)
        if groups:
            _, top_indices = groups[-1]
            top_atoms = len(top_indices)
            surface_energy_est = top_atoms / max(total, 1)
            top_layer_z_std = float(np.std(z_vals[top_indices]))
        else:
            surface_energy_est = 1.0
            top_layer_z_std = 0.0

        z_center = (np.max(z_vals) + np.min(z_vals)) / 2.0
        dipole = float(np.sum(z_vals - z_center))

        layer_comps_list = []
        for _, indices in groups:
            layer_count = Counter(slab[int(i)].symbol for i in indices)
            layer_comps_list.append("".join(
                f"{s}{c}" if c > 1 else s
                for s, c in sorted(layer_count.items())
            ))
        layer_compositions = " | ".join(layer_comps_list)

        fix_mode = "TTF" if self.fix_z_only else "FFF"
        fix_n = self.fix_layers if self.fix_layers is not None else ((self.n_layers + 1) // 2)

        return {
            "n_layers": n_layers,
            "n_atoms": total,
            "composition": composition,
            "layer_compositions": layer_compositions,
            "composition_deviation": round(composition_deviation, 6),
            "surface_energy_est": round(surface_energy_est, 6),
            "dipole_moment": round(dipole, 6),
            "top_layer_z_std": round(top_layer_z_std, 6),
            "fix_layers": fix_n,
            "fix_mode": fix_mode,
            "vacuum_top": round(self.vacuum_top, 3),
            "vacuum_bottom": round(self.vacuum_bottom, 3),
        }

    def _write_summary(
        self,
        all_slabs: list,
        outdir: Path,
    ) -> Path:
        """Write summary CSV for all generated slabs."""
        csv_path = outdir / "summary.csv"
        fieldnames = [
            "filename", "termination_id", "gap_position_z",
            "n_layers", "n_atoms", "composition", "layer_compositions",
            "composition_deviation", "surface_energy_est", "dipole_moment",
            "top_layer_z_std", "fix_layers", "fix_mode",
            "vacuum_top", "vacuum_bottom",
        ]

        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for slab_path, slab, gap_z, n_layers, term_id in all_slabs:
                metrics = self._compute_metrics(slab, n_layers, term_id)
                row = {
                    "filename": slab_path.name,
                    "termination_id": term_id,
                    "gap_position_z": round(float(gap_z), 6),
                }
                row.update(metrics)
                writer.writerow(row)

        logging.info(f"Summary saved to {csv_path}")
        return csv_path
