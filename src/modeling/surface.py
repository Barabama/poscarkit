# src/modeling/surface.py

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import csv
import logging

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

        # 3. Build new cell: a, b unchanged; c = thickness + vacuum
        a_vec = xform_cell[0].copy()
        b_vec = xform_cell[1].copy()
        c_dir = xform_cell[2]
        c_dir_norm = c_dir / np.linalg.norm(c_dir)

        z_cart = np.array([a.coord[2] for a in all_atoms])
        slab_thickness = float(np.max(z_cart) - np.min(z_cart))
        c_length = slab_thickness + self.vacuum_bottom + self.vacuum_top
        c_vec = c_dir_norm * c_length

        new_cell = np.array([a_vec, b_vec, c_vec])

        # 4. Translate atoms: bottom of slab at z = vacuum_bottom
        z_min = float(np.min(z_cart))
        z_shift = self.vacuum_bottom - z_min
        for atom in all_atoms:
            atom.coord[2] += z_shift

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
                      from unique Cartesian z-values (fallback for backward
                      compatibility).
        """
        if self.fix_layers is not None:
            n_fix = self.fix_layers
        else:
            n_fix = (self.n_layers + 1) // 2

        coords = slab.get_coords(direct=False)
        z_vals = coords[:, 2]
        z_unique = np.sort(np.unique(np.round(z_vals, decimals=self.precision)))

        slab_layers = n_layers if n_layers is not None else max(1, len(z_unique))

        if n_fix > slab_layers:
            raise ValueError(
                f"Cannot fix {n_fix} layers in a {slab_layers}-layer slab. "
                f"Reduce --fix-layers."
            )

        if slab_layers <= 1:
            for atom in slab:
                atom.constr = ["T", "T", "T"]
            return slab

        fix_threshold = z_unique[min(n_fix, len(z_unique) - 1)] + 0.01
        fix_constr = ["T", "T", "F"] if self.fix_z_only else ["F", "F", "F"]
        free_constr = ["T", "T", "T"]

        for atom in slab:
            atom_z = atom.coord[2]
            if atom_z <= fix_threshold:
                atom.constr = fix_constr.copy()
            else:
                atom.constr = free_constr.copy()

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

        z_unique = np.sort(np.unique(np.round(z_vals, decimals=self.precision)))
        if len(z_unique) >= 2:
            top_threshold = z_unique[-2] + 0.01
            top_atoms = int(np.sum(z_vals >= top_threshold))
            surface_energy_est = top_atoms / max(total, 1)
        else:
            surface_energy_est = 1.0

        z_center = (np.max(z_vals) + np.min(z_vals)) / 2.0
        dipole = float(np.sum(z_vals - z_center))

        if len(z_unique) >= 2:
            top_threshold = z_unique[-2] + 0.01
            top_zs = z_vals[z_vals >= top_threshold]
            top_layer_z_std = float(np.std(top_zs)) if len(top_zs) > 0 else 0.0
        else:
            top_layer_z_std = 0.0

        z_rounded = np.round(z_vals, decimals=self.precision)
        z_unique_rounded = np.sort(np.unique(z_rounded))
        layer_comps_list = []
        for zu in z_unique_rounded:
            mask = np.abs(z_rounded - zu) < 0.5 * (10 ** (-self.precision))
            layer_syms = [
                slab[i].symbol for i, m in enumerate(mask) if m
            ]
            layer_count = Counter(layer_syms)
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
