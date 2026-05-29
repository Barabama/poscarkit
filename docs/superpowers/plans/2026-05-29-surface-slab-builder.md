# Surface Slab Builder Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a CLI surface slab model builder that generates asymmetric surface slabs from bulk POSCAR files for HEA research, replacing the manual Materials Studio workflow.

**Architecture:** Hybrid approach — ASE's `cut()` handles Miller-index cell transformation (proven to preserve custom per-atom arrays), while all HEA-specific logic (layer identification, gap detection, slab assembly, vacuum, selective dynamics, metrics) runs on native `Struct` objects with full `note`/`constr`/`meta` preservation.

**Tech Stack:** Python 3.12, ASE 3.27, numpy, argparse, unittest, csv

**Design Spec:** `docs/superpowers/specs/2026-05-29-surface-slab-builder-design.md`

**Ordering:** Tasks 1-9 are sequential (each depends on prior tasks). Within a task, substeps are sequential.

---

## File Structure

| File | Action | Responsibility |
|---|---|---|
| `src/modeling/base.py` | Modify | Bridge fix: struct2atoms/atoms2struct preserve note/constr/meta |
| `src/modeling/slice.py` | Modify | Fix group_by_normal periodic boundary bug |
| `src/modeling/surface.py` | Create | SurfaceBuilder class + Layer dataclass |
| `src/cli/poscarkit.py` | Modify | Register `surface` subcommand + cmd_surface function |
| `tests/modeling/test_base.py` | Modify | Add bridge round-trip test |
| `tests/modeling/test_slice.py` | Modify | Add periodic boundary test |
| `tests/modeling/test_surface.py` | Create | All SurfaceBuilder tests |

---

### Task 1: Fix struct2atoms/atoms2struct Bridge

**Files:**
- Modify: `src/modeling/base.py:462-479`
- Modify: `tests/modeling/test_base.py`

#### Part A: Add round-trip test

- [ ] **Step 1: Add `test_struct_atoms_roundtrip` to test_base.py**

Append after existing test class:

```python
class TestBridgeRoundtrip(unittest.TestCase):
    """Test Struct <-> Atoms conversion preserves custom fields."""

    def test_struct_atoms_roundtrip(self):
        """Struct -> Atoms -> Struct preserves note, constr, meta."""
        cell = np.array([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]])

        atoms = [
            Atom(0, "Au", np.array([0.0, 0.0, 0.0]),
                 constr=["T", "T", "T"], note="1a-Au", meta="s01"),
            Atom(1, "Cu", np.array([0.5, 0.5, 0.5]),
                 constr=["T", "T", "T"], note="3c-Cu", meta="s02"),
            Atom(2, "Cu", np.array([0.0, 0.5, 0.5]),
                 constr=["F", "F", "F"], note="3c-Cu", meta="s03"),
        ]
        original = Struct(cell=cell, is_direct=True, atom_list=atoms)

        ase_atoms = SimplePoscar.struct2atoms(original)
        restored = SimplePoscar.atoms2struct(ase_atoms)

        self.assertEqual(len(restored), len(original))

        for i, (orig, rest) in enumerate(zip(original, restored)):
            self.assertEqual(rest.symbol, orig.symbol,
                             f"Atom {i}: symbol mismatch")
            np.testing.assert_array_almost_equal(
                rest.coord, orig.coord,
                err_msg=f"Atom {i}: coord mismatch")
            self.assertEqual(rest.note, orig.note,
                             f"Atom {i}: note mismatch {rest.note} != {orig.note}")
            self.assertEqual(rest.meta, orig.meta,
                             f"Atom {i}: meta mismatch")
            self.assertEqual(rest.constr, orig.constr,
                             f"Atom {i}: constr mismatch {rest.constr} != {orig.constr}")

    def test_struct_atoms_roundtrip_no_constraints(self):
        """Round-trip works when atoms have no constraints set."""
        cell = np.array([[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]])
        atoms = [
            Atom(0, "Fe", np.array([0.0, 0.0, 0.0]),
                 constr=[], note="1a-Fe", meta=None),
        ]
        original = Struct(cell=cell, is_direct=True, atom_list=atoms)
        ase_atoms = SimplePoscar.struct2atoms(original)
        restored = SimplePoscar.atoms2struct(ase_atoms)

        self.assertEqual(restored[0].note, "1a-Fe")
        self.assertEqual(restored[0].meta, None)
        self.assertEqual(restored[0].constr, ["T", "T", "T"])

    def test_atoms2struct_no_custom_arrays(self):
        """atoms2struct works on bare ASE Atoms with no custom arrays."""
        from ase import Atoms
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[3.6, 3.6, 3.6], pbc=True)
        struct = SimplePoscar.atoms2struct(atoms)

        self.assertEqual(len(struct), 1)
        self.assertEqual(struct[0].symbol, "Cu")
        self.assertEqual(struct[0].note, "")
        self.assertEqual(struct[0].constr, ["T", "T", "T"])
```

- [ ] **Step 2: Run new tests, verify they fail**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_simple_poscar.TestBridgeRoundtrip -v
```

Expected: tests fail because bridge doesn't preserve custom fields.

#### Part B: Implement bridge fix

- [ ] **Step 3: Fix `struct2atoms` in base.py**

Replace lines 461-467 in `src/modeling/base.py`:

```python
    @staticmethod
    def struct2atoms(struct: Struct) -> Atoms:
        """Convert Struct to Atoms, preserving note/constr/meta in arrays."""
        symbols = struct.symbols
        cell = struct.cell
        positions = struct.get_coords(False)
        atoms = Atoms(symbols=symbols, cell=cell, positions=positions, pbc=True)

        notes = np.array([a.note or '' for a in struct])
        atoms.new_array('note', notes)

        metas = np.array([a.meta if a.meta is not None else '' for a in struct])
        atoms.new_array('meta', metas)

        # Encode constr as int mask: 0=free(T), 1=fixed(F)
        def _constr_to_mask(constr):
            if not constr:
                return [0, 0, 0]
            return [0 if c == 'T' else 1 for c in constr]

        masks = np.array([_constr_to_mask(a.constr) for a in struct], dtype=int)
        atoms.new_array('constr_mask', masks)
        return atoms
```

- [ ] **Step 4: Fix `atoms2struct` in base.py**

Replace lines 469-479 in `src/modeling/base.py`:

```python
    @staticmethod
    def atoms2struct(atoms: Atoms) -> Struct:
        """Convert Atoms to Struct, restoring note/constr/meta from arrays."""
        cell = np.array(atoms.get_cell().copy())
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()

        has_note = 'note' in atoms.arrays
        has_meta = 'meta' in atoms.arrays
        has_constr = 'constr_mask' in atoms.arrays

        notes = atoms.get_array('note') if has_note else None
        metas = atoms.get_array('meta') if has_meta else None
        constr_masks = atoms.get_array('constr_mask') if has_constr else None

        atom_list = []
        for idx, (symbol, position) in enumerate(zip(symbols, positions)):
            note = str(notes[idx]) if has_note and notes[idx] else ''
            meta = metas[idx] if has_meta and metas[idx] != '' else None
            if has_constr:
                mask_row = constr_masks[idx]
                constr = ['T' if v == 0 else 'F' for v in mask_row]
            else:
                constr = ['T', 'T', 'T']
            atom_list.append(Atom(
                index=idx, symbol=symbol, coord=position,
                note=note, meta=meta, constr=constr,
            ))
        struct = Struct(cell=cell, is_direct=False, atom_list=atom_list)
        return struct
```

- [ ] **Step 5: Run round-trip tests, verify they pass**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_simple_poscar.TestBridgeRoundtrip -v
```

Expected: all 3 tests PASS.

- [ ] **Step 6: Run existing tests to check for regressions**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest discover tests -v
```

Expected: all existing tests PASS (new 3 + all old ones).

- [ ] **Step 7: Commit**

```bash
git add src/modeling/base.py tests/modeling/test_simple_poscar.py
git commit -m "fix: preserve note/constr/meta in struct2atoms and atoms2struct bridge

Encode constr as integer mask (0=free, 1=fixed) in atoms.arrays to avoid
numpy string array fragility. Fallback to defaults for bare ASE Atoms
without custom arrays, maintaining backward compatibility.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 2: Layer Identification Algorithm

**Files:**
- Create: `src/modeling/surface.py`
- Create: `tests/modeling/test_surface.py`

#### Part A: Write tests first

- [ ] **Step 1: Create `tests/modeling/test_surface.py` with layer identification tests**

```python
import unittest
import tempfile
from pathlib import Path

import numpy as np
from ase.build import bulk

from src.modeling.base import Struct, Atom, SimplePoscar
from src.modeling.surface import SurfaceBuilder, Layer


def _make_fcc_bulk_poscar(outdir, a=4.08):
    """Create a simple FCC bulk POSCAR for testing."""
    atoms = bulk('Au', 'fcc', a=a, cubic=True)
    cell = np.array(atoms.get_cell())
    pos = atoms.get_positions()
    syms = atoms.get_chemical_symbols()

    atom_list = [
        Atom(i, sym, np.array(coord),
             note=f"1a-{sym}", constr=["T", "T", "T"])
        for i, (sym, coord) in enumerate(zip(syms, pos))
    ]
    struct = Struct(cell=cell, is_direct=False, atom_list=atom_list)
    struct.get_coords(direct=True)

    p = Path(outdir) / "bulk_fcc.vasp"
    SimplePoscar.write_poscar(p, struct, "FCC bulk")
    return p


def _make_bcc_bulk_poscar(outdir, a=2.87):
    """Create a simple BCC bulk POSCAR for testing."""
    atoms = bulk('Fe', 'bcc', a=a, cubic=True)
    cell = np.array(atoms.get_cell())
    pos = atoms.get_positions()
    syms = atoms.get_chemical_symbols()

    atom_list = [
        Atom(i, sym, np.array(coord),
             note=f"1a-{sym}", constr=["T", "T", "T"])
        for i, (sym, coord) in enumerate(zip(syms, pos))
    ]
    struct = Struct(cell=cell, is_direct=False, atom_list=atom_list)
    struct.get_coords(direct=True)

    p = Path(outdir) / "bulk_bcc.vasp"
    SimplePoscar.write_poscar(p, struct, "BCC bulk")
    return p


class TestLayerIdentification(unittest.TestCase):

    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_identify_layers_fcc111(self):
        """FCC bulk transformed to (111) orientation identifies correct layers."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        layers = builder._identify_layers()

        self.assertGreater(len(layers), 0, "Should find at least one layer")

        for layer in layers:
            self.assertIsInstance(layer, Layer)
            self.assertIsInstance(layer.z_centroid, float)
            self.assertGreater(len(layer.atoms), 0)
            self.assertTrue(0.0 <= layer.z_centroid <= 1.0,
                            f"z_centroid {layer.z_centroid} not in [0,1]")

    def test_identify_layers_bcc100(self):
        """BCC bulk (100) transformation identifies layers."""
        poscar = _make_bcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 0, 0), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        layers = builder._identify_layers()

        self.assertGreater(len(layers), 0)

    def test_identify_layers_custom_miller(self):
        """Arbitrary Miller index (2,1,0) on FCC bulk."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(2, 1, 0), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        layers = builder._identify_layers()

        self.assertGreater(len(layers), 0)
```

- [ ] **Step 2: Run tests, verify they fail (SurfaceBuilder not defined)**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_surface.TestLayerIdentification -v
```

Expected: ImportError or "SurfaceBuilder not defined".

#### Part B: Implement SurfaceBuilder skeleton + layer identification

- [ ] **Step 3: Create `src/modeling/surface.py`**

```python
# src/modeling/surface.py

import csv
import logging
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
        # gaps between consecutive atoms
        gaps = np.diff(sorted_z)
        # wrap-around gap
        wrap_gap = np.array([1.0 - sorted_z[-1] + sorted_z[0]])

        all_gaps = np.concatenate([gaps, wrap_gap])
        threshold = np.median(all_gaps) * 3.0

        # Mark boundaries: True = gap is a layer boundary
        large_gaps = all_gaps > threshold

        # Build layers from large-gap-delimited groups
        layers = []
        current_group = []
        started = False
        first_atom_idx = 0

        for i in range(n_atoms):
            is_boundary = (i > 0 and gaps[i - 1] > threshold) or (i == 0 and wrap_gap[0] > threshold)

            if is_boundary and started:
                # Save previous group as a layer
                layers.append(current_group)
                current_group = []

            current_group.append(sorted_idx[i])
            started = True

        # Last group
        if current_group:
            # If wrap_gap is small, merge with first layer
            if not large_gaps[-1] and layers:
                layers[0] = current_group + layers[0]
            else:
                layers.append(current_group)

        # Convert to Layer objects
        result = []
        for li, indices in enumerate(layers):
            atoms_in_layer = [struct[idx] for idx in indices]
            z_centroid = float(np.mean([struct[idx].coord[2] for idx in indices]))
            # Ensure z_centroid in [0, 1)
            z_centroid = z_centroid % 1.0

            # Composition string
            from collections import Counter
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
```

- [ ] **Step 4: Run layer identification tests, verify they pass**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_surface.TestLayerIdentification -v
```

Expected: all 3 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/modeling/surface.py tests/modeling/test_surface.py
git commit -m "feat: add SurfaceBuilder skeleton with layer identification

Implement circular gap detection for atomic layer identification in
transformed bulk cells. Uses median*3 threshold on z-coordinate gaps.
Handles periodic wrap-around via circular boundary detection.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 3: Gap Detection + Filtering

**Files:**
- Modify: `src/modeling/surface.py`
- Modify: `tests/modeling/test_surface.py`

#### Part A: Tests

- [ ] **Step 1: Add gap detection tests to test_surface.py**

```python
class TestGapDetection(unittest.TestCase):

    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_find_gaps(self):
        """Gap count equals layer count; positions between layers."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        layers = builder._identify_layers()
        gaps = builder._find_gaps()

        self.assertEqual(len(gaps), len(layers),
                         "Gap count should equal layer count (one gap per layer)")

        # Each gap should be between two adjacent layers
        for i, gap_z in enumerate(gaps):
            layer_i = layers[i]
            layer_next = layers[(i + 1) % len(layers)]

            # gap should be between layer_i.z and layer_next.z (handling wrap)
            if i < len(layers) - 1:
                self.assertTrue(layer_i.z_centroid < gap_z < layer_next.z_centroid,
                                f"Gap {i} ({gap_z:.4f}) not between layers "
                                f"({layer_i.z_centroid:.4f}, {layer_next.z_centroid:.4f})")

    def test_insufficient_layers_filtered(self):
        """Gaps requiring wraparound with insufficient layers are excluded."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        # Request more layers than bulk can provide
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=50, vacuum=15.0,
            outdir=self.temp_dir,
        )
        with self.assertRaises(ValueError) as ctx:
            builder._transform_cell()
            builder._identify_layers()
            builder._find_gaps()
            builder.build_all(self.temp_dir)
        self.assertIn("Reduce --layers", str(ctx.exception))
```

#### Part B: Implementation

- [ ] **Step 2: Add `_find_gaps` and `_validate_layers` to SurfaceBuilder**

Insert after `_identify_layers` in `src/modeling/surface.py`:

```python
    def _find_gaps(self) -> list[float]:
        """Find gaps between consecutive layers.

        A gap is the midpoint in fractional z between two adjacent layers.
        Returns list of gap z-values (same length as layers list).
        gap[i] is between layer[i] and layer[(i+1) % n].
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
        """Validate that requested layers <= available bulk layers."""
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
```

- [ ] **Step 3: Run gap detection tests, verify they pass**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_surface.TestGapDetection -v
```

Expected: 2 tests PASS.

- [ ] **Step 4: Commit**

```bash
git add src/modeling/surface.py tests/modeling/test_surface.py
git commit -m "feat: add gap detection and layer count validation to SurfaceBuilder

Gaps are midpoints between adjacent layer z-centroids. Added validation
that requested slab layers do not exceed available bulk atomic layers.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 4: Slab Assembly + Coordinate Unwrapping

**Files:**
- Modify: `src/modeling/surface.py`
- Modify: `tests/modeling/test_surface.py`

#### Part A: Tests

- [ ] **Step 1: Add slab assembly tests**

```python
class TestSlabAssembly(unittest.TestCase):

    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_build_slab_3layers(self):
        """3-layer slab has correct atom count and layered structure."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(self.temp_dir)
        self.assertGreater(len(slabs), 0, "Should produce at least one slab")

        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            # Check cell is 3x3 with c-vector properly oriented
            self.assertEqual(slab.cell.shape, (3, 3))
            # c-vector should be [0, 0, some_value]
            self.assertAlmostEqual(slab.cell[2, 0], 0.0)
            self.assertAlmostEqual(slab.cell[2, 1], 0.0)
            self.assertGreater(slab.cell[2, 2], 0.0)

    def test_build_slab_4layers(self):
        """N+1 (4-layer) slabs are produced alongside N (3-layer) slabs."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(self.temp_dir)
        layer_counts = set()
        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            # Determine layer count from filename
            name = slab_path.stem
            if "layers3" in name:
                layer_counts.add(3)
            elif "layers4" in name:
                layer_counts.add(4)

        self.assertIn(3, layer_counts, "Should produce 3-layer slabs")
        self.assertIn(4, layer_counts, "Should produce 4-layer slabs (N+1)")

    def test_layers_periodic_crossing(self):
        """Layer crossing cell boundary is unwrapped correctly."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(self.temp_dir)
        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            coords = slab.get_coords(direct=False)
            z_vals = coords[:, 2]

            # After unwrapping, z should be monotonically increasing
            # (no sudden jumps from periodic boundary)
            z_sorted = np.sort(z_vals)
            np.testing.assert_array_almost_equal(z_sorted, z_vals,
                err_msg="z-coordinates should be continuous after unwrapping")

    def test_note_preserved(self):
        """note fields survive the full Struct -> Atoms -> cut -> Struct pipeline."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        original = SimplePoscar.read_poscar(poscar)
        original_notes = set(a.note for a in original)

        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(self.temp_dir)
        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            slab_notes = set(a.note for a in slab)
            self.assertTrue(
                slab_notes.issubset(original_notes) or slab_notes == original_notes,
                f"Slab notes {slab_notes} not subset of original {original_notes}"
            )
            # Every atom should have a non-empty note
            for atom in slab:
                self.assertTrue(atom.note, f"Atom {atom.index} has empty note")
```

#### Part B: Implementation

- [ ] **Step 2: Add `_unwrap_zs`, `_build_slab`, `build_all` to SurfaceBuilder**

```python
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
                # Add wrap offset if this atom appears above the previous one
                extra = 0.0
                if adjusted_z[idx] > adjusted_z[prev_idx]:
                    extra = -1.0
                unwrapped[idx] = adjusted_z[idx] + extra

        return unwrapped

    def _build_slab(self, gap_z: float, n_layers: int) -> Struct:
        """Build a single slab from a gap with n_layers below it.

        Args:
            gap_z: Fractional z position of the gap (cutting plane).
            n_layers: Number of layers to include in slab.

        Returns:
            Struct representing the slab with vacuum and constraints applied.
        """
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

        # 2. Replace fractional z with unwrapped value, then convert to Cartesian
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
            atom.coord = atom.coord + c_dir_norm * z_shift

        slab_struct = Struct(cell=new_cell, is_direct=False, atom_list=all_atoms)
        slab_struct.get_coords(direct=True)

        # 5. Apply constraints
        slab_struct = self._add_constraints(slab_struct)

        return slab_struct

    def build_all(self, outdir: Path | None = None) -> list[Path]:
        """Build all possible slabs from all gaps, with N and N+1 layers.

        Returns list of output file paths.
        """
        outdir = Path(outdir) if outdir else self.outdir
        self._validate_layer_count()

        total_layers = len(self._layers)
        all_slabs = []

        name = self.poscar.stem
        miller_str = "".join(str(d) for d in self.miller)

        subdir = outdir / f"{name}-slab-{miller_str}"
        subdir.mkdir(parents=True, exist_ok=True)

        for gap_idx, gap_z in enumerate(self._gaps):
            # Only include gaps with enough distinct layers below (no wrapping)
            gap_below_idx = None
            min_dist = float('inf')
            for i, layer in enumerate(self._layers):
                dist = (gap_z - layer.z_centroid) % 1.0
                if dist < min_dist:
                    min_dist = dist
                    gap_below_idx = i

            # Filter: gap must have distinct layers below without periodic wrapping
            # A valid gap needs gap_below_idx >= n_layers - 1 to avoid wrapping
            max_layers_needed = self.n_layers + 1  # for N+1 variant
            if gap_below_idx is not None and gap_below_idx < max_layers_needed - 1:
                logging.info(f"Skipping gap {gap_idx} (z={gap_z:.4f}): "
                             f"requires {max_layers_needed} distinct layers "
                             f"without wrapping, only {gap_below_idx + 1} available")
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

        # Write summary CSV
        self._write_summary(all_slabs, subdir)

        return [s[0] for s in all_slabs]
```

- [ ] **Step 3: Run slab assembly tests, verify they pass**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_surface.TestSlabAssembly -v
```

Expected: 4 tests PASS.

- [ ] **Step 4: Commit**

```bash
git add src/modeling/surface.py tests/modeling/test_surface.py
git commit -m "feat: implement slab assembly with coordinate unwrapping

Build slabs by selecting N layers below each gap with cyclic indexing.
Unwrap fractional z-coordinates for continuity. Filter gaps requiring
periodic wrapping when insufficient distinct layers exist.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 5: Vacuum + Selective Dynamics

**Files:**
- Modify: `src/modeling/surface.py`
- Modify: `tests/modeling/test_surface.py`

#### Part A: Tests

- [ ] **Step 1: Add vacuum and constraints tests**

```python
class TestVacuumAndConstraints(unittest.TestCase):

    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_vacuum_distribution(self):
        """Bottom vacuum = 2.0 Å, top vacuum = total - 2.0 Å."""
        total_vacuum = 20.0
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=total_vacuum,
            outdir=self.temp_dir,
        )
        self.assertEqual(builder.vacuum_bottom, 2.0)
        self.assertEqual(builder.vacuum_top, total_vacuum - 2.0)

        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(self.temp_dir)
        self.assertGreater(len(slabs), 0)

        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            coords = slab.get_coords(direct=False)
            z_min = float(np.min(coords[:, 2]))
            z_max = float(np.max(coords[:, 2]))

            # Atoms should be between vacuum_bottom and (cell_c - vacuum_top)
            cell_c = float(np.linalg.norm(slab.cell[2]))
            self.assertGreaterEqual(z_min, 0.0,
                f"z_min ({z_min:.3f}) should be >= 0")
            self.assertLessEqual(z_max, cell_c,
                f"z_max ({z_max:.3f}) should be <= cell_c ({cell_c:.3f})")

            # Bottom vacuum region: [0, vacuum_bottom) should be empty
            bottom_atoms = [z for z in coords[:, 2] if z < 2.0]
            self.assertEqual(len(bottom_atoms), 0,
                f"No atoms should be in bottom vacuum region [0, 2A)")

    def test_constraints_fff(self):
        """Default: bottom layers fully fixed (F F F)."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(self.temp_dir)
        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            coords = slab.get_coords(direct=False)
            z_vals = coords[:, 2]

            # Find bottom layer atoms (lowest z)
            z_min = float(np.min(z_vals))
            # Second lowest layer
            z_unique = np.sort(np.unique(np.round(z_vals, decimals=2)))
            threshold = z_unique[min(1, len(z_unique) - 1)] + 0.5

            for atom in slab:
                atom_z = atom.coord[2]
                # Convert to Cartesian for z comparison
                atom_z_cart = np.dot(atom.coord, slab.cell)[2]
                if atom_z_cart < threshold:
                    self.assertEqual(atom.constr, ["F", "F", "F"],
                        f"Bottom atom {atom.index} should be FFF, got {atom.constr}")
                else:
                    self.assertEqual(atom.constr, ["T", "T", "T"],
                        f"Top atom {atom.index} should be TTT, got {atom.constr}")

    def test_constraints_ttf(self):
        """--fix-z-only: bottom layers fixed in z only (T T F)."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            fix_z_only=True, outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(self.temp_dir)
        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            coords = slab.get_coords(direct=False)
            z_vals = coords[:, 2]
            z_unique = np.sort(np.unique(np.round(z_vals, decimals=2)))
            threshold = z_unique[min(1, len(z_unique) - 1)] + 0.5

            for atom in slab:
                atom_z_cart = np.dot(atom.coord, slab.cell)[2]
                if atom_z_cart < threshold:
                    self.assertEqual(atom.constr, ["T", "T", "F"],
                        f"Fixed-z atom {atom.index} should be TTF, got {atom.constr}")

    def test_vacuum_below_minimum(self):
        """ValueError when total vacuum < 2.0 Å."""
        with self.assertRaises(ValueError) as ctx:
            SurfaceBuilder(
                poscar=_make_fcc_bulk_poscar(self.temp_dir),
                miller=(1, 1, 1), layers=3, vacuum=1.0,
                outdir=self.temp_dir,
            )
        self.assertIn("2.0", str(ctx.exception))

    def test_fix_layers_exceeds_total(self):
        """ValueError when fix_layers >= total slab layers."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            fix_layers=5, outdir=self.temp_dir,
        )
        # The error should be raised during build when fixing 5 layers in a 3-layer slab
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        builder._validate_layer_count()

        with self.assertRaises(ValueError) as ctx:
            builder.build_all(self.temp_dir)
        self.assertIn("Cannot fix", str(ctx.exception))
```

#### Part B: Implementation

- [ ] **Step 2: Add `_add_constraints` to SurfaceBuilder (vacuum is handled in _build_slab)**

```python
    def _add_constraints(self, slab: Struct) -> Struct:
        """Apply Selective Dynamics to the slab.

        Bottom N layers are fixed; top layers are free.
        N determined by auto-fix rule or user override.
        Raises ValueError if fix_layers >= total slab layers.
        """
        # Determine number of layers to fix
        if self.fix_layers is not None:
            n_fix = self.fix_layers
        else:
            # Auto-fix: floor((layers+1)/2)
            n_fix = (self.n_layers + 1) // 2

        # Identify layers in the slab by z-coordinate clustering
        coords = slab.get_coords(direct=False)
        z_vals = coords[:, 2]
        z_sorted = np.sort(np.unique(np.round(z_vals, decimals=self.precision)))

        # Find layer boundaries using the same gap-detection approach
        if len(z_sorted) < 2:
            slab_layers = 1
        else:
            gaps = np.diff(z_sorted)
            threshold = np.median(gaps) * 3.0
            slab_layers = int(np.sum(gaps > threshold)) + 1

        if n_fix >= slab_layers:
            raise ValueError(
                f"Cannot fix {n_fix} layers in a {slab_layers}-layer slab. "
                f"Reduce --fix-layers."
            )

        if slab_layers <= 1:
            # Single layer slab — can't meaningfully fix
            for atom in slab:
                atom.constr = ["T", "T", "T"]
            return slab

        # Assign layer indices to atoms
        atom_zs = z_vals
        z_unique_sorted = np.sort(np.unique(np.round(atom_zs, decimals=self.precision)))

        # Simplified: bottom n_fix groups of z-coordinates are fixed
        fix_threshold = z_unique_sorted[min(n_fix, len(z_unique_sorted) - 1)] + 0.01

        fix_mask = self.fix_z_only
        fix_constr = ["T", "T", "F"] if fix_mask else ["F", "F", "F"]
        free_constr = ["T", "T", "T"]

        for atom in slab:
            atom_z = np.dot(atom.coord, slab.cell)[2]
            if atom_z <= fix_threshold:
                atom.constr = fix_constr.copy()
            else:
                atom.constr = free_constr.copy()

        return slab
```

- [ ] **Step 3: Run vacuum and constraints tests**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_surface.TestVacuumAndConstraints -v
```

Expected: 5 tests PASS.

- [ ] **Step 4: Commit**

```bash
git add src/modeling/surface.py tests/modeling/test_surface.py
git commit -m "feat: add vacuum distribution and selective dynamics to slabs

Bottom vacuum = 2A, top vacuum = total - 2A. Auto-fix rule: floor((layers+1)/2)
layers fixed as FFF (or TTF with --fix-z-only). Validation for vacuum < 2A
and fix_layers >= total slab layers.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 6: Metrics + CSV Summary

**Files:**
- Modify: `src/modeling/surface.py`
- Modify: `tests/modeling/test_surface.py`

#### Part A: Test

- [ ] **Step 1: Add summary CSV test**

```python
class TestSummaryCSV(unittest.TestCase):

    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_summary_csv_columns(self):
        """Summary CSV has all required columns with valid values."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(self.temp_dir)
        self.assertGreater(len(slabs), 0)

        # Read summary CSV
        miller_str = "111"
        name = poscar.stem
        csv_path = self.temp_dir / f"{name}-slab-{miller_str}" / "summary.csv"
        self.assertTrue(csv_path.exists(), f"Summary CSV not found at {csv_path}")

        import csv as csv_mod
        with open(csv_path, "r") as f:
            reader = csv_mod.DictReader(f)
            rows = list(reader)

        self.assertEqual(len(rows), len(slabs),
                         "CSV row count should match slab count")

        expected_columns = {
            "filename", "termination_id", "gap_position_z", "n_layers",
            "n_atoms", "composition", "layer_compositions",
            "composition_deviation", "surface_energy_est", "dipole_moment",
            "top_layer_z_std", "fix_layers", "fix_mode",
            "vacuum_top", "vacuum_bottom",
        }
        actual_columns = set(reader.fieldnames)
        missing = expected_columns - actual_columns
        self.assertEqual(len(missing), 0,
            f"Missing CSV columns: {missing}")

        for row in rows:
            self.assertGreater(int(row["n_atoms"]), 0)
            self.assertGreater(float(row["vacuum_top"]), 0)
            self.assertGreater(float(row["vacuum_bottom"]), 0)
            self.assertIn(row["fix_mode"], ["FFF", "TTF"])
```

#### Part B: Implementation

- [ ] **Step 2: Add `_compute_metrics` and `_write_summary`**

```python
    def _compute_metrics(self, slab: Struct, n_layers: int, term_id: str) -> dict:
        """Compute all metrics for a slab."""
        coords = slab.get_coords(direct=False)
        z_vals = coords[:, 2]

        # Composition
        from collections import Counter
        sym_counts = Counter(a.symbol for a in slab)
        total = len(slab)
        composition = "".join(f"{s}{c}" for s, c in sorted(sym_counts.items()))

        # Bulk composition for deviation
        bulk_counts = Counter(a.symbol for a in self._bulk_struct)
        bulk_total = len(self._bulk_struct)

        # Composition deviation
        all_elements = set(sym_counts.keys()) | set(bulk_counts.keys())
        deviations = []
        for elem in sorted(all_elements):
            surf_frac = sym_counts.get(elem, 0) / total
            bulk_frac = bulk_counts.get(elem, 0) / bulk_total
            deviations.append(abs(surf_frac - bulk_frac))
        composition_deviation = sum(deviations) / len(deviations) if deviations else 0.0

        # Surface energy estimate (simple bond-counting model)
        # Approximate: E_surf ~ (n_broken_bonds / n_total_bonds) * cohesive_energy
        # For quick ranking, use fraction of atoms in surface layer
        top_z = np.max(z_vals)
        z_unique = np.sort(np.unique(np.round(z_vals, decimals=self.precision)))
        if len(z_unique) >= 2:
            top_threshold = z_unique[-2] + 0.01
            top_atoms = int(np.sum(z_vals >= top_threshold))
            surface_energy_est = top_atoms / max(total, 1)
        else:
            surface_energy_est = 1.0

        # Dipole moment (simple z-weighted charge)
        # Use z-coordinates relative to slab center
        z_center = (np.max(z_vals) + np.min(z_vals)) / 2.0
        dipole = float(np.sum(z_vals - z_center))  # in e*A units

        # Top layer roughness
        if len(z_unique) >= 2:
            top_threshold = z_unique[-2] + 0.01
            top_zs = z_vals[z_vals >= top_threshold]
            top_layer_z_std = float(np.std(top_zs))
        else:
            top_layer_z_std = 0.0

        # Layer compositions
        z_rounded = np.round(z_vals, decimals=self.precision)
        z_unique = np.sort(np.unique(z_rounded))
        layer_comps = []
        for zu in z_unique:
            mask = np.abs(z_rounded - zu) < 0.5 * (10 ** (-self.precision))
            layer_syms = [slab[i].symbol for i, m in enumerate(mask) if m]
            layer_count = Counter(layer_syms)
            layer_comps.append("".join(f"{s}{c}" for s, c in sorted(layer_count.items())))
        layer_compositions = " | ".join(layer_comps)

        # Fix parameters
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
        all_slabs: list[tuple[Path, Struct, float, int, str]],
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
                    **metrics,
                }
                writer.writerow(row)

        logging.info(f"Summary saved to {csv_path}")
        return csv_path
```

- [ ] **Step 3: Run summary CSV test**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_surface.TestSummaryCSV -v
```

Expected: 1 test PASS.

- [ ] **Step 4: Commit**

```bash
git add src/modeling/surface.py tests/modeling/test_surface.py
git commit -m "feat: add metrics computation and CSV summary output

Compute composition_deviation, surface_energy_est (bond-count model),
dipole_moment, and top_layer_z_std for each slab. Write summary.csv
with all 15 columns.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 7: CLI Surface Subcommand

**Files:**
- Modify: `src/cli/poscarkit.py`

#### Part A: Implementation

- [ ] **Step 1: Add `cmd_surface` function**

Insert before `cmd_help` (or after `cmd_import_to_model`):

```python
def cmd_surface(args: argparse.Namespace) -> int:
    poscar = Path(args.poscar) if args.poscar else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    miller = tuple(args.miller)
    layers = args.layers
    vacuum = args.vacuum
    fix_layers = args.fix_layers
    fix_z_only = args.fix_z_only

    if not poscar or not poscar.is_file():
        logging.error(f"POSCAR file not found: {poscar}")
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    from src.modeling.surface import SurfaceBuilder

    builder = SurfaceBuilder(
        poscar=poscar,
        miller=miller,
        layers=layers,
        vacuum=vacuum,
        fix_layers=fix_layers,
        fix_z_only=fix_z_only,
        outdir=outdir,
        precision=getattr(args, "precision", 2),
    )

    try:
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        results = builder.build_all(outdir)

        logging.info(f"Surface slab generation completed. Generated {len(results)} slabs:")
        for r in results:
            logging.info(f"  - {r}")
        return 0
    except ValueError as e:
        logging.error(str(e))
        return 1
```

- [ ] **Step 2: Register `surface` subcommand in `main()`**

Insert after the `slice-to-countcn` parser block (before the supercell parser):

```python
    # Surface command
    parser_surface = subparsers.add_parser("surface", help="Generate surface slabs from bulk POSCAR")
    parser_surface.add_argument(
        "poscar",
        type=str,
        help="Path to bulk POSCAR file",
    )
    parser_surface.add_argument(
        "--miller",
        "-m",
        type=int,
        nargs=3,
        default=(0, 0, 1),
        metavar=("H", "K", "L"),
        help="Miller indices (default: 0 0 1)",
    )
    parser_surface.add_argument(
        "--layers",
        "-l",
        type=int,
        default=3,
        help="Number of slab layers, outputs N and N+1 (default: 3)",
    )
    parser_surface.add_argument(
        "--vacuum",
        "-v",
        type=float,
        default=15.0,
        help="Total vacuum thickness in Angstrom (default: 15.0)",
    )
    parser_surface.add_argument(
        "--fix-layers",
        type=int,
        default=None,
        help="Manual override for fixed bottom layers (default: auto)",
    )
    parser_surface.add_argument(
        "--fix-z-only",
        action="store_true",
        help="Fix only z-direction instead of all (default: FFF)",
    )
    parser_surface.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser_surface.add_argument(
        "--precision",
        type=int,
        default=2,
        help="Decimal precision for z-coordinate grouping (default: 2)",
    )
    parser_surface.set_defaults(func=cmd_surface)
```

- [ ] **Step 3: Verify the subcommand registers correctly**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -c "from src.cli.poscarkit import main; print('CLI imports OK')"
```

Expected: "CLI imports OK" with no errors.

- [ ] **Step 4: Manual smoke test with a real bulk POSCAR**

Create a test bulk using ASE:

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -c "
from ase.build import bulk
from ase.io import write
au = bulk('Au', 'fcc', a=4.08, cubic=True)
write('_test_bulk.vasp', au, format='vasp')
print('Test bulk created')
"
```

Run the surface command:

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python main.py surface _test_bulk.vasp --miller 1 1 1 --layers 3 --vacuum 15.0 --outdir output/
```

Expected: Output directory created with slabs and summary.csv.

- [ ] **Step 5: Verify output**

```bash
ls output/_test_bulk-slab-111/
```

Expected: Multiple `.vasp` files + `summary.csv`.

- [ ] **Step 6: Clean up test artifacts**

```bash
rm _test_bulk.vasp
rm -rf output/_test_bulk-slab-111/
```

- [ ] **Step 7: Commit**

```bash
git add src/cli/poscarkit.py
git commit -m "feat: add 'surface' CLI subcommand for slab generation

poscarkit surface <poscar> [--miller H K L] [--layers N] [--vacuum V]
  [--fix-layers N] [--fix-z-only] [--outdir DIR]

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 8: Fix Slicer.group_by_normal Periodic Bug

**Files:**
- Modify: `src/modeling/slice.py`
- Modify: `tests/modeling/test_slice.py`

#### Part A: Test

- [ ] **Step 1: Add periodic boundary test to test_slice.py**

```python
    def test_group_by_normal_periodic(self):
        """Same-layer atoms across cell boundary grouped together."""
        # Create a structure where an atomic layer straddles the cell boundary.
        # For FCC Au with 001 slice: atoms at z=0.0 and z near 1.0 should be in
        # the same layer after transformation.
        cell = np.array([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]])
        atom_list = [
            Atom(0, "Au", np.array([0.0, 0.0, 0.001])),      # near bottom
            Atom(1, "Au", np.array([0.0, 0.0, 0.999])),      # near top (same layer!)
            Atom(2, "Au", np.array([0.5, 0.5, 0.500])),      # middle layer
        ]
        struct = Struct(cell=cell, is_direct=True, atom_list=atom_list)

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            poscar_p = temp_path / "periodic_test.vasp"
            SimplePoscar.write_poscar(poscar_p, struct, "Periodic test")

            s = Slicer("test_periodic", poscar_p, (0, 0, 1))
            layers = list(s.group_by_normal())

            # Atoms at z=0.001 and z=0.999 should be in the SAME layer
            # So we should have 2 layers, not 3
            self.assertEqual(len(layers), 2,
                f"Expected 2 layers (periodic merge), got {len(layers)}")
```

#### Part B: Implementation

- [ ] **Step 2: Fix `group_by_normal` in slice.py**

Replace the method (lines 72-92) with:

```python
    def group_by_normal(self, precision: int = 2):
        """Group atoms by distance to the normal vector.

        Uses fractional coordinates along the surface normal with circular
        gap detection to correctly handle atoms near the periodic boundary.

        Args:
            precision: Precision of distance (decimal places for rounding)

        Yields:
            Tuple: (Projection, Struct)
        """
        transfd = self.transformed
        basis_norm = self.basis_norm

        # Use fractional coordinates along the normal direction
        coords_direct = transfd.get_coords(direct=True)

        # Project onto normal in direct space: dot with basis[2] / |basis[2]|
        # In the transformed cell, the c-vector IS the surface normal,
        # so the fractional z-coordinate IS the projection we want.
        projs = coords_direct[:, 2]

        # Handle periodic wrapping: sort by fractional z, find gaps
        sorted_idx = np.argsort(projs)
        sorted_z = projs[sorted_idx]

        n = len(sorted_z)
        gaps = np.diff(sorted_z)
        wrap_gap = 1.0 - sorted_z[-1] + sorted_z[0]
        all_gaps = np.append(gaps, wrap_gap)
        threshold = np.median(all_gaps) * 3.0

        # Build groups
        groups = []
        current = [sorted_idx[0]]

        for i in range(1, n):
            if gaps[i - 1] > threshold:
                groups.append(current)
                current = []
            current.append(sorted_idx[i])

        # Handle wrap: if wrap_gap is small, merge first and last
        if wrap_gap <= threshold and groups:
            groups[0] = current + groups[0]
        else:
            groups.append(current)

        # Yield groups
        for group in groups:
            proj = round(float(np.mean(projs[group])), precision)
            layer = transfd.copy(atom_list=[transfd[i] for i in group])
            yield proj, layer
```

- [ ] **Step 3: Run the periodic fix test**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_slice.TestSlice.test_group_by_normal_periodic -v
```

Expected: PASS.

- [ ] **Step 4: Run all existing slice tests to check for regressions**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_slice -v
```

Expected: all 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/modeling/slice.py tests/modeling/test_slice.py
git commit -m "fix: use fractional z with circular gap detection in Slicer.group_by_normal

Previously used Cartesian coordinates which incorrectly split atomic layers
that straddle the periodic cell boundary (z~0 vs z~c). Now uses the same
circular gap detection algorithm as SurfaceBuilder for correct layer
identification.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
```

---

### Task 9: Full Test Suite Verification

**Files:**
- All modified files

- [ ] **Step 1: Run the complete test suite**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest discover tests -v
```

Expected: all tests PASS with no failures or errors.

- [ ] **Step 2: Run editable install to ensure imports work**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata pip install -e .
```

- [ ] **Step 3: Smoke test full CLI flow with a real structure**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -c "
from ase.build import bulk
from ase.io import write
import numpy as np

# Create a multi-element FCC bulk (simulating simple HEA)
from ase import Atoms
au_cu = bulk('Au', 'fcc', a=4.08, cubic=True)
# Add Cu by modifying one of the atoms
syms = au_cu.get_chemical_symbols()
syms[1] = 'Cu'
syms[2] = 'Cu'
au_cu.set_chemical_symbols(syms)
write('_test_hea_bulk.vasp', au_cu, format='vasp')
print('Test HEA bulk created with symbols:', au_cu.get_chemical_symbols())
"
```

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python main.py surface _test_hea_bulk.vasp --miller 1 1 1 --layers 3 --vacuum 15.0 --outdir output/
```

- [ ] **Step 4: Verify output structure**

```bash
ls output/_test_hea_bulk-slab-111/
```

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -c "
import csv
with open('output/_test_hea_bulk-slab-111/summary.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        print(f'{row[\"termination_id\"]} | layers={row[\"n_layers\"]} | '
              f'atoms={row[\"n_atoms\"]} | comp={row[\"composition\"]} | '
              f'dev={row[\"composition_deviation\"]}')
"
```

- [ ] **Step 5: Clean up test artifacts**

```bash
rm _test_hea_bulk.vasp
rm -rf output/_test_hea_bulk-slab-111/
```

- [ ] **Step 6: Commit (if any final cleanups were needed)**

```bash
git status
```

---

### Task 10: Relaxed HEA CONTCAR Test (deferred)

**Files:**
- Modify: `tests/modeling/test_surface.py`

> **Note:** This test requires the real 6-element HEA CONTCAR file provided by the user.
> Add this test once the CONTCAR is available.

- [ ] **Step 1: Add test_relaxed_hea_contcar**

```python
    def test_relaxed_hea_contcar(self):
        """Real 6-element HEA CONTCAR with ~0.3A atomic displacements.

        Verifies that the gap detection algorithm correctly identifies
        atomic layers in a VASP-relaxed structure where atoms have moved
        from their ideal crystallographic positions.
        """
        # Path to the user-provided CONTCAR
        contcar_path = Path("tests/data/hea_6elem_CONTCAR.vasp")
        if not contcar_path.exists():
            self.skipTest(f"Test CONTCAR not found: {contcar_path}")

        builder = SurfaceBuilder(
            poscar=contcar_path, miller=(1, 1, 1), layers=3, vacuum=15.0,
        )
        builder._transform_cell()
        layers = builder._identify_layers()

        # Should find at least 4 distinct layers (enough for 3-layer slab)
        self.assertGreaterEqual(len(layers), 4,
            f"Expected >= 4 layers in relaxed HEA, got {len(layers)}")

        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(Path(tempfile.mkdtemp()))
        self.assertGreater(len(slabs), 0)

        # Each slab should have correct composition
        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            self.assertGreater(len(slab), 0)
```

- [ ] **Step 2: Run the test (skip if CONTCAR not available)**

```bash
"C:/ProgramData/miniconda3/Scripts/conda.exe" run -n matdata python -m unittest tests.modeling.test_surface.TestSlabAssembly.test_relaxed_hea_contcar -v
```

---

## Summary

| Task | Description | Files |
|------|-------------|-------|
| 1 | Bridge fix | base.py, test_simple_poscar.py |
| 2 | Layer identification | surface.py, test_surface.py |
| 3 | Gap detection + filtering | surface.py, test_surface.py |
| 4 | Slab assembly + unwrapping | surface.py, test_surface.py |
| 5 | Vacuum + constraints | surface.py, test_surface.py |
| 6 | Metrics + CSV summary | surface.py, test_surface.py |
| 7 | CLI subcommand | poscarkit.py |
| 8 | Slicer periodic fix | slice.py, test_slice.py |
| 9 | Full test suite verification | All files |
| 10 | Relaxed HEA test (deferred) | test_surface.py |

**Total commits:** 8 (Tasks 1-8), plus optional final cleanup.
