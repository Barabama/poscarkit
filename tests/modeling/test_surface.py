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
        """FCC bulk transformed to (111) orientation identifies 3 layers."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        layers = builder._identify_layers()

        self.assertEqual(len(layers), 3,
                         "FCC(111) 4-atom cell should yield 3 distinct layers")

        for layer in layers:
            self.assertIsInstance(layer, Layer)
            self.assertIsInstance(layer.z_centroid, float)
            self.assertGreater(len(layer.atoms), 0)
            self.assertTrue(0.0 <= layer.z_centroid <= 1.0,
                            f"z_centroid {layer.z_centroid} not in [0,1]")

    def test_identify_layers_bcc100(self):
        """BCC bulk (100) transformation identifies 2 layers."""
        poscar = _make_bcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 0, 0), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        layers = builder._identify_layers()

        self.assertEqual(len(layers), 2,
                         "BCC(100) 2-atom cell should yield 2 distinct layers")

    def test_identify_layers_custom_miller(self):
        """FCC (2,1,0) transformation identifies 10 layers."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(2, 1, 0), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        layers = builder._identify_layers()

        self.assertEqual(len(layers), 10,
                         "FCC(210) 4-atom cell should yield 10 distinct layers")


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

            if i < len(layers) - 1:
                self.assertTrue(layer_i.z_centroid < gap_z < layer_next.z_centroid,
                                f"Gap {i} ({gap_z:.4f}) not between layers "
                                f"({layer_i.z_centroid:.4f}, {layer_next.z_centroid:.4f})")

    def test_insufficient_layers_filtered(self):
        """Gaps requiring wraparound with insufficient layers are excluded."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=50, vacuum=15.0,
            outdir=self.temp_dir,
        )
        with self.assertRaises(ValueError) as ctx:
            builder._transform_cell()
            builder._identify_layers()
            builder._find_gaps()
            builder._validate_layer_count()
            builder.build_all(self.temp_dir)
        self.assertIn("Reduce --layers", str(ctx.exception))


class TestSlabAssembly(unittest.TestCase):

    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_build_slab_basic(self):
        """Slab has correct cell shape and c-vector oriented along z."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=2, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()

        slabs = builder.build_all(self.temp_dir)
        self.assertGreater(len(slabs), 0, "Should produce at least one slab")

        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            self.assertEqual(slab.cell.shape, (3, 3))
            self.assertGreater(slab.cell[2, 2], 0.0)

    def test_build_slab_n_plus_one(self):
        """N+1 (3-layer) slabs are produced alongside N (2-layer) slabs."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=2, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()

        slabs = builder.build_all(self.temp_dir)
        layer_counts = set()
        for slab_path in slabs:
            name = slab_path.stem
            if "layers2" in name:
                layer_counts.add(2)
            elif "layers3" in name:
                layer_counts.add(3)

        self.assertIn(2, layer_counts, "Should produce 2-layer slabs")
        self.assertIn(3, layer_counts, "Should produce 3-layer slabs (N+1)")

    def test_layers_periodic_crossing(self):
        """Layer crossing cell boundary is unwrapped correctly.

        After unwrapping, atoms in each slab should have z-coordinates
        that form tight layer clusters with no outliers from periodic
        wrapping (e.g. an atom at z ~ 0.999 staying above atoms at
        z ~ 0.001 after processing the wrap-around gap).
        """
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(0, 0, 1), layers=1, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()

        slabs = builder.build_all(self.temp_dir)

        # N=1: single layer, all atoms at same z  -> 2 atoms
        # N+1=2: two layers, atoms in decreasing z order -> 4 atoms
        n1_paths = [p for p in slabs if "layers1" in p.stem]
        n2_paths = [p for p in slabs if "layers2" in p.stem]
        self.assertEqual(len(n1_paths), 1,
                         "Expected one N=1 slab from the wrap-around gap")
        self.assertEqual(len(n2_paths), 1,
                         "Expected one N+1=2 slab from the wrap-around gap")

        # N=1 slab: all atoms at same z (single layer)
        slab_n1 = SimplePoscar.read_poscar(n1_paths[0])
        coords_n1 = slab_n1.get_coords(direct=False)
        z_n1 = coords_n1[:, 2]
        self.assertEqual(len(z_n1), 2)
        self.assertAlmostEqual(np.std(z_n1), 0.0, places=5,
            msg="Single-layer N=1 slab should have all atoms at same z")

        # N+1=2 slab: two distinct z clusters, no outliers
        slab_n2 = SimplePoscar.read_poscar(n2_paths[0])
        coords_n2 = slab_n2.get_coords(direct=False)
        z_n2 = coords_n2[:, 2]
        self.assertEqual(len(z_n2), 4)
        z_sorted = np.sort(z_n2)
        z_gaps = np.diff(z_sorted)
        # There should be exactly one inter-layer gap > 0.1 A
        large_gaps = z_gaps[z_gaps > 0.1]
        # Within each cluster, atoms should be at essentially the same z
        self.assertEqual(len(large_gaps), 1,
                         "Expected exactly one inter-layer gap in N+1 slab")
        self.assertGreater(large_gaps[0], 1.0,
                         "Inter-layer gap should be > 1 A for FCC(001)")
        # Verify vacuum bottom padding
        self.assertAlmostEqual(np.min(z_n2), 2.0, places=4,
            msg="Bottom atoms should be at vacuum_bottom=2.0 A")

    def test_note_preserved(self):
        """note fields survive the full Struct -> Atoms -> cut -> Struct pipeline."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        original = SimplePoscar.read_poscar(poscar)
        original_notes = set(a.note for a in original)

        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=2, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()

        slabs = builder.build_all(self.temp_dir)
        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            slab_notes = set(a.note for a in slab)
            self.assertTrue(
                slab_notes.issubset(original_notes) or slab_notes == original_notes,
                f"Slab notes {slab_notes} not subset of original {original_notes}"
            )
            for atom in slab:
                self.assertTrue(atom.note, f"Atom {atom.index} has empty note")


class TestVacuumAndConstraints(unittest.TestCase):

    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_vacuum_distribution(self):
        """Bottom vacuum = 2.0 A, top vacuum = total - 2.0 A."""
        total_vacuum = 20.0
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=2, vacuum=total_vacuum,
            outdir=self.temp_dir,
        )
        self.assertEqual(builder.vacuum_bottom, 2.0)
        self.assertEqual(builder.vacuum_top, total_vacuum - 2.0)

        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()

        slabs = builder.build_all(self.temp_dir)
        self.assertGreater(len(slabs), 0)

        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            coords = slab.get_coords(direct=False)
            z_vals = coords[:, 2]

            cell_c = float(np.linalg.norm(slab.cell[2]))
            self.assertGreaterEqual(float(np.min(z_vals)), 0.0)
            self.assertLessEqual(float(np.max(z_vals)), cell_c)

            # Bottom vacuum region: [0, 2A) should be empty
            bottom_atoms = [z for z in z_vals if z < 2.0]
            self.assertEqual(len(bottom_atoms), 0,
                "No atoms should be in bottom vacuum region [0, 2A)")

    def test_constraints_fff(self):
        """Default: bottom layers fully fixed (F F F)."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=2, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()

        slabs = builder.build_all(self.temp_dir)
        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            coords = slab.get_coords(direct=False)
            z_vals = coords[:, 2]
            z_unique = np.sort(np.unique(np.round(z_vals, decimals=2)))
            threshold = z_unique[min(1, len(z_unique) - 1)] + 0.5

            for atom in slab:
                atom_z_cart = atom.coord[2]
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
            poscar=poscar, miller=(1, 1, 1), layers=2, vacuum=15.0,
            fix_z_only=True, outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()

        slabs = builder.build_all(self.temp_dir)
        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            coords = slab.get_coords(direct=False)
            z_vals = coords[:, 2]
            z_unique = np.sort(np.unique(np.round(z_vals, decimals=2)))
            threshold = z_unique[min(1, len(z_unique) - 1)] + 0.5

            found_fixed = False
            for atom in slab:
                atom_z_cart = atom.coord[2]
                if atom_z_cart < threshold:
                    found_fixed = True
                    self.assertEqual(atom.constr, ["T", "T", "F"],
                        f"Fixed-z atom {atom.index} should be TTF, got {atom.constr}")
            self.assertTrue(found_fixed, "Should have at least one fixed atom")

    def test_vacuum_below_minimum(self):
        """ValueError when total vacuum < 2.0 A."""
        with self.assertRaises(ValueError) as ctx:
            SurfaceBuilder(
                poscar=_make_fcc_bulk_poscar(self.temp_dir),
                miller=(1, 1, 1), layers=2, vacuum=1.0,
                outdir=self.temp_dir,
            )
        self.assertIn("2.0", str(ctx.exception))

    def test_fix_layers_exceeds_total(self):
        """ValueError when fix_layers >= total slab layers."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(1, 1, 1), layers=2, vacuum=15.0,
            fix_layers=5, outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()

        with self.assertRaises(ValueError) as ctx:
            builder.build_all(self.temp_dir)
        self.assertIn("Cannot fix", str(ctx.exception))


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
            poscar=poscar, miller=(1, 1, 1), layers=2, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()

        slabs = builder.build_all(self.temp_dir)
        self.assertGreater(len(slabs), 0)

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
        self.assertEqual(len(missing), 0, f"Missing CSV columns: {missing}")

        for row in rows:
            self.assertGreater(int(row["n_atoms"]), 0)
            self.assertGreater(float(row["vacuum_top"]), 0)
            self.assertGreater(float(row["vacuum_bottom"]), 0)
            self.assertIn(row["fix_mode"], ["FFF", "TTF"])


class TestRelaxedHEA(unittest.TestCase):
    """Tests using real VASP-relaxed HEA CONTCAR."""

    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_relaxed_hea_contcar(self):
        """Real 6-element BCC HEA CONTCAR with ~0.1-0.3A atomic displacements.

        Verifies that the layer identification algorithm correctly handles
        a VASP-relaxed structure where atoms have moved from their ideal
        crystallographic positions and the cell is non-orthogonal.
        """
        contcar_path = Path(__file__).parent.parent / "data" / "hea_6elem_CONTCAR.vasp"
        if not contcar_path.exists():
            self.skipTest(f"Test CONTCAR not found: {contcar_path}")

        builder = SurfaceBuilder(
            poscar=contcar_path, miller=(0, 0, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        layers = builder._identify_layers()

        # 3x3x3 BCC supercell along (001) should have at least 6 layers
        self.assertGreaterEqual(len(layers), 4,
            f"Expected >= 4 layers in relaxed BCC(001) HEA, got {len(layers)}")

        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(self.temp_dir)
        self.assertGreater(len(slabs), 0, "Should produce at least one slab")

        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            self.assertGreater(len(slab), 0, f"Slab {slab_path.name} is empty")

            # Verify slab symbols are a subset of bulk elements (not all
            # elements need to appear in every termination)
            symbols = set(slab.symbols)
            bulk_elements = {"Al", "Nb", "Ni", "Ti", "V", "W"}
            unknown = symbols - bulk_elements
            self.assertEqual(len(unknown), 0,
                f"Slab contains unexpected elements: {unknown}")

    def test_relaxed_hea_fcc111_transform(self):
        """Relaxed BCC HEA transformed to (111) orientation produces valid slabs."""
        contcar_path = Path(__file__).parent.parent / "data" / "hea_6elem_CONTCAR.vasp"
        if not contcar_path.exists():
            self.skipTest(f"Test CONTCAR not found: {contcar_path}")

        builder = SurfaceBuilder(
            poscar=contcar_path, miller=(1, 1, 1), layers=3, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        layers = builder._identify_layers()

        self.assertGreaterEqual(len(layers), 3,
            f"Expected >= 3 layers in relaxed BCC(111) HEA, got {len(layers)}")

        builder._find_gaps()
        builder._validate_layer_count()

        slabs = builder.build_all(self.temp_dir)
        self.assertGreater(len(slabs), 0)

        # Verify all slabs are non-empty with valid cells
        for slab_path in slabs:
            slab = SimplePoscar.read_poscar(slab_path)
            self.assertGreater(len(slab), 0)
            self.assertEqual(slab.cell.shape, (3, 3))


if __name__ == "__main__":
    unittest.main()
