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
            poscar=poscar, miller=(0, 0, 1), layers=1, vacuum=15.0,
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
            # c-vector should be oriented purely in z for (001)
            self.assertAlmostEqual(slab.cell[2, 0], 0.0)
            self.assertAlmostEqual(slab.cell[2, 1], 0.0)
            self.assertGreater(slab.cell[2, 2], 0.0)

    def test_build_slab_n_plus_one(self):
        """N+1 (2-layer) slabs are produced alongside N (1-layer) slabs."""
        poscar = _make_fcc_bulk_poscar(self.temp_dir)
        builder = SurfaceBuilder(
            poscar=poscar, miller=(0, 0, 1), layers=1, vacuum=15.0,
            outdir=self.temp_dir,
        )
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()

        slabs = builder.build_all(self.temp_dir)
        layer_counts = set()
        for slab_path in slabs:
            name = slab_path.stem
            if "layers1" in name:
                layer_counts.add(1)
            elif "layers2" in name:
                layer_counts.add(2)

        self.assertIn(1, layer_counts, "Should produce 1-layer slabs")
        self.assertIn(2, layer_counts, "Should produce 2-layer slabs (N+1)")

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


if __name__ == "__main__":
    unittest.main()
