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


if __name__ == "__main__":
    unittest.main()
