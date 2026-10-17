import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.modeling.slice import Slicer
from src.modeling.base import Struct, Atom, SimplePoscar
from src.modeling.layers import DEFAULT_LAYER_GAP, layer_indices, scaled_basis


class TestSlice(unittest.TestCase):
    """Test cases for slice module."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a simple cubic structure for testing
        cell = np.array([[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]])

        atom_list = [
            Atom(index=0, symbol="Al", coord=np.array([0.0, 0.0, 0.0])),
            Atom(index=1, symbol="Al", coord=np.array([0.5, 0.5, 0.0])),
            Atom(index=2, symbol="Al", coord=np.array([0.5, 0.0, 0.5])),
            Atom(index=3, symbol="Al", coord=np.array([0.0, 0.5, 0.5])),
        ]

        self.test_struct = Struct(cell=cell, is_direct=True, atom_list=atom_list)

        # Create temporary directory and POSCAR file
        self.temp_dir = Path(tempfile.mkdtemp())
        self.poscar_path = self.temp_dir / "POSCAR"
        SimplePoscar.write_poscar(self.poscar_path, self.test_struct, "Test structure")

    def tearDown(self):
        """Tear down test fixtures."""
        # Clean up temporary directory
        import shutil

        shutil.rmtree(self.temp_dir)

    def test_slice_001_direction(self):
        """Test slicing along (001) direction."""
        # Create slicer for (001) direction
        s = Slicer("test_001", self.poscar_path, (0, 0, 1))

        # Test that basis vectors are correctly calculated
        expected_basis = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        np.testing.assert_array_equal(s.basis, expected_basis)

        # Test slice2files method
        with tempfile.TemporaryDirectory() as output_dir:
            outdir = Path(output_dir)
            layer_files = s.slice2files(outdir)

            # Check that we got a list of layer files
            self.assertIsInstance(layer_files, list)
            self.assertTrue(all(isinstance(f, Path) for f in layer_files))
            self.assertGreater(len(layer_files), 0)

            # Check that output directory was created based on first file
            result_dir = layer_files[0].parent
            self.assertTrue(result_dir.exists())

            # Check that transformed structure file was created
            transformed_file = result_dir / "Transformed-(001).vasp"
            self.assertTrue(transformed_file.exists())

            # Check that layer files were created
            layer_files_in_dir = list(result_dir.glob("Transformed-(001)-layer*.vasp"))
            self.assertGreater(len(layer_files_in_dir), 0)

            # Check that corresponding image files were created
            image_files = list(result_dir.glob("Transformed-(001)-layer*.png"))
            self.assertEqual(len(image_files), len(layer_files_in_dir))

    def test_slice_110_direction(self):
        """Test slicing along (110) direction."""
        # Create slicer for (110) direction
        s = Slicer("test_110", self.poscar_path, (1, 1, 0))

        # Test that basis vectors are correctly calculated
        expected_basis = np.array([[0, 0, 1], [0.5, -0.5, 0], [1, 1, 0]])
        np.testing.assert_array_equal(s.basis, expected_basis)

        # Test slice2files method
        with tempfile.TemporaryDirectory() as output_dir:
            outdir = Path(output_dir)
            layer_files = s.slice2files(outdir)

            # Check that we got a list of layer files
            self.assertIsInstance(layer_files, list)
            self.assertTrue(all(isinstance(f, Path) for f in layer_files))
            self.assertGreater(len(layer_files), 0)

            # Check that output directory was created based on first file
            result_dir = layer_files[0].parent
            self.assertTrue(result_dir.exists())

            # Check that transformed structure file was created
            transformed_file = result_dir / "Transformed-(110).vasp"
            self.assertTrue(transformed_file.exists())

            # Check that layer files were created
            layer_files_in_dir = list(result_dir.glob("Transformed-(110)-layer*.vasp"))
            self.assertGreater(len(layer_files_in_dir), 0)

            # Check that corresponding image files were created
            image_files = list(result_dir.glob("Transformed-(110)-layer*.png"))
            self.assertEqual(len(image_files), len(layer_files_in_dir))

    def test_slice_111_direction(self):
        """Test slicing along (111) direction."""
        # Create slicer for (111) direction
        s = Slicer("test_111", self.poscar_path, (1, 1, 1))

        # Test that basis vectors are correctly calculated
        expected_basis = np.array([[0.5, -0.5, 0], [0.5, 0.5, -1], [1, 1, 1]])
        np.testing.assert_array_equal(s.basis, expected_basis)

        # Test slice2files method
        with tempfile.TemporaryDirectory() as output_dir:
            outdir = Path(output_dir)
            layer_files = s.slice2files(outdir)

            # Check that we got a list of layer files
            self.assertIsInstance(layer_files, list)
            self.assertTrue(all(isinstance(f, Path) for f in layer_files))
            self.assertGreater(len(layer_files), 0)

            # Check that output directory was created based on first file
            result_dir = layer_files[0].parent
            self.assertTrue(result_dir.exists())

            # Check that transformed structure file was created
            transformed_file = result_dir / "Transformed-(111).vasp"
            self.assertTrue(transformed_file.exists())

            # Check that layer files were created
            layer_files_in_dir = list(result_dir.glob("Transformed-(111)-layer*.vasp"))
            self.assertGreater(len(layer_files_in_dir), 0)

            # Check that corresponding image files were created
            image_files = list(result_dir.glob("Transformed-(111)-layer*.png"))
            self.assertEqual(len(image_files), len(layer_files_in_dir))

    def test_group_by_normal(self):
        """Test grouping atoms by normal vector."""
        s = Slicer("test_grouping", self.poscar_path, (0, 0, 1))

        # Group atoms by normal
        layers = list(s.group_by_normal())

        # Check that we get layers
        self.assertGreater(len(layers), 0)

        # Check that each layer has a projection value and a Struct
        for proj, layer in layers:
            self.assertIsInstance(proj, (int, float))
            self.assertIsInstance(layer, Struct)

    def test_group_by_normal_periodic(self):
        """Same-layer atoms across cell boundary grouped together."""
        cell = np.array([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]])
        atom_list = [
            Atom(0, "Au", np.array([0.0, 0.0, 0.001])),
            Atom(1, "Au", np.array([0.0, 0.0, 0.999])),
            Atom(2, "Au", np.array([0.5, 0.5, 0.500])),
        ]
        struct = Struct(cell=cell, is_direct=True, atom_list=atom_list)

        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            poscar_p = temp_path / "periodic_test.vasp"
            SimplePoscar.write_poscar(poscar_p, struct, "Periodic test")

            s = Slicer("test_periodic", poscar_p, (0, 0, 1))

            # The cut must keep the atom sitting on the upper boundary.
            self.assertEqual(len(s.transformed), 3,
                "cut() dropped a boundary atom instead of keeping the periodic image")

            layers = list(s.group_by_normal())

            self.assertEqual(len(layers), 2,
                f"Expected 2 layers (periodic merge), got {len(layers)}")
            sizes = sorted(len(layer) for _, layer in layers)
            self.assertEqual(sizes, [1, 2],
                f"Expected the boundary layer to hold 2 atoms, got {sizes}")


class TestLayerClustering(unittest.TestCase):
    """Regression tests for gap-based layer identification.

    These target the failure mode where a thermally broadened layer is split
    by a fixed-width decimal bin: the layer count must not depend on where the
    crystal sits relative to the rounding grid.
    """

    def _distorted_structure(self, n_layers=40, jitter=0.01, shift=0.0):
        """Cubic cell with ``n_layers`` planes, each broadened by +/-jitter (A)."""
        length = float(n_layers)  # 1 A interlayer spacing
        cell = np.eye(3) * length
        atoms = []
        for k in range(n_layers):
            for i in range(2):
                for j in range(2):
                    sign = 1.0 if (i + j) % 2 else -1.0
                    z = (k / n_layers) + shift + sign * jitter / length
                    atoms.append(Atom(
                        index=len(atoms),
                        symbol="Cu",
                        coord=np.array([i / 2.0, j / 2.0, z]),
                    ))
        return Struct(cell=cell, is_direct=True, atom_list=atoms)

    def test_layer_count_is_translation_invariant(self):
        """Rigidly shifting a structure must not change its layer count."""
        for shift in (0.0, 0.001, 0.003, 0.005):
            struct = self._distorted_structure(n_layers=40, shift=shift)
            coords = struct.get_coords(direct=False)
            groups = layer_indices(coords, struct.cell, gap=DEFAULT_LAYER_GAP)
            self.assertEqual(
                len(groups), 40,
                f"shift={shift}: expected 40 layers, got {len(groups)}",
            )
            self.assertTrue(all(len(idx) == 4 for _, idx in groups))

    def test_gap_clustering_keeps_broadened_layer_together(self):
        """A 0.02 A wide layer must not be split by a 0.3 A gap threshold."""
        struct = self._distorted_structure(n_layers=40, jitter=0.01)
        coords = struct.get_coords(direct=False)
        groups = layer_indices(coords, struct.cell, gap=DEFAULT_LAYER_GAP)

        self.assertEqual(len(groups), 40)
        # Every atom belongs to exactly one layer.
        assigned = np.concatenate([idx for _, idx in groups])
        self.assertEqual(len(assigned), len(coords))
        self.assertEqual(len(np.unique(assigned)), len(coords))

    def test_basis_is_commensurate(self):
        """Scaled bases must span an integer number of cells."""
        expected = {
            (0, 0, 1): 1.0,
            (1, 1, 0): 1.0,
            (1, 1, 1): 1.5,
        }
        for miller, det in expected.items():
            basis = scaled_basis(miller)
            self.assertAlmostEqual(
                abs(float(np.linalg.det(basis))), det, places=10,
                msg=f"{miller}: basis is not commensurate",
            )
            # The third vector is the surface normal, so a and b must be
            # perpendicular to it.
            self.assertAlmostEqual(float(np.dot(basis[0], basis[2])), 0.0, places=10)
            self.assertAlmostEqual(float(np.dot(basis[1], basis[2])), 0.0, places=10)

    def test_cut_conserves_atoms_for_all_orientations(self):
        """cut() must not silently drop atoms on the cell boundary."""
        # 2x2x2 conventional FCC cell: 32 atoms.
        a = 4.0
        basis_atoms = np.array([
            [0.0, 0.0, 0.0], [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5], [0.0, 0.5, 0.5],
        ])
        atom_list = []
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for base in basis_atoms:
                        atom_list.append(Atom(
                            index=len(atom_list),
                            symbol="Cu",
                            coord=(base + np.array([i, j, k])) / 2.0,
                        ))
        cell = np.eye(3) * (2 * a)
        struct = Struct(cell=cell, is_direct=True, atom_list=atom_list)

        with tempfile.TemporaryDirectory() as temp_dir:
            poscar_p = Path(temp_dir) / "fcc.vasp"
            SimplePoscar.write_poscar(poscar_p, struct, "FCC 2x2x2")

            for miller in [(0, 0, 1), (1, 1, 0), (1, 1, 1)]:
                s = Slicer("fcc", poscar_p, miller)
                expected = int(round(
                    len(struct) * abs(np.linalg.det(scaled_basis(miller)))
                ))
                self.assertEqual(
                    len(s.transformed), expected,
                    f"{miller}: cut lost atoms ({len(s.transformed)} != {expected})",
                )
                layers = list(s.group_by_normal())
                total = sum(len(layer) for _, layer in layers)
                self.assertEqual(total, expected,
                    f"{miller}: layers do not cover the transformed structure")
