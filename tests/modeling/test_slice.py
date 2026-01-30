import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.modeling.slice import Slicer
from src.modeling.base import Struct, Atom, SimplePoscar


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
            transformed_file = result_dir / "Transformed(001).vasp"
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
        expected_basis = np.array([[0, 0, -1], [-1, 1, 0], [1, 1, 0]])
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
            transformed_file = result_dir / "Transformed(110).vasp"
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
        expected_basis = np.array([[1, 1, -2], [-1, 1, 0], [1, 1, 1]])
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
            transformed_file = result_dir / "Transformed(111).vasp"
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
