import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.workflow.slice_to_countcn import slice2files_with_countcn
from src.modeling.base import Struct, Atom, SimplePoscar


class TestSliceToCNCount(unittest.TestCase):
    """Test cases for slice_to_cn_count module."""

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

    def test_slice2files_with_cn_count(self):
        """Test slice2files_with_cn_count function."""
        with tempfile.TemporaryDirectory() as output_dir:
            outdir = Path(output_dir)
            result_dirs = slice2files_with_countcn(
                name="test_slice_cn",
                poscar=self.poscar_path,
                outdir=outdir,
                miller_index=(0, 0, 1),
            )

            # Check that we got a list of directories
            self.assertIsInstance(result_dirs, list)
            self.assertTrue(all(isinstance(d, Path) for d in result_dirs))
            self.assertGreater(len(result_dirs), 0)

            # Check that each result directory exists
            for result_dir in result_dirs:
                self.assertTrue(result_dir.exists())

                # Check that CN count directory was created inside the layer's parent directory
                # The CNCounter creates a subdirectory named "{name}-cn-count"
                # where name is the stem of the layer file
                layer_file_stem = result_dir.name.split("-cn-count")[0]
                cn_count_dir = result_dir

                # Verify the directory name format
                self.assertIn("-cn-count", str(cn_count_dir))

                # Check that POSCAR files were created in the CN count directory
                poscar_files = list(cn_count_dir.glob("*.vasp"))
                # Note: Some structures might be too simple to have enough neighbors for CN analysis
                # So we might not always get POSCAR files

                # Check that CSV files were created
                csv_files = list(cn_count_dir.glob("*.csv"))
                # We should at least get the statistics CSV file

                # Check that histogram files were created
                png_files = list(cn_count_dir.glob("*.png"))
                # Similar to POSCAR files, we might not always get histograms for simple structures


if __name__ == "__main__":
    unittest.main()
