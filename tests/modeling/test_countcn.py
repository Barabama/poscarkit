import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.modeling.countcn import CNCounter
from src.modeling.base import Struct, Atom, SimplePoscar


class TestCNCounter(unittest.TestCase):
    """Test cases for CNCounter class."""

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
        
        self.counter = CNCounter(name="test", poscar=self.poscar_path)

    def test_countCN2files(self):
        """Test countCN2files method."""
        # Run the countCN2files method with by_ase=False (default)
        output_dir = self.counter.countCN2files(outdir=self.temp_dir, cutoff_mult=1.1, parallel=2, by_ase=False)
        
        # Check that output directory was created
        self.assertTrue(output_dir.exists())
        self.assertTrue(output_dir.is_dir())
        
        # Check that some output files were created
        csv_file = output_dir / "test-d1nn-cn-count.csv"
        self.assertTrue(csv_file.exists())
        
        png_files = list(output_dir.glob("*.png"))
        self.assertGreater(len(png_files), 0)
        
        poscar_files = list(output_dir.glob("*.vasp"))
        # Note: depending on the structure, there may or may not be POSCAR files generated
        
    def test_countCN2files_ase(self):
        """Test countCN2files method with ASE."""
        # Run the countCN2files method with by_ase=True
        output_dir = self.counter.countCN2files(outdir=self.temp_dir, cutoff_mult=1.1, parallel=2, by_ase=True)
        
        # Check that output directory was created
        self.assertTrue(output_dir.exists())
        self.assertTrue(output_dir.is_dir())
        
        # Check that some output files were created
        csv_file = output_dir / "test-d1nn-cn-count.csv"
        self.assertTrue(csv_file.exists())
        
    def tearDown(self):
        """Clean up test fixtures."""
        # Clean up temporary directory
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()