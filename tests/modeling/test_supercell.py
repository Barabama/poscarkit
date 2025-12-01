import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.modeling.supercell import make_supercell, unitcell2file, supercell2file
from src.modeling.base import Struct, Atom
from src.modeling.base import SimplePoscar

poscar = """SYSTEM=Au1Cu3-FCC
 1.0000000000000000
    3.7740000000000000    0.0000000000000000    0.0000000000000000
    0.0000000000000000    3.7740000000000000    0.0000000000000000
    0.0000000000000000    0.0000000000000000    3.7740000000000000
  Au  Cu
   1   3
Selective dynamics
Direct
 0.0000000000000000 0.0000000000000000 0.0000000000000000 T T T # 1a-Au-#1
 0.0000000000000000 0.5000000000000000 0.5000000000000000 T T T # 3c-Cu-#1
 0.5000000000000000 0.0000000000000000 0.5000000000000000 T T T # 3c-Cu-#2
 0.5000000000000000 0.5000000000000000 0.0000000000000000 T T T # 3c-Cu-#3
"""

structure_info = {
    "cell": [[3.774, 0.0, 0.0], [0.0, 3.774, 0.0], [0.0, 0.0, 3.774]],
    "1a": {"atoms": ("Au", [[0.0, 0.0, 0.0]]), "sofs": {"V": 1.0}},
    "3c": {
        "atoms": ("Cu", [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
        "sofs": {"Co": 4.444e-1, "Ni": 0.4444, "V": 0.1112},
    },
}


class TestSupercell(unittest.TestCase):
    def test_simple_supercell_workflow(self):
        """Test simple supercell workflow including unitcell creation, supercell generation and file saving."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Step 1: Create unitcell file from structure info
            unitcell_file = unitcell2file(structure_info, temp_path)
            self.assertTrue(unitcell_file.exists(), "Unitcell file should be created")

            # Step 2: Generate supercell from unitcell file
            factors = (2, 2, 2)
            supercell_file = supercell2file(unitcell_file, temp_path, factors)
            self.assertTrue(supercell_file.exists(), "Supercell file should be created")

            # Step 3: Read supercell and check atom count
            # Read the unitcell to determine initial atom count
            unitcell_struct = SimplePoscar.read_poscar(unitcell_file)
            initial_atoms = len(unitcell_struct)

            # Read supercell structure and verify atom count
            supercell_struct = SimplePoscar.read_poscar(supercell_file)
            expected_atoms = initial_atoms * factors[0] * factors[1] * factors[2]
            self.assertEqual(
                len(supercell_struct),
                expected_atoms,
                f"Supercell should have {expected_atoms} atoms",
            )

            # Step 4: Verify basic properties of the supercell
            # Read the supercell file to check it was created correctly
            with open(supercell_file, "r") as f:
                content = f.read()

            # Check that the comment contains supercell information
            self.assertIn(
                "Supercell-2x2x2",
                content,
                "Supercell file should contain correct comment",
            )

            # Check that the file is not empty
            self.assertGreater(len(content), 0, "Supercell file should not be empty")
