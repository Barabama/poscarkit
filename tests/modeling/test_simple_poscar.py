import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.modeling.base import Struct, Atom, SimplePoscar


poscar1 = """SYSTEM=Structure1
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

poscar2 = """SYSTEM=Structure2
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

poscar3 = """SYSTEM=Structure3
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
 0.5000000000000000 0.5000000000000000 0.2500000000000000 T T T # 3c-Cu-#3
"""

poscar4 = """SYSTEM=Structure4
 1.0000000000000000
    3.7740000000000000    0.0000000000000000    0.0000000000000000
    0.0000000000000000    3.7740000000000000    0.0000000000000000
    0.0000000000000000    0.0000000000000000    3.7740000000000000
  Au  Cu
   2   2
Selective dynamics
Direct
 0.0000000000000000 0.0000000000000000 0.0000000000000000 T T T # 1a-Au-#1
 0.0000000000000000 0.5000000000000000 0.5000000000000000 T T T # 3c-Cu-#1
 0.5000000000000000 0.0000000000000000 0.5000000000000000 T T T # 3c-Cu-#2
 0.5000000000000000 0.5000000000000000 0.0000000000000000 T T T # 1a-Au-#2
"""

poscar5 = """SYSTEM=Structure5
 1.0000000000000000
    3.7740000000000000    0.0000000000000000    0.0000000000000000
    0.0000000000000000    3.7740000000000000    0.0000000000000000
    0.0000000000000000    0.0000000000000000    3.7740000000000000
  Au  Cu
   2   2
Selective dynamics
Direct
 0.0000000000000000 0.0000000000000000 0.0000000000000000 T T T # 1a-Au-#1
 0.0000000000000000 0.5000000000000000 0.5000000000000000 T T T # 3c-Cu-#1
 0.5000000000000000 0.0000000000000000 0.5000000000000000 T T T # 3c-Cu-#2
 0.5000000000000000 0.5000000000000000 0.0000000000000000 T T T # 3c-Cu-#2
"""


class TestSimplePoscar(unittest.TestCase):
    def test_compare_poscar_identical(self):
        """Test comparing two identical POSCAR files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            poscar1_file = temp_path / "POSCAR1.vasp"
            poscar2_file = temp_path / "POSCAR2.vasp"
            
            poscar1_file.write_text(poscar1)
            poscar2_file.write_text(poscar2)
            
            struct1 = SimplePoscar.read_poscar(poscar1_file)
            struct2 = SimplePoscar.read_poscar(poscar2_file)
            
            flag, msg = struct1.compare(struct2)
            self.assertTrue(flag, "Identical structures should be equal")
    
    def test_compare_poscar_different_coords(self):
        """Test comparing two POSCAR files with different coordinates."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            poscar1_file = temp_path / "POSCAR1.vasp"
            poscar3_file = temp_path / "POSCAR3.vasp"
            
            poscar1_file.write_text(poscar1)
            poscar3_file.write_text(poscar3)
            
            struct1 = SimplePoscar.read_poscar(poscar1_file)
            struct3 = SimplePoscar.read_poscar(poscar3_file)
            
            flag, msg = struct1.compare(struct3)
            self.assertFalse(flag, "Structures with different coordinates should not be equal")
    
    def test_compare_poscar_different_symbols(self):
        """Test comparing two POSCAR files with different symbol counts."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            poscar1_file = temp_path / "POSCAR1.vasp"
            poscar4_file = temp_path / "POSCAR4.vasp"
            
            poscar1_file.write_text(poscar1)
            poscar4_file.write_text(poscar4)
            
            struct1 = SimplePoscar.read_poscar(poscar1_file)
            struct4 = SimplePoscar.read_poscar(poscar4_file)
            
            flag, msg = struct1.compare(struct4)
            self.assertFalse(flag, "Structures with different symbol counts should not be equal")
    
    def test_merge_poscar(self):
        """Test merging two POSCAR files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # 创建两个具有不同坐标的 POSCAR 文件
            # POSCAR 1: 原始结构
            poscar1_content = """SYSTEM=Structure1
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
            
            # POSCAR 2: 具有完全不同坐标的结构
            poscar2_content = """SYSTEM=Structure2
 1.0000000000000000
    3.7740000000000000    0.0000000000000000    0.0000000000000000
    0.0000000000000000    3.7740000000000000    0.0000000000000000
    0.0000000000000000    0.0000000000000000    3.7740000000000000
  Ag  Ni
   2   2
Selective dynamics
Direct
 0.2500000000000000 0.2500000000000000 0.2500000000000000 T T T # 1a-Ag-#1
 0.2500000000000000 0.7500000000000000 0.7500000000000000 T T T # 3c-Ni-#1
 0.7500000000000000 0.2500000000000000 0.7500000000000000 T T T # 3c-Ni-#2
 0.7500000000000000 0.7500000000000000 0.2500000000000000 T T T # 1a-Ag-#2
"""
            
            poscar1_file = temp_path / "POSCAR1.vasp"
            poscar2_file = temp_path / "POSCAR2.vasp"
            
            poscar1_file.write_text(poscar1_content)
            poscar2_file.write_text(poscar2_content)
            
            struct1 = SimplePoscar.read_poscar(poscar1_file)
            struct2 = SimplePoscar.read_poscar(poscar2_file)
            
            initial_atoms = len(struct1) + len(struct2)
            
            merged_file = SimplePoscar.merge_poscar(poscar1_file, poscar2_file, temp_path)
            
            expected_merged_file = temp_path / "POSCAR-merged-POSCAR1-POSCAR2.vasp"
            self.assertEqual(merged_file, expected_merged_file, "Returned file path should match expected path")
            self.assertTrue(merged_file.exists(), "Merged file should be created")
            
            merged_struct = SimplePoscar.read_poscar(merged_file)
            self.assertEqual(
                len(merged_struct),
                initial_atoms,
                f"Merged structure should have {initial_atoms} atoms"
            )
    
    def test_separate_poscar_by_note(self):
        """Test separating a POSCAR file by note."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            poscar1_file = temp_path / "POSCAR1.vasp"
            poscar1_file.write_text(poscar1)
            
            struct = SimplePoscar.read_poscar(poscar1_file)
            initial_atoms = len(struct)
            
            separated_files = SimplePoscar.separate_poscar(poscar1_file, temp_path, key="note")
            
            group1_file = temp_path / "POSCAR-group-1a-Au.vasp"
            group2_file = temp_path / "POSCAR-group-3c-Cu.vasp"
            
            self.assertEqual(len(separated_files), 2, "Should return exactly 2 separated files")
            self.assertIn(group1_file, separated_files, "Group 1a-Au file should be in returned list")
            self.assertIn(group2_file, separated_files, "Group 3c-Cu file should be in returned list")
            for file in separated_files:
                self.assertTrue(file.exists(), f"Returned file {file} should exist")
            
            self.assertTrue(group1_file.exists(), "Group 1a-Au file should be created")
            self.assertTrue(group2_file.exists(), "Group 3c-Cu file should be created")
            
            group1_struct = SimplePoscar.read_poscar(group1_file)
            group2_struct = SimplePoscar.read_poscar(group2_file)
            
            total_separated = len(group1_struct) + len(group2_struct)
            self.assertEqual(
                total_separated,
                initial_atoms,
                f"Total separated atoms should be {initial_atoms}"
            )
    
    def test_separate_poscar_by_symbol(self):
        """Test separating a POSCAR file by symbol."""
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            poscar1_file = temp_path / "POSCAR1.vasp"
            poscar1_file.write_text(poscar1)
            
            struct = SimplePoscar.read_poscar(poscar1_file)
            initial_atoms = len(struct)
            
            separated_files = SimplePoscar.separate_poscar(poscar1_file, temp_path, key="symbol")
            
            au_file = temp_path / "POSCAR-group-Au.vasp"
            cu_file = temp_path / "POSCAR-group-Cu.vasp"
            
            self.assertEqual(len(separated_files), 2, "Should return exactly 2 separated files")
            self.assertIn(au_file, separated_files, "Au group file should be in returned list")
            self.assertIn(cu_file, separated_files, "Cu group file should be in returned list")
            for file in separated_files:
                self.assertTrue(file.exists(), f"Returned file {file} should exist")
            
            self.assertTrue(au_file.exists(), "Au group file should be created")
            self.assertTrue(cu_file.exists(), "Cu group file should be created")
            
            au_struct = SimplePoscar.read_poscar(au_file)
            cu_struct = SimplePoscar.read_poscar(cu_file)
            
            self.assertEqual(len(au_struct), 1, "Au group should have 1 atom")
            self.assertEqual(len(cu_struct), 3, "Cu group should have 3 atoms")
            
            total_separated = len(au_struct) + len(cu_struct)
            self.assertEqual(
                total_separated,
                initial_atoms,
                f"Total separated atoms should be {initial_atoms}"
            )
