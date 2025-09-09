import unittest
import numpy as np
import os
import tempfile
import shutil
from PoscarTools.AtomSupercell import make_supercell, supercell2file, unitcell2file
from PoscarTools.SimplePoscar import Atom, Atoms


class TestAtomSupercell(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # 创建临时目录
        self.test_dir = tempfile.mkdtemp()
        
        # 创建基本的晶胞结构用于测试
        cell = np.array([[1.0, 0.0, 0.0], 
                         [0.0, 1.0, 0.0], 
                         [0.0, 0.0, 1.0]])
        
        atom_list = [
            Atom(index=0, symbol='Si', coord=np.array([0.0, 0.0, 0.0])),
            Atom(index=1, symbol='Si', coord=np.array([0.5, 0.5, 0.5]))
        ]
        
        self.test_atoms = Atoms(cell=cell, is_direct=True, atom_list=atom_list)
        
        # 创建测试POSCAR文件
        self.test_poscar_path = os.path.join(self.test_dir, "test_POSCAR")
        with open(self.test_poscar_path, 'w') as f:
            f.write("""Test Structure
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Si
2
Direct
0.000000 0.000000 0.000000
0.500000 0.500000 0.500000
""")

    def tearDown(self):
        """Clean up after each test method."""
        # 删除临时目录
        shutil.rmtree(self.test_dir)

    def test_make_supercell_basic(self):
        """测试基本的超胞生成功能"""
        factors = (2, 2, 2)
        supercell_atoms = make_supercell(self.test_atoms, factors)
        
        # 检查原子数量是否正确 (2 * 2 * 2 * 2 = 16)
        self.assertEqual(len(supercell_atoms), 16)
        
        # 检查晶胞是否正确扩展
        expected_cell = np.array([[2.0, 0.0, 0.0], 
                                  [0.0, 2.0, 0.0], 
                                  [0.0, 0.0, 2.0]])
        np.testing.assert_array_almost_equal(supercell_atoms.cell, expected_cell)

    def test_make_supercell_unequal_factors(self):
        """测试不同方向上不同超胞因子"""
        factors = (1, 2, 3)
        supercell_atoms = make_supercell(self.test_atoms, factors)
        
        # 检查原子数量是否正确 (1 * 2 * 3 * 2 = 12)
        self.assertEqual(len(supercell_atoms), 12)
        
        # 检查晶胞是否正确扩展
        expected_cell = np.array([[1.0, 0.0, 0.0], 
                                  [0.0, 2.0, 0.0], 
                                  [0.0, 0.0, 3.0]])
        np.testing.assert_array_almost_equal(supercell_atoms.cell, expected_cell)

    def test_supercell2file(self):
        """测试超胞生成并保存到文件"""
        factors = (2, 2, 2)
        output_path = supercell2file(self.test_poscar_path, self.test_dir, factors)
        
        # 检查输出文件是否存在
        self.assertTrue(os.path.exists(output_path))
        
        # 检查文件名是否正确
        self.assertIn("POSCAR-supercell-2x2x2.vasp", output_path)

    def test_unitcell2file(self):
        """测试单胞结构信息保存到POSCAR文件"""
        struct_info = {
            "cell": [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            "1a": {"atoms": ["Si", [[0.0, 0.0, 0.0]]]}
        }
        
        output_path = unitcell2file(struct_info, self.test_dir)
        
        # 检查输出文件是否存在
        self.assertTrue(os.path.exists(output_path))
        
        # 检查文件名是否正确
        self.assertIn("POSCAR-unitcell.vasp", output_path)

    def test_unitcell2file_diagonal_cell(self):
        """测试对角晶胞的单胞结构信息保存"""
        struct_info = {
            "cell": [10.0, 10.0, 10.0],  # 对角形式
            "1a": {"atoms": ["Si", [[0.0, 0.0, 0.0]]]}
        }
        
        output_path = unitcell2file(struct_info, self.test_dir)
        
        # 检查输出文件是否存在
        self.assertTrue(os.path.exists(output_path))


if __name__ == '__main__':
    unittest.main()