import unittest
import numpy as np
import os
import tempfile
import shutil
from PoscarTools.AtomAllocate import allocate_atoms, allocate2files, _integer_fracts
from PoscarTools.SimplePoscar import Atom, Atoms


class TestAtomAllocate(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures before each test method."""
        # 创建临时目录
        self.test_dir = tempfile.mkdtemp()
        
        # 创建用于测试的原子结构
        cell = np.array([[1.0, 0.0, 0.0], 
                         [0.0, 1.0, 0.0], 
                         [0.0, 0.0, 1.0]])
        
        # 创建带注释的原子列表，模拟超胞结构
        atom_list = [
            Atom(index=0, symbol='Si', coord=np.array([0.0, 0.0, 0.0]), note="1a-Si"),
            Atom(index=1, symbol='Si', coord=np.array([0.5, 0.5, 0.5]), note="1a-Si"),
            Atom(index=2, symbol='Si', coord=np.array([0.0, 0.5, 0.5]), note="1a-Si"),
            Atom(index=3, symbol='Si', coord=np.array([0.5, 0.0, 0.5]), note="1a-Si")
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
4
Direct
0.000000 0.000000 0.000000
0.500000 0.500000 0.500000
0.000000 0.500000 0.500000
0.500000 0.000000 0.500000
""")

    def tearDown(self):
        """Clean up after each test method."""
        # 删除临时目录
        shutil.rmtree(self.test_dir)

    def test_integer_fracts(self):
        """测试将小数形式的分数转换为整数形式的分数"""
        fracts = {"Au": 0.5, "Cu": 0.3, "Ni": 0.2}
        factors = (2, 2, 2)
        multi = 4  # 4原子位点
        
        result = _integer_fracts(fracts, factors, multi)
        
        # 检查总原子数是否正确 (4 * 2 * 2 * 2 = 32)
        self.assertEqual(sum(result.values()), 32)
        
        # 检查各元素的数量是否合理
        self.assertIn("Au", result)
        self.assertIn("Cu", result)
        self.assertIn("Ni", result)

    def test_allocate_atoms_without_fracts(self):
        """测试不指定分数的原子分配"""
        result = allocate_atoms(self.test_atoms, site_fracts=None, seed=42)
        
        # 检查返回结果类型
        self.assertIsInstance(result, dict)
        
        # 检查位点数量
        self.assertEqual(len(result), 1)
        self.assertIn("1a", result)
        
        # 检查原子数量
        allocated_atoms = result["1a"]
        self.assertEqual(len(allocated_atoms), 4)

    def test_allocate_atoms_with_fracts(self):
        """测试指定分数的原子分配"""
        site_fracts = {"1a": {"Au": 2, "Cu": 2}}
        result = allocate_atoms(self.test_atoms, site_fracts=site_fracts, seed=42)
        
        # 检查位点
        self.assertIn("1a", result)
        
        # 检查原子符号是否正确分配
        allocated_atoms = result["1a"]
        symbols = [atom.symbol for atom in allocated_atoms]
        self.assertIn("Au", symbols)
        self.assertIn("Cu", symbols)

    def test_allocate_atoms_reproducibility(self):
        """测试使用相同种子的可重现性"""
        site_fracts = {"1a": {"Au": 2, "Cu": 2}}
        
        result1 = allocate_atoms(self.test_atoms, site_fracts=site_fracts, seed=42)
        result2 = allocate_atoms(self.test_atoms, site_fracts=site_fracts, seed=42)
        
        # 检查两次结果是否相同
        self.assertEqual(len(result1), len(result2))
        self.assertEqual(len(result1["1a"]), len(result2["1a"]))
        
        symbols1 = [atom.symbol for atom in result1["1a"]]
        symbols2 = [atom.symbol for atom in result2["1a"]]
        self.assertEqual(symbols1, symbols2)

    def test_allocate2files(self):
        """测试原子分配并保存到文件"""
        struct_info = {
            "1a": {"sofs": {"Au": 0.5, "Cu": 0.5}}
        }
        
        output_paths = allocate2files(
            filepath=self.test_poscar_path,
            outdir=self.test_dir,
            factors=(1, 1, 1),
            struct_info=struct_info,
            seeds=[42]
        )
        
        # 检查输出文件列表
        self.assertIsInstance(output_paths, list)
        self.assertGreater(len(output_paths), 0)
        
        # 检查文件是否存在
        for path in output_paths:
            self.assertTrue(os.path.exists(path))


if __name__ == '__main__':
    unittest.main()