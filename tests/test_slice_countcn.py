# tests/test_slice_countcn.py

from workflows.sliceandcountcn import slice2file_with_cn
import os
import sys
import unittest
import tempfile
import shutil

# 添加项目根目录到Python路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


class TestSliceAndCountCN(unittest.TestCase):

    def setUp(self):
        """测试前准备"""
        # 创建临时目录用于测试
        self.test_dir = tempfile.mkdtemp()

        # 创建一个简单的POSCAR测试文件（仅用于演示）
        self.test_poscar = os.path.join(self.test_dir, "POSCAR")
        with open(self.test_poscar, 'w') as f:
            f.write("""Test Structure
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
Si
8
Direct
0.000000 0.000000 0.000000
0.250000 0.250000 0.250000
0.500000 0.500000 0.500000
0.750000 0.750000 0.750000
0.000000 0.500000 0.500000
0.500000 0.000000 0.500000
0.500000 0.500000 0.000000
0.250000 0.250000 0.750000
""")

    def tearDown(self):
        """测试后清理"""
        # 删除临时目录
        shutil.rmtree(self.test_dir)

    def test_slice_and_count_cn_creation(self):
        """测试slice_and_count_cn函数是否能正常导入和创建"""
        # 注意：由于需要真实的晶体结构文件和复杂的计算，这里仅测试函数是否存在
        self.assertTrue(callable(slice2file_with_cn))

    def test_workflow_structure(self):
        """测试工作流模块结构"""
        # 检查模块是否正确导入
        try:
            from workflows.sliceandcountcn import slice2file_with_cn, merge_layer_cn_results
            self.assertTrue(callable(slice2file_with_cn))
            self.assertTrue(callable(merge_layer_cn_results))
        except ImportError:
            self.fail("无法从workflows.sliceandcountcn导入所需函数")


if __name__ == '__main__':
    unittest.main()
