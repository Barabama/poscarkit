import unittest
import sys
import os

# 添加项目根目录到Python路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

def create_test_suite():
    """创建测试套件，包含所有测试"""
    suite = unittest.TestSuite()
    
    # 导入并添加各个测试模块
    try:
        from tests.PoscarTools.test_atom_supercell import TestAtomSupercell
        suite.addTest(unittest.makeSuite(TestAtomSupercell))
    except ImportError as e:
        print(f"警告: 无法导入TestAtomSupercell: {e}")
    
    try:
        from tests.PoscarTools.test_atom_allocate import TestAtomAllocate
        suite.addTest(unittest.makeSuite(TestAtomAllocate))
    except ImportError as e:
        print(f"警告: 无法导入TestAtomAllocate: {e}")
    
    try:
        from tests.PoscarTools.test_atom_slice import TestAtomSlice
        suite.addTest(unittest.makeSuite(TestAtomSlice))
    except ImportError as e:
        print(f"警告: 无法导入TestAtomSlice: {e}")
    
    try:
        from tests.test_slice_countcn import TestSliceAndCountCN
        suite.addTest(unittest.makeSuite(TestSliceAndCountCN))
    except ImportError as e:
        print(f"警告: 无法导入TestSliceAndCountCN: {e}")
        
    try:
        from tests.workflows.test_modeling import TestModelingWorkflow
        suite.addTest(unittest.makeSuite(TestModelingWorkflow))
    except ImportError as e:
        print(f"警告: 无法导入TestModelingWorkflow: {e}")
        
    try:
        from tests.workflows.test_sliceandcountcn import TestSliceAndCountCNWorkflow
        suite.addTest(unittest.makeSuite(TestSliceAndCountCNWorkflow))
    except ImportError as e:
        print(f"警告: 无法导入TestSliceAndCountCNWorkflow: {e}")
        
    # 移除重复的TestSliceAndCountCN导入
    
    return suite


def run_all_tests():
    """运行所有测试"""
    suite = create_test_suite()
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_all_tests()
    # 退出码表示测试是否全部通过
    exit(0 if success else 1)