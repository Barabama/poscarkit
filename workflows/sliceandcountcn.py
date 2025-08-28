# workflows/sliceandcountcn.py

import os
import logging
from typing import List

from PoscarTools.AtomSlice import slice2file
from PoscarTools.AtomCountCN import countCN2files
from PoscarTools.SimplePoscar import SimplePoscar


def slice_and_count_cn(filepath: str, outdir: str, miller_index: tuple[int, int, int]) -> List[str]:
    """
    对POSCAR文件进行切片, 然后对每个切片层进行配位数统计
    
    Args:
        filepath: 输入POSCAR文件路径
        outdir: 输出目录
        miller_index: 晶面指数，用于切片方向
        
    Returns:
        List[str]: 所有层的配位数统计结果目录路径列表
    """
    # 先进行切片操作
    slice_outdir = slice2file(filepath=filepath, outdir=outdir, miller_index=miller_index)
    
    # 查找所有切片层文件
    layer_files = []
    for file in os.listdir(slice_outdir):
        if file.startswith("POSCAR-convert") and "layer" in file and file.endswith(".vasp"):
            layer_files.append(os.path.join(slice_outdir, file))
    
    # 对每个切片层进行配位数统计
    cn_result_dirs = []
    for layer_file in sorted(layer_files):
        logging.info(f"Processing layer file: {layer_file}")
        
        # 对该层进行配位数统计
        cn_outdir = countCN2files(filepath=layer_file, outdir=slice_outdir)
        cn_result_dirs.append(cn_outdir)
        
    logging.info(f"Completed slice and count CN for {len(layer_files)} layers")
    return cn_result_dirs


def merge_layer_cn_results(slice_outdir: str, cn_result_dirs: List[str]) -> str:
    """
    合并所有层的配位数统计结果
    
    Args:
        slice_outdir: 切片操作的输出目录
        cn_result_dirs: 各层配位数统计结果目录列表
        
    Returns:
        str: 合并结果的目录路径
    """
    # 创建合并结果目录
    merged_outdir = os.path.join(slice_outdir, "merged_cn_results")
    os.makedirs(merged_outdir, exist_ok=True)
    
    # TODO: 实现合并逻辑，例如：
    # 1. 合并所有层的CSV文件
    # 2. 生成汇总的直方图
    # 3. 生成各层配位数对比图
    
    return merged_outdir