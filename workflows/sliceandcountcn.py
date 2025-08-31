# workflows/sliceandcountcn.py

import logging
import os
import shutil

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

from PoscarTools.AtomSlice import _normalize, _convert, _get_basis, group_by_normal, plot_layer
from PoscarTools.AtomCountCN import detect_cutoff_distance, generate_poscar_map, \
    merge_by_cn, merge_by_symbol, plot_histogram
from PoscarTools.SimplePoscar import SimplePoscar


def slice2file_with_cn(filepath: str, outdir: str, miller_index: tuple[int, int, int]) -> list[str]:
    """
    对POSCAR文件进行切片, 同时对每个切片层进行配位数统计

    Args:
        filepath: 输入POSCAR文件路径
        outdir: 输出目录
        miller_index: 晶面指数，用于切片方向

    Returns:
        List[str]: 所有层的配位数统计结果目录路径列表
    """

    miller_index_str = "".join(str(d) for d in miller_index)

    # 创建输出目录
    dirname = f"{os.path.splitext(os.path.basename(filepath))[0]}"
    outdir = os.path.join(outdir, f"{dirname}-({miller_index_str})-sliced")
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)

    # 读取POSCAR
    atoms = SimplePoscar.read_poscar(filepath)
    symbols_str = "".join(s for s, c in atoms.symbol_count)
    logging.debug(atoms)

    # 获取基向量, 将晶面指数视为法线
    basis = _get_basis(miller_index)
    logging.info(f"基向量: {basis}")

    # 沿基向量转换原子
    basis_n = np.array([_normalize(v) for v in basis])
    new_atoms = _convert(atoms, basis_n)
    output = os.path.join(outdir, f"POSCAR-convert({miller_index_str})-{symbols_str}.vasp")
    comment = f"Convert({miller_index_str})-{symbols_str}"
    SimplePoscar.write_poscar(filepath=output, atoms=new_atoms, comment=comment)

    # 按法线对原子进行分组
    # basis = np.array([(1, 0, 0), (0, 1, 0), (0, 0, 1)])
    layers = [ls for ls in group_by_normal(atoms=new_atoms, basis=basis_n)]
    num_layers = len(layers)
    logging.info(f"找到层数 {num_layers}")

    # 将层保存为POSCAR并绘制层
    l = len(str(num_layers))
    for i, (proj, layer) in enumerate(tqdm(layers, desc="处理原子层",
                                           total=num_layers, ncols=80), start=1):
        logging.debug(f"第 {i:0{l}d} 层 投影={proj:.4f}")
        logging.debug(f"原子层 {layer}")

        # 将层保存到POSCAR文件
        output = os.path.join(outdir, f"POSCAR-convert({miller_index_str})-layer{i:0{l}d}.vasp")
        comment = f"Convert({miller_index_str})-Layer{i:0{l}d}"
        SimplePoscar.write_poscar(filepath=output, atoms=layer, comment=comment)

        # 对切片层进行配位数统计

        # 创建输出目录
        outd = os.path.join(outdir, f"{os.path.splitext(os.path.basename(output))[0]}-cn")
        if os.path.exists(outd):
            shutil.rmtree(outd)
        os.makedirs(outd, exist_ok=True)

        # 检测截断距离
        cut_off = detect_cutoff_distance(new_atoms)
        logging.info(f"自动检测到的截断距离 {cut_off:.3f} Å")

        # 生成POSCAR映射
        cn_file_map, symbol_df, symbol_cn_freq, pair_counts = generate_poscar_map(
            atoms=layer, cut_off=cut_off, outdir=outd)

        # 按元素合并文件
        merge_by_symbol(cn_file_map=cn_file_map, outdir=outd)

        # 按配位数合并文件
        merge_by_cn(cn_file_map=cn_file_map, outdir=outd)

        # 将CN数据写入CSV
        all_df = pd.concat(symbol_df.values(), ignore_index=True)
        output = os.path.join(outd, "cn-counts.csv")
        all_df.to_csv(output, index=False)
        logging.info(f"配位数计数已保存到 {output}")

        # 绘制直方图
        plot_histogram(symbol_cn_freq=symbol_cn_freq, outdir=outd, pair_counts=pair_counts)

        # 按基向量绘制层
        imgname = os.path.join(outdir, f"{comment}.png")
        plot_layer(layer=layer, basis=basis, title=comment, filepath=imgname, pair_counts=pair_counts)

    logging.info(f"结果已保存在 {outdir}")
    return outdir


def merge_layer_cn_results(slice_outdir: str, cn_result_dirs: list[str]) -> str:
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
