# AtomCountCN.py

import logging
import os
import shutil
from collections import defaultdict
from typing import TypedDict, Dict, List, Set, Tuple, Any, Optional

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm

from .SimplePoscar import Atom, Atoms, SimplePoscar
from .Utils import color_map


class CNData(TypedDict):
    center: Atom
    neighbors: list[Atom]
    cn: int


def detect_cutoff_distance(atoms: Atoms, sample_size: int = 1000) -> float:
    """
    基于距离分布分析自动检测第一近邻距离

    Args:
        atoms: Atoms对象
        sample_size: 用于距离计算的原子采样数量, 避免内存问题

    Returns:
        float: 检测到的第一近邻的截断距离
    """
    coords = atoms.cartesian_coords

    # 对大系统使用采样
    if len(coords) > sample_size:
        logging.info(f"进行截断距离检测采样数 {sample_size}/{len(coords)}")
        indices = np.random.choice(len(coords), size=sample_size, replace=False)
        coords = coords[indices]

    # 使用pdist计算成对距离
    distances = pdist(coords)

    # 过滤掉非常小的距离（同一原子或数值误差）
    distances = distances[distances > 0.1]

    if len(distances) == 0:
        raise ValueError("未找到有效距离")

    # 排序距离以找到第一组相似距离（第一近邻）
    sorted_distances = np.sort(distances)

    # 对于小系统, 使用最小距离
    if len(sorted_distances) < 100:
        return float(sorted_distances[0])

    # 对于较大系统, 使用更复杂的方法
    # 查看前5%的距离
    subset_size = min(len(sorted_distances) // 20, 500)
    subset_distances = sorted_distances[:subset_size]

    # 计算梯度以找到增加率显著变化的位置
    diffs = np.diff(subset_distances)

    # 寻找距离增加率显著增加的第一个位置
    # 这表明从第一近邻到更高近邻的过渡
    threshold = np.mean(diffs) + np.std(diffs)  # 显著增加的阈值

    # 找到第一个与下一个距离差距显著较大的距离
    cutoff_index = 0
    for i in range(len(diffs)):
        if diffs[i] > threshold:
            cutoff_index = i
            break

    # 如果未找到显著差距, 则使用统计方法
    if cutoff_index == 0:
        # 使用最小距离或小百分位数
        cutoff_distance = np.percentile(subset_distances, 5)
    else:
        cutoff_distance = subset_distances[cutoff_index]

    return float(cutoff_distance)


def calculate_nearest_neighbors(atoms: Atoms, cut_off: float) -> Tuple[Dict[Atom, List[Atom]], int]:
    """
    给定的截断距离内计算每个原子的最近邻及其配位数, 使用KDTree提高效率并减少内存使用

    Args:
        atoms: 包含原子位置的Atoms对象
        cut_off: 用于确定最近邻的截断距离

    Returns:
        tuple: 包含原子到其邻居映射和对数的元组
    """
    nn_map: Dict[Atom, List[Atom]] = defaultdict(list)
    pair_count = 0

    # 获取笛卡尔坐标
    coords = atoms.cartesian_coords
    # 构建KDTree
    tree = KDTree(coords)
    for i, coord in enumerate(tqdm(coords, desc="搜索最近邻", ncols=80)):
        # 在rqstd + 容差范围内搜索邻居
        neighbor_indices = tree.query_ball_point(coord, r=1.2 * cut_off)

        # 计算距离
        distances = np.linalg.norm(coords[neighbor_indices] - coord, axis=1)

        # 过滤邻居
        atom_i = atoms[i]
        tolerance = 0.1 * cut_off
        for j, dist in zip(neighbor_indices, distances):
            if dist < tolerance:
                continue
            if (dist - cut_off) > tolerance:
                continue

            atom_j = atoms[j]
            nn_map[atom_i].append(atom_j)
            if atom_j not in nn_map:
                pair_count += 1

    return nn_map, pair_count


def generate_poscar_map(atoms: Atoms, cut_off: float, outdir: str
                        ) -> Tuple[Dict[Tuple[str, int], str], Dict[str, Any], Dict[str, List[int]], Dict[str, int]]:
    """
    生成POSCAR文件映射关系

    Args:
        atoms: 包含原子位置的Atoms对象
        cut_off: 用于确定最近邻的截断距离
        outdir: 生成文件的输出目录

    Returns:
        tuple: 包含以下内容的元组：
        - (元素, 配位数) 到文件路径的映射
        - 用于生成CSV的元素数据字典
        - 用于直方图的元素配位数频率
        - 元素到对数的映射
    """
    cn_file_map = {}  # map {(symbol, cn): filepath}
    symbol_df = {}
    symbol_cn_freq = {}
    symbol_pair_counts = {}  # 新增：存储每个元素的pair_count

    for symbol, subatoms in atoms.group_atoms():
        # 获取最近邻
        nn_map, pair_count = calculate_nearest_neighbors(subatoms, cut_off)
        logging.info(f"{symbol}-{symbol} 对数 {pair_count}")
        symbol_pair_counts[symbol] = pair_count  # 保存pair_count

        # 收集CN数据映射
        cn_data_map = defaultdict(list)
        for center, neighbors in nn_map.items():
            cn = len(neighbors)
            cn_data_map[cn].append(CNData(center=center, neighbors=neighbors, cn=cn))

        # 获取每个配位数的原子并生成POSCAR文件
        # nn_atoms_set: Set[Atom] = set()
        for cn, cn_data_list in cn_data_map.items():
            cn_atom_sets: Set[Atom] = set()
            for cn_data in cn_data_list:
                sets = {cn_data["center"], *cn_data["neighbors"]}
                cn_atom_sets.update(sets)

            # nn_atoms_set.update(cn_atom_sets)
            cn_atoms = atoms.copy(atom_list=list(cn_atom_sets))
            logging.debug(f"{symbol}*-{cn}{symbol} 配位原子 {cn_atoms}")

            # 保存到文件
            output = os.path.join(outdir, f"POSCAR-d1nn-{symbol}-{cn}.vasp")
            comment = f"CoordinationNumber-{symbol}-{cn}"
            SimplePoscar.write_poscar(filepath=output, atoms=cn_atoms, comment=comment)

            # 添加到映射
            cn_file_map[(symbol, cn)] = output

        # # 获取最近邻原子
        # nn_atoms = atoms.copy(atom_list=list(nn_atoms_set))
        # logging.debug(f"最近邻原子: {nn_atoms}")

        # # 保存到文件
        # output = os.path.join(outdir, f"POSCAR-d1nn-{symbol}.vasp")
        # comment = f"最近邻-{symbol}-对数={pair_count}"
        # SimplePoscar.write_poscar(filepath=output, atoms=nn_atoms, comment=comment)
        
        # 收集原子对数据
        
        
        # 收集CN数据到df
        cn_counts = [[f"{symbol}*-{cn}{symbol}", len(d)] for cn, d in cn_data_map.items()]
        cn_df = pd.DataFrame(data=cn_counts, columns=["CN", "Count"])
        symbol_df[symbol] = cn_df

        # 收集CN数据到列表
        cn_freq = [d["cn"] for ds in cn_data_map.values() for d in ds]
        symbol_cn_freq[symbol] = cn_freq

    return cn_file_map, symbol_df, symbol_cn_freq, symbol_pair_counts


def merge_by_symbol(cn_file_map: Dict[Tuple[str, int], str], outdir: str):
    """
    按元素合并POSCAR文件

    Args:
        cn_file_map: (元素, 配位数) 到文件路径的映射
        outdir: 生成文件的输出目录
    """
    symbol_files = defaultdict(list)
    for (symbol, cn), filepaths in cn_file_map.items():
        symbol_files[symbol].append(filepaths)

    # 为每个元素合并CN文件
    for symbol, filepaths in symbol_files.items():
        if len(filepaths) == 1:
            # 如果只有一个文件, 则无需合并
            continue

        # 合并该元素的所有CN文件
        merged_atoms: Optional[Atoms] = None
        for filepath in filepaths:
            add_atoms = SimplePoscar.read_poscar(filepath)
            if merged_atoms is None:
                merged_atoms = add_atoms.copy()
            else:
                merged_atoms.extend(add_atoms)

        # 将合并的原子保存到文件
        output = os.path.join(outdir, f"POSCAR-d1nn-{symbol}.vasp")
        comment = f"All coordination numbers for symbol {symbol}"
        if merged_atoms is not None:
            SimplePoscar.write_poscar(filepath=output, atoms=merged_atoms, comment=comment)


def merge_by_cn(cn_file_map: Dict[Tuple[str, int], str], outdir: str):
    """
    按配位数合并POSCAR文件

    Args:
        cn_file_map: (元素, 配位数) 到文件路径的映射
        outdir: 生成文件的输出目录
    """
    cn_files = defaultdict(list)
    for (symbol, cn), filepaths in cn_file_map.items():
        cn_files[cn].append(filepaths)

    # 为每个配位数合并POSCAR文件
    for cn, filepaths in cn_files.items():
        if len(filepaths) == 1:
            # 如果只有一个文件, 则无需合并
            continue

        # 合并该配位数的所有CN文件
        merged_atoms: Optional[Atoms] = None
        for filpath in filepaths:
            add_atoms = SimplePoscar.read_poscar(filpath)
            if merged_atoms is None:
                merged_atoms = add_atoms.copy()
            else:
                merged_atoms.extend(add_atoms)

        # 将合并的原子保存到文件
        output = os.path.join(outdir, f"POSCAR-d1nn-all-{cn}.vasp")
        comment = f"All atoms with coordination number {cn}"
        if merged_atoms is not None:
            SimplePoscar.write_poscar(filepath=output, atoms=merged_atoms, comment=comment)


def plot_histogram(symbol_cn_freq: Dict[str, List[int]], outdir: str, symbol_pair_counts: Optional[Dict[str, int]] = None):
    """
    绘制配位数直方图

    Args:
        symbol_cn_freq: 将元素映射到其配位数频率的字典
        outdir: 生成文件的输出目录
        symbol_pair_counts: 可选，将元素映射到其原子对数的字典
    """
    # 绘制配位数直方图
    symbols, cn_freqs = zip(*symbol_cn_freq.items())
    colors = [color_map[s] for s in symbols]
    
    # 准备图例标签，如果提供了pair_counts则包含在内
    if symbol_pair_counts:
        labels = [f"{symbol} (pairs: {symbol_pair_counts.get(symbol, 0)})" for symbol in symbols]
    else:
        labels = list(symbols)
    
    plt.figure(figsize=(10, 6))
    plt.hist(cn_freqs, bins=list(range(0, 12)), alpha=0.5, label=labels, color=colors, edgecolor="black")
    plt.legend()
    plt.title("Histogram of Coordination Numbers")
    plt.xlabel("Coordinate Number")
    plt.ylabel("Frequency")
    plt.xticks(range(0, 12))  # Define ticks
    plt.grid(True, linestyle='--', alpha=0.7)
    output = os.path.join(outdir, "cn-histogram.png")
    plt.savefig(output)
    plt.close()
    logging.info(f"配位数直方图已保存到 {output}")


def countCN2files(filepath: str, outdir: str) -> str:
    """
    计算配位数并为每种配位类型生成POSCAR文件

    Args:
        filepath: 输入POSCAR文件的路径
        outdir: 生成文件的输出目录

    Returns:
        包含输出文件的目录路径
    """
    # 读取POSCAR
    atoms = SimplePoscar.read_poscar(filepath)
    logging.debug(f"原子 {atoms}")

    # 检测截断距离
    cut_off = detect_cutoff_distance(atoms)
    logging.info(f"自动检测到的截断距离 {cut_off:.3f} Å")

    # 创建输出目录
    outdir = os.path.join(outdir, f"{os.path.splitext(os.path.basename(filepath))[0]}-cn")
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)

    # 生成POSCAR映射
    cn_file_map, symbol_df, symbol_cn_freq, symbol_pair_counts = generate_poscar_map(
        atoms=atoms, cut_off=cut_off, outdir=outdir)

    # 按元素合并文件
    merge_by_symbol(cn_file_map=cn_file_map, outdir=outdir)

    # 按配位数合并文件
    merge_by_cn(cn_file_map=cn_file_map, outdir=outdir)

    # 将CN数据写入CSV
    all_df = pd.concat(symbol_df.values(), ignore_index=True)
    output = os.path.join(outdir, "cn-counts.csv")
    all_df.to_csv(output, index=False)
    logging.info(f"配位数计数已保存到 {output}")

    # 绘制直方图
    plot_histogram(symbol_cn_freq, outdir, symbol_pair_counts)

    return outdir
