# AtomCountCN.py

import logging
import os
import shutil
from collections import defaultdict
from dataclasses import dataclass
from typing import Any, TypedDict

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm

from .SimplePoscar import Atom, Atoms, SimplePoscar
from .Utils import color_map

_hatch_patterns = ["//", "\\\\", "||", "--", "++", "xx", "oo", "O)", "..", "**"]


@dataclass
class CNData:
    symbols: tuple[str, str]  # 中心元素, 邻居元素
    center: Atom
    neighbors: list[Atom]
    cn: int     # 配位数


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


def calculate_nearest_neighbors(atoms: Atoms, cut_off: float
                                ) -> tuple[list[CNData], dict[frozenset, int]]:
    """
    给定的截断距离内计算每个原子的最近邻及其配位数, 使用KDTree提高效率并减少内存使用

    Args:
        atoms: 包含原子位置的Atoms对象
        cut_off: 用于确定最近邻的截断距离

    Returns:
        tuple: 
            - cndata_list: list[CNData], 配位数据列表
            - pair_counts: {pair: count}, 每种原子对的数量
    """
    nn_map: dict[Atom, dict[str, list[Atom]]] = defaultdict(lambda: defaultdict(list))
    pair_counts: dict[frozenset, int] = defaultdict(int)

    coords = atoms.cartesian_coords  # 获取笛卡尔坐标
    tree = KDTree(coords)  # 构建KDTree
    tolerance = 0.1 * cut_off

    for i, coord in enumerate(tqdm(coords, desc="搜索最近邻", ncols=80)):
        atom_i = atoms[i]
        symbol_i = atom_i.symbol

        # 在 cut_off + tolerance 范围内搜索邻居
        neighbor_indices = tree.query_ball_point(coord, r=1.2 * cut_off)
        for j in neighbor_indices:
            if i == j:  # 跳过自己
                continue

            dist = np.linalg.norm(coords[j] - coord)
            if dist < tolerance or dist > (cut_off + tolerance):
                # 跳过超出范围的邻居
                continue

            atom_j = atoms[j]
            symbol_j = atom_j.symbol

            # 邻居j添加到中心i的邻居列表中
            nn_map[atom_i][symbol_j].append(atom_j)

            # 更新对数统计
            pair_key = frozenset([symbol_i, symbol_j])
            pair_counts[pair_key] += 1

    # defaultdict转换为dict
    cndata_list = []
    for atom_ct, nn_data in nn_map.items():
        s_ct = atom_ct.symbol
        for s_nb, nn_list in nn_data.items():
            if len(nn_list) <= 0:
                continue
            cndata_list.append(CNData(symbols=(s_ct, s_nb),
                                      center=atom_ct,
                                      neighbors=nn_list,
                                      cn=len(nn_list)))
    pair_counts = dict(pair_counts)

    return cndata_list, pair_counts


def generate_poscar(atoms: Atoms, cndata_list: list[CNData], outdir: str) -> dict[tuple[str, str, int], str]:
    """
    生成每个(中心, 邻居, 配位数)的POSCAR文件

    Args:
        atoms: Atoms对象
        cndata_list: list[CNData] CNData对象列表
        outdir: 输出目录

    Returns:
        dict: {(symbol_ct, symbol_nb, cn): filepath, ...} 文件路径字典
    """
    # 收集相关原子, 中心和邻居
    cndata_dict = defaultdict(set)
    for cndata in cndata_list:
        s_ct, s_nb = cndata.symbols
        cn = cndata.cn
        cndata_dict[(s_ct, s_nb, cn)].add(cndata.center)
        cndata_dict[(s_ct, s_nb, cn)].update(cndata.neighbors)

    cn_file_map = {}
    for (s_ct, s_nb, cn), subatom_set in cndata_dict.items():
        # 创建子结构
        subatoms = atoms.copy(atom_list=list(subatom_set))
        logging.info(f"{s_ct}*-{cn}{s_nb} 配位原子 {subatoms}")

        # 保存POSCAR文件
        output = os.path.join(outdir, f"POSCAR-d1nn-{s_ct}-{s_nb}-cn{cn}.vasp")
        comment = f"CoordinationNumber-{s_ct}-{s_nb}-{cn}"
        SimplePoscar.write_poscar(filepath=output, atoms=subatoms, comment=comment)
        cn_file_map[(s_ct, s_nb, cn)] = output

    return cn_file_map


def merge_by_cn(cn_file_map: dict[tuple[str, str, int], str], outdir: str):

    cn_files = defaultdict(list)
    for (s_ct, s_nb, cn), filepath in cn_file_map.items():
        cn_files[cn].append(filepath)

    for cn, filepaths in cn_files.items():
        if len(filepaths) == 1:
            # 单个文件无需合并
            continue

        atoms_merged: Atoms | None = None
        for filepath in filepaths:
            atoms_add = SimplePoscar.read_poscar(filepath)
            if atoms_merged is None:
                atoms_merged = atoms_add
            else:
                atoms_merged.extend(atoms_add)
        # 保存合并后的POSCAR文件
        output = os.path.join(outdir, f"POSCAR-d1nn-all-cn{cn}.vasp")
        comment = f"All atoms with coordination number {cn}"
        if atoms_merged is not None:
            SimplePoscar.write_poscar(filepath=output, atoms=atoms_merged, comment=comment)


def save_dataframe(cndata_list: list[CNData], outdir: str):

    # 分组配位数据
    cndata_dict = defaultdict(list)
    for cndata in cndata_list:
        s_ct, s_nb = cndata.symbols
        cn = cndata.cn
        cndata_dict[(s_ct, s_nb, cn)].append(cndata)
    cndata_dict = {k: v for k, v in sorted(cndata_dict.items(), key=lambda x: x[0])}

    # 配位数据保存CSV文件
    data = defaultdict(list)
    for (s_ct, s_nb, cn), cndata_list in cndata_dict.items():
        data["CN"].append(f"{s_ct}*-{cn}{s_nb}")
        data["Count"].append(len(cndata_list))

    output = os.path.join(outdir, "cn-count.csv")
    df = pd.DataFrame(data)
    df.to_csv(output, index=False)
    logging.info(f"配位数据已保存到 {output}")


def plot_histogram_faceted(cndata_list: list[CNData], pair_counts: dict[frozenset, int], outdir: str):

    symbols = sorted(list(set(d.symbols[0] for d in cndata_list)))
    n = len(symbols)

    # 分组配位数据
    cn_stats = defaultdict(list)
    for cndata in cndata_list:
        cn_stats[cndata.symbols].append(cndata.cn)

    symbol2hatch = {symbol: _hatch_patterns[i % len(_hatch_patterns)] for i, symbol in enumerate(symbols)}

    # n x n 子图网格
    fig, axes = plt.subplots(n, n, figsize=(4*n, 4*n), sharex=True, sharey=True)
    if n == 1:
        axes = np.array([axes])  # 兼容单图

    for i, s_ct in enumerate(symbols):
        for j, s_nb in enumerate(symbols):
            ax = axes[i, j]
            key = (s_ct, s_nb)
            data = cn_stats.get(key, [])
            if len(data) <= 0:
                continue

            ax.hist(data, bins=range(0, max(data) + 1), alpha=0.7,
                    edgecolor="black", hatch=symbol2hatch[s_nb], align="left")
            pairs = pair_counts.get(frozenset([s_ct, s_nb]), 0)
            ax.text(0.95, 0.95, f"{s_ct}-{s_nb} pairs: {pairs}", transform=ax.transAxes,
                    ha="right", va="top", bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

            if i == 0:  # 第一行设置列标题
                ax.set_title(f"Neighbor: {s_nb}")
            if j == 0:  # 第一列设置行标题
                ax.set_ylabel(f"Center: {s_ct}")

            ax.set_xlim(0, 12)  # 最大配位数 12
            ax.grid(True, alpha=0.3)

    # 调整布局并保存
    fig.suptitle("Coordination Number Distribution", fontsize=16)
    plt.tight_layout()
    output = os.path.join(outdir, "cn-histogram-faceted.png")
    plt.savefig(output)
    plt.close()
    logging.info(f"配位数分面直方图已保存到 {output}")


def plot_histogram_stacked(cndata_list: list[CNData], pair_counts: dict[frozenset, int], outdir: str):
    # 统计配位数据
    cn_stats = defaultdict(lambda: defaultdict(list))
    for cndata in cndata_list:
        s_ct, s_nb = cndata.symbols
        cn_stats[s_ct][s_nb].append(cndata.cn)
    cn_stats = {k: v for k, v in sorted(cn_stats.items(), key=lambda x: x[0])}

    # 每个中心绘制堆叠图
    for s_ct, data_nb in cn_stats.items():
        data_nb = {k: v for k, v in sorted(data_nb.items(), key=lambda x: x[0])}
        s_nb_list = list(data_nb.keys())
        data = [data_nb[s_nb] for s_nb in s_nb_list]
        cn_max = max(max(d) for d in data) if data else 12
        pairs = [f"{s_ct}-{s_nb} pairs: {pair_counts.get(frozenset([s_ct, s_nb]), 0)}"
                 for s_nb in s_nb_list]

        hatches = [_hatch_patterns[i % len(_hatch_patterns)] for i in range(len(s_nb_list))]

        plt.figure(figsize=(10, 6))
        plt.hist(data, bins=range(0, cn_max + 1), label=pairs, alpha=0.7,
                 edgecolor="black", hatch=hatches, stacked=True, align="left")
        plt.title(f"Coordination Number of Center {s_ct}")
        plt.xlabel("Coordination Number")
        plt.ylabel("Frequency")
        plt.legend(title="Neighbor")
        plt.grid(True, alpha=0.3)
        plt.xticks(range(0, cn_max + 1))

        output = os.path.join(outdir, f"cn-histogram-{s_ct}.png")
        plt.savefig(output, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"配位数堆叠直方图已保存到 {output}")


def plot_heatmap(cndata_list: list[CNData], outdir: str):

    # 分组配位数据
    cndata_dict = defaultdict(list)
    for cndata in cndata_list:
        cndata_dict[cndata.symbols].append(cndata)
    cndata_dict = {k: v for k, v in sorted(cndata_dict.items(), key=lambda x: x[0])}

    cn_avg = defaultdict(dict)
    for (s_ct, s_nb), cndata_list in cndata_dict.items():
        cn_avg[s_ct][s_nb] = np.mean([d.cn for d in cndata_list])

    df = pd.DataFrame(cn_avg).fillna(0)  # center 为行, neighbor 为列

    plt.figure(figsize=(8, 8))
    sns.heatmap(df, annot=True, fmt=".1f", cmap="YlGnBu", cbar_kws={"label": "Average CN"})
    plt.title("Average Coordination Number Heatmap")
    plt.xlabel("Center Atom")
    plt.ylabel("Neighbor Atom")
    plt.tight_layout()
    output = os.path.join(outdir, "cn-heatmap.png")
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"配位热力图已保存到 {output}")


def countCN2files(filepath: str, outdir: str) -> str:
    """
    计算配位数并为每种配位类型生成POSCAR文件

    Args:
        filepath: 输入POSCAR文件的路径
        outdir: 生成文件的输出目录

    Returns:
        包含输出文件的目录路径
    """
    # 创建输出目录
    outdir = os.path.join(outdir, f"{os.path.splitext(os.path.basename(filepath))[0]}-cn")
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)

    # 读取POSCAR
    atoms = SimplePoscar.read_poscar(filepath)
    logging.debug(f"原子 {atoms}")

    # 检测截断距离
    cut_off = detect_cutoff_distance(atoms)
    logging.info(f"自动检测到的截断距离 {cut_off:.3f} Å")

    # 搜索最近邻
    cndata_list, pair_counts = calculate_nearest_neighbors(atoms=atoms, cut_off=cut_off)

    # 生成POSCAR
    cn_file_map = generate_poscar(atoms=atoms, cndata_list=cndata_list, outdir=outdir)
    merge_by_cn(cn_file_map=cn_file_map, outdir=outdir)

    # 将CN数据写入CSV
    save_dataframe(cndata_list=cndata_list, outdir=outdir)

    # 绘制直方图
    plot_histogram_faceted(cndata_list=cndata_list, pair_counts=pair_counts, outdir=outdir)
    plot_histogram_stacked(cndata_list=cndata_list, pair_counts=pair_counts, outdir=outdir)
    plot_heatmap(cndata_list=cndata_list, outdir=outdir)

    return outdir
