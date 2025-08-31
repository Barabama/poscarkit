# AtomSlice.py

import logging
import os
import shutil
from collections import defaultdict
from itertools import groupby

import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm

from .SimplePoscar import Atoms, SimplePoscar
from .Utils import color_map

basis_map = {
    (0, 0, 1): [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    (1, 1, 0): [(0, 0, -1), (-1, 1, 0), (1, 1, 0)],
    (1, 1, 1): [(1, 1, -2), (-1, 1, 0), (1, 1, 1)],
}


def _normalize(vector: np.ndarray) -> np.ndarray:
    return vector / np.linalg.norm(vector)


def _get_basis(miller_index: tuple[int, int, int]) -> list[np.ndarray]:
    """通过晶面指数找到基向量

    Args:
        miller_index: 平面的晶面指数
    Returns:
        list: 3个基向量
    """
    if miller_index in basis_map:
        basis = [np.array(v) for v in basis_map[miller_index]]
    else:
        n = np.array(miller_index)
        # 找到与晶面指数垂直的两个基向量
        t0 = np.array([1, 0, 0]) if abs(n[0]) < abs(n[1]) else np.array([0, 1, 0])
        b1 = np.cross(n, t0)
        b2 = np.cross(n, b1)
        basis = [b1, b2, n]
    return basis


def _convert(atoms: Atoms, basis: np.ndarray) -> Atoms:
    # from .AtomSupercell import make_supercell
    from ase.build.tools import cut
    ase_atoms = SimplePoscar.to_ase_atoms(atoms)
    a, b, c = basis
    converted = cut(ase_atoms, a, b, c)
    new_atoms = SimplePoscar.from_ase_atoms(converted)
    # new_atoms = make_supercell(atoms=new_atoms, factors=(1, 1, 1))  # 标准化
    return new_atoms


def group_by_normal(atoms: Atoms, basis: np.ndarray, precision: int = 2):
    """按基向量法线方向的投影距离对原子进行分组

    Args:
        atoms: Atoms对象
        basis: 基向量
        precision: 四舍五入的小数位数默认为6
    Yields:
        tuple: 投影, 层
    """
    # 计算并四舍五入投影
    coords = atoms.cartesian_coords
    projs = np.dot(coords, basis[2])  # 在法线上的投影
    projs = np.round(projs, precision)

    # 根据四舍五入的投影对原子进行排序
    sorted_indices = np.argsort(projs)
    for proj, group in groupby(sorted_indices, key=lambda x: projs[x]):
        layer = atoms.copy(atom_list=[atoms[i] for i in group])
        yield proj, layer


def plot_layer(layer: Atoms, basis: list[np.ndarray], title: str, filepath: str,
               pair_counts: dict[str, int] | None = None):
    """按基向量绘制层

    Args:
        layer: 层的Atoms对象
        basis: 基向量
        title: 图表标题
        filepath: 保存图表的文件路径
    """
    # 计算在法线上的投影以获得投影坐标
    layer = layer.sort()
    coords = layer.cartesian_coords
    b1, b2, n = [_normalize(v) for v in layer.cell]
    n_projs = np.dot(coords, n)  # 在法线上的投影
    p_projs = coords - np.outer(n_projs, n)  # 在平面上的投影
    # xs = np.dot(p_projs, b1)  # 在b1上的分量
    # ys = np.dot(p_projs, b2)  # 在b2上的分量
    proj_coords = np.column_stack((np.dot(p_projs, b1), np.dot(p_projs, b2)))

    # 按元素对投影坐标进行分组
    symbol_coords = defaultdict(list)
    for atom, coord in zip(layer, proj_coords):
        symbol_coords[atom.symbol].append(coord)

    # 获取基向量的范围
    x_min, x_max = 0.0, np.linalg.norm(layer.cell[0])
    y_min, y_max = 0.0, np.linalg.norm(layer.cell[1])
    x_margin = (x_max - x_min) * 0.1
    y_margin = (y_max - y_min) * 0.1

    # 使用投影坐标绘制层
    plt.figure(figsize=(6, 6))
    for symbol, coords in symbol_coords.items():
        color = color_map.get(symbol, "#FF00FF")
        x, y = zip(*coords)
        # 准备图例标签，如果提供了pair_counts则包含在内
        label = symbol if not (pair_counts and symbol in pair_counts) else \
            f"{symbol}-{symbol} pairs: {pair_counts[symbol]}"
        plt.scatter(x, y, marker="o", s=10, color=color, alpha=1.0, label=label)

    plt.title(title)
    plt.xlabel(f"[{' '.join(str(v) for v in basis[0])}] Coordinate (Å)")
    plt.ylabel(f"[{' '.join(str(v) for v in basis[1])}] Coordinate (Å)")
    # plt.axis("equal")
    plt.grid()
    plt.legend(title="Symbols", bbox_to_anchor=(1, 1), loc="upper left")
    plt.xlim(-x_margin, x_max + x_margin)
    plt.ylim(-y_margin, y_max + y_margin)
    # plt.tight_layout(rect=[0, 0, 1, 0])
    plt.savefig(filepath, bbox_inches="tight")
    plt.close()


def slice2file(filepath: str, outdir: str, miller_index: tuple[int, int, int]) -> str:
    """通过晶面指数切片POSCAR"""
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

        # 按基向量绘制层
        imgname = os.path.join(outdir, f"{comment}.png")
        plot_layer(layer=layer, basis=basis, title=comment, filepath=imgname)
        # break  # 用于测试

    logging.info(f"结果已保存在 {outdir}")
    return outdir
