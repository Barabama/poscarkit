# AtomSupercell.py

import logging
import os

import numpy as np
from tqdm import tqdm

from .SimplePoscar import Atom, Atoms, SimplePoscar


def _clean_matrix(matrix: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """通过将小值设置为零来清理矩阵

    Args:
        matrix: 输入矩阵
        eps: 小于阈值`eps`的值将被设置为零

    Returns:
        np.ndarray: 清理后的矩阵
    """
    matrix = np.array(matrix)  # 确保输入是NumPy数组
    matrix[np.abs(matrix) < eps] = 0  # 将小值设置为零
    return matrix


def make_supercell(atoms: Atoms, factors: tuple[int, int, int]) -> Atoms:
    """制作超胞

    Args:
        atoms: Atoms对象
        factors: 超胞因子

    Returns:
        Atoms: 超胞原子
    """
    atoms = atoms.copy()
    # 原始坐标
    atoms.switch_coords(direct=True)  # 确保在直接坐标中
    coords = atoms.direct_coords

    # 超胞索引
    n, m, p = factors
    i, j, k = np.mgrid[0:n, 0:m, 0:p]  # 形状(n, m, p)
    indices = np.stack([i, j, k], axis=-1)  # 形状(n, m, p, 3)
    indices = indices.reshape(-1, 3)  # 形状(n*m*p, 3)

    # 广播
    super_coords = (coords[:, np.newaxis, :] + indices[np.newaxis, :, :]) / [n, m, p]
    super_coords = super_coords.reshape(-1, 3)  # 形状(N * n * m * p, 3)
    super_coords = super_coords % 1.0

    # 新原子
    matrix = np.array([[n, 0, 0], [0, m, 0], [0, 0, p]])
    new_cell = _clean_matrix(np.dot(atoms.cell, matrix))
    new_atoms = Atoms(cell=new_cell, is_direct=True)
    for idx, coord in enumerate(tqdm(super_coords, ncols=80, desc="生成超胞")):
        atom = atoms[idx // (n * m * p)]
        note = atom.note if atom.note else atom.symbol
        new_atoms.append(Atom(index=idx, symbol=atom.symbol, coord=coord,
                              constr=atom.constr, note=note, meta=atom.meta))

    return new_atoms


def unitcell2file(struct_info: dict[str, dict], outdir: str) -> str:
    """单胞结构信息保存到POSCAR文件"""
    # 获取晶胞信息
    if "cell" not in struct_info:
        raise ValueError(f"cell不在结构信息 {struct_info}")
    cell = np.array(struct_info.pop("cell"))
    if cell.size == 3:
        diagonal_matrix = np.zeros((3, 3))
        np.fill_diagonal(diagonal_matrix, cell)
        cell = diagonal_matrix
    if cell.shape != (3, 3):
        raise ValueError(f"cell {cell} 不是3x3矩阵")

    # 获取原子信息
    atom_list = []
    for idx, (site, data) in enumerate(struct_info.items()):
        if "atoms" not in data:
            raise ValueError(f"atoms不在结构信息, 位于 {site}")
        symbol, coords = data["atoms"]
        coords = np.array(coords)
        note = f"{site}-{symbol}"
        atom_list.extend([Atom(index=idx, symbol=symbol, coord=coord, note=note)
                          for coord in coords])

    atoms = Atoms(cell=cell, is_direct=True, atom_list=atom_list)
    output = os.path.join(outdir, "POSCAR-unitcell.vasp")
    comment = "Unitcell"
    SimplePoscar.write_poscar(filepath=output, atoms=atoms, comment=comment)

    logging.info(f"单胞已保存到 {output}")
    return output


def supercell2file(filepath: str, outdir: str, factors: tuple[int, int, int]) -> str:
    """制作超胞并保存到POSCAR文件"""
    # 读取原始POSCAR
    atoms = SimplePoscar.read_poscar(filepath)
    logging.debug(f"原子 {atoms}")

    # 制作超胞
    new_atoms = make_supercell(atoms, factors)
    logging.debug(f"超胞 {new_atoms}")

    # ### 使用ASE.make_supercell()制作超胞
    # from ase.build import make_supercell as ASEmake_supercell
    # ase_atoms = SimplePoscar.to_ase_atoms(atoms)
    # a, b, c = factors
    # matrix = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])
    # supercell = ASEmake_supercell(ase_atoms, matrix, order="atom-major")
    # new_atoms = SimplePoscar.from_ase_atoms(supercell)

    # 保存到文件
    factors_str = "x".join(str(f) for f in factors)
    output = os.path.join(outdir, f"POSCAR-supercell-{factors_str}.vasp")
    comment = f"Supercell-{factors_str}"
    SimplePoscar.write_poscar(filepath=output, atoms=new_atoms, comment=comment)

    logging.info(f"超胞已保存到 {output}")
    return output