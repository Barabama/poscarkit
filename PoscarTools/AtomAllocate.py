# AtomAllocate.py

import logging
import os
import random
import re
from collections import defaultdict

import numpy as np
from tqdm import tqdm

from .SimplePoscar import Atoms, SimplePoscar


def _integer_fracts(fracts: dict[str, float], factors: tuple[int, int, int], multi: int) -> dict:
    """将小数形式的分数转换为整数形式的分数

    Args:
        fracts: 元素和小数形式分数的字典
        factors: 超胞因子
        multi: 亚晶格的多重性

    Returns:
        dict: 元素和整数形式分数的字典
    """
    factor = np.prod(factors)
    # Super_fracts 是每个 fract * multi * factor
    super_fracts = {s: f * multi * factor for s, f in fracts.items()}
    # 四舍五入得到整数
    rounded_fracts = {s: round(f) for s, f in super_fracts.items()}
    total_rounded = sum(rounded_fracts.values())
    target_total = multi * factor

    # 调整四舍五入误差
    if total_rounded != target_total:
        # 计算差异并根据接近下一个整数的程度进行调整
        diffs = [(s, float(abs(f - round(f))), 1 if f - round(f) > 0 else -1)
                 for s, f in super_fracts.items()]
        diffs.sort(key=lambda x: x[1], reverse=True)
        # 确定需要的调整量
        adjustment = target_total - total_rounded

        # 应用调整
        for i in range(abs(adjustment)):
            symbol, decimal, direction = diffs[i]
            rounded_fracts[symbol] += direction

    return rounded_fracts


def allocate_atoms(atoms: Atoms, site_fracts: dict[str, dict[str, int]] | None = None,
                   seed: int | None = None) -> dict[str, Atoms]:
    """根据亚点阵的分数分配原子

    Args:
        atoms: 待分配的原子
        site_fracts: 字典 {site: {symbol: fractions}
        seed: 用于洗牌的随机种子, 以确保可重现性
    Returns:
        dict: 字典 {site: 已分配的原子}
    """
    if seed is not None:
        random.seed(seed)

    site_subatoms = {}
    pbar = tqdm(total=len(atoms), ncols=80, desc="分配原子")
    for note, subatoms in atoms.group_atoms(key="note"):
        match = re.search(r"(\d+[a-z])-([A-Za-z]+)", note)
        if not match:
            raise ValueError(f"无法识别注释 {note} 来确定亚点阵, 位于({subatoms})")
        site, symbol = match.groups()
        for i, atom in enumerate(subatoms):
            atom.index = i

        # 在同一位点内打乱原子顺序
        sub_list = subatoms.atom_list
        random.shuffle(sub_list)

        if site_fracts is None:
            symbols = [atom.symbol for atom in sub_list]
        elif site in site_fracts:
            symbols = [s for s, f in site_fracts[site].items() for i in range(f)]
        else:
            raise ValueError(f"占位分数 {site_fracts} 中未找到亚点阵 {site}")

        # 分配元素和元数据
        slsl = len(str(len(sub_list)))
        for idx, (symbol, atom) in enumerate(zip(symbols, sub_list), start=1):
            atom.symbol = symbol
            atom.meta = f"{idx:0{slsl}d}"
            pbar.update(1)

        subatoms = subatoms.copy(atom_list=sub_list)
        site_subatoms.update({site: subatoms})
    pbar.close()
    return site_subatoms


def allocate2files(filepath: str, outdir: str, factors: tuple[int, int, int],
                   struct_info: dict[str, dict] | None = None,
                   seeds: list[int | None] = [None]) -> list[str]:
    """根据亚点阵的分数分配原子

    Args:
        filepath: POSCAR文件路径
        outdir: 输出目录
        factors: 超胞因子
        struct_info: 结构信息
        seeds: 用于洗牌的种子, 默认为[None]

    Returns:
        list: 输出文件路径列表
    """
    # 读取POSCAR
    atoms = SimplePoscar.read_poscar(filepath)
    logging.debug(f"原子: {atoms}")

    # 生成整数形式的位点分数
    if struct_info is None:
        site_fracts = None
    else:
        struct_info = struct_info.copy()
        struct_info.pop("cell", None)
        site_fracts = defaultdict(dict)
        for site, data in struct_info.items():
            if "sofs" not in data:
                raise ValueError(f"亚点阵 {site} 没有sofs数据")
            fracts: dict = data["sofs"]
            if abs(sum(fracts.values()) - 1) > 1e-6:
                raise ValueError(f"亚点阵 {site} 的分数之和不接近于1, 请检查配置")
            site_fracts[site] = _integer_fracts(fracts=fracts, factors=factors, multi=int(site[0]))
        logging.info(f"占位分数 {dict(site_fracts)}")

    outputs = []
    sl = len(seeds)
    ssl = len(str(sl))
    for t, seed in enumerate(seeds, start=1):
        logging.info(f"分配 {t}/{sl}")
        new_atoms = atoms.copy(clean=True)
        site_subatoms = allocate_atoms(atoms=atoms, site_fracts=site_fracts, seed=seed)
        for site, subatoms in site_subatoms.items():
            new_atoms.extend(subatoms)
            symbol_str = "".join(s for s, c in subatoms.symbol_count)
            output = os.path.join(outdir, f"POSCAR-allocate{t:0{ssl}d}-{site}-{symbol_str}.vasp")
            comment = f"Allocated-seed={seed}-{site}-{symbol_str}"
            SimplePoscar.write_poscar(filepath=output, atoms=subatoms, comment=comment)

        # 保存到文件
        logging.debug(f"已分配原子 {new_atoms}")
        symbol_str = "".join(s for s, c in new_atoms.symbol_count)
        output = os.path.join(outdir, f"POSCAR-allocate{t:0{ssl}d}-{symbol_str}.vasp")
        comment = f"Allocated-seed={seed}-{symbol_str}"
        SimplePoscar.write_poscar(filepath=output, atoms=new_atoms, comment=comment)
        outputs.append(output)
        logging.info(f"POSCAR已保存到 {output}")

    return outputs