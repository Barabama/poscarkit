# workflows/modeling.py

import logging
import os

from PoscarTools.AtomSupercell import unitcell2file, supercell2file
from PoscarTools.AtomAllocate import allocate2files


def run_modeling(filepath: str, outdir: str, supercell_factors: tuple[int, int, int],
                 struct_info: dict[str, dict] | None = None,
                 shuffle_seeds: list[int | None] = [None]) -> list[str]:
    """
    运行建模工作流：首先生成超胞，然后分配原子

    Args:
        filepath: 输入POSCAR文件路径
        outdir: 输出目录
        supercell_factors: 超胞因子 (x, y, z)
        struct_info: 结构信息（用于生成单胞）
        shuffle_seeds: 随机种子列表，用于原子分配

    Returns:
        list[str]: 分配后的POSCAR文件路径列表
    """
    # 处理超胞生成
    if not filepath:
        logging.warning("No filepath provided. Using structure info instead.")
        if struct_info is None:
            raise ValueError("Either filepath or struct_info must be provided")
        filepath = unitcell2file(struct_info=struct_info, outdir=outdir)
    else:
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"No such file: {filepath}")

    # 生成超胞
    supercell_filepath = supercell2file(filepath=filepath, outdir=outdir, factors=supercell_factors)
    logging.info(f"Supercell generated: {supercell_filepath}")

    # 分配原子
    allocated_filepaths = allocate2files(filepath=supercell_filepath, outdir=outdir,
                                         factors=supercell_factors, struct_info=struct_info,
                                         seeds=shuffle_seeds)
    logging.info(f"Atom allocation completed. Generated files: {allocated_filepaths}")

    return allocated_filepaths
