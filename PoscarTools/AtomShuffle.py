"""shuffle.py"""

import logging
import os
import random
from typing import Optional

from tqdm import tqdm

from .SimplePoscar import Atoms, SimplePoscar


def shuffle_atoms(atoms: Atoms, symbol_sites: dict[str, str], seed: Optional[int] = None):
    """Shuffle atoms within the same species

    Args:
        atoms (Atoms): Atoms to shuffle
        symbol_sites (dict[str, str]): Dictionary mapping symbol to sublattice site
        seed (int, optional): Random seed for reproducibility
    """
    if seed is not None:
        random.seed(seed)

    pbar = tqdm(total=len(atoms), ncols=80, desc="Shuffling atoms")

    new_atoms = atoms.copy(clean=True)
    for symbol, symbol_atoms in atoms.group_atoms:
        site = symbol_sites.get(symbol, "xx")

        # Shuffle atoms within the same site
        random.shuffle(symbol_atoms)

        # Add comments with site
        l = len(str(len(symbol_atoms)))
        for idx, atom in enumerate(symbol_atoms, start=1):
            atom.comment = f"{site}-{symbol}-#{atom.index+1:0{l}d} {idx:0{l}d}"
            pbar.update(1)

        new_atoms.extend(symbol_atoms)
    pbar.close()
    return new_atoms


def shuffle2files(filepath: str, structure: dict[str, dict],
                  seeds: list[int | None] = [None]) -> list[str]:
    """Shuffle atoms and save to series of files"""
    # Generate symbol_sites mapping
    info = structure.copy()
    info.pop("cell")
    symbol_sites = {}  # e.g. {'Ag': '1a', 'Cu': '3c'}
    for site, value in info.items():
        symbol, _ = value["atoms"]
        symbol_sites[symbol] = site

    # Read POSCAR
    poscar = SimplePoscar()
    atoms = poscar.read_poscar(filepath)

    # Shuffle atoms for each time
    outputs = []
    for t, seed in enumerate(seeds, start=1):
        sl = len(seeds)
        logging.info(f"\nWorking... {t}/{sl}")
        new_atoms = shuffle_atoms(atoms.copy(), symbol_sites, seed)

        # Save to file
        output = f"{os.path.splitext(filepath)[0]}-r{t:0{len(str(sl))}d}.vasp"
        poscar.write_poscar(output, new_atoms)
        logging.info(f"POSCAR Saved to {output}")
        outputs.append(output)

    return outputs
