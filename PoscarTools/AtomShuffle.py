# AtomShuffle.py

import logging
import os
import random

from tqdm import tqdm

from .SimplePoscar import Atoms, read_poscar, write_poscar


def shuffle_atoms(atoms: Atoms, symbol_sites: dict[str, str], seed: int | None = None):
    """Shuffle atoms within the same species

    Args:
        atoms (Atoms): Atoms to shuffle
        symbol_sites (dict[str, str]): Dictionary mapping symbol to sublattice site
        seed (int | None): Random seed for reproducibility
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
    info.pop("cell", None)
    symbol_sites = {}  # e.g. {'Au': '1a', 'Cu': '3c'}
    for site, value in info.items():
        symbol, coords = value["atoms"]
        symbol_sites[symbol] = site
    logging.debug(f"Symbol sites: {symbol_sites}")

    # Read POSCAR
    atoms = read_poscar(filepath)
    logging.debug(f"Atoms: {atoms}")

    # Shuffle atoms for each time
    output = os.path.splitext(filepath)[0]
    if not os.path.exists(output):
        os.makedirs(output)

    # Check if sites are defined
    for s, c in atoms.symbol_count:
        if s not in symbol_sites:
            logging.warning(f"Unknown site for {s}")

    outputs = []
    for t, seed in enumerate(seeds, start=1):
        sl = len(seeds)
        logging.info(f"Shuffling {t}/{sl}")
        new_atoms = shuffle_atoms(atoms.copy(), symbol_sites, seed)
        logging.debug(f"Shuffled: {new_atoms}")

        # Save to file
        filename = os.path.join(output, f"POSCAR-r{t:0{len(str(sl))}d}.vasp")
        write_poscar(filename, new_atoms)
        outputs.append(filename)

    logging.info(f"POSCAR Saved to {output}")

    return outputs
