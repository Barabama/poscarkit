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
    """Convert decimal fractions to integer fractions .

    Args:
        fracts (dict): Dictionary of symbol and decimal fractions.
        factors (tuple[int, int, int]): Factors for supercell.
        multi (int): Multiplicity of sublattice.

    Returns:
        dict: Dictionary of symbol and integer fractions.
    """
    factor = np.prod(factors)
    # Super_fracts is every fract * multi * factor
    super_fracts = {s: f * multi * factor for s, f in fracts.items()}
    # Round the values to get integers
    rounded_fracts = {s: round(f) for s, f in super_fracts.items()}
    total_rounded = sum(rounded_fracts.values())
    target_total = multi * factor

    # Adjusting rounding errors
    if total_rounded != target_total:
        # Calculate differences and adjust based on closeness to the next integer
        diffs = [(s, float(abs(f - round(f))), 1 if f - round(f) > 0 else -1)
                 for s, f in super_fracts.items()]
        diffs.sort(key=lambda x: x[1], reverse=True)
        # Determine the adjustment needed
        adjustment = target_total - total_rounded

        # Apply adjustments
        for i in range(abs(adjustment)):
            symbol, decimal, direction = diffs[i]
            rounded_fracts[symbol] += direction

    return rounded_fracts


def allocate_atoms(atoms: Atoms, site_fracts: dict[str, dict[str, int]] | None = None,
                   seed: int | None = None) -> dict[str, Atoms]:
    """Allocate atoms according to site of fractions.

    Args:
        atoms (Atoms): Atoms to allocated.
        site_fracts (dict[str, dict[str, int]]): Dict {site: {symbol: fractions}.
        seed (int | None): Random seed to shuffle for reproducibility.
    Returns:
        dict[str, Atoms]: Dict {site: Allocated subatoms}.
    """
    if seed is not None:
        random.seed(seed)

    site_subatoms = {}
    pbar = tqdm(total=len(atoms), ncols=80, desc="Allocating atoms")
    for note, subatoms in atoms.group_atoms(key="note"):
        match = re.search(r"(\d+[a-z])-([A-Za-z]+)", note)
        if not match:
            raise ValueError(f"Unknown note({note}) to recognize the site.")
        site, symbol = match.groups()
        for i, atom in enumerate(subatoms):
            atom.index = i

        # Shuffle atoms within the same site
        sub_list = subatoms.atom_list
        random.shuffle(sub_list)

        if site_fracts is None:
            symbols = [atom.symbol for atom in sub_list]
        elif site in site_fracts:
            symbols = [s for s, f in site_fracts[site].items() for i in range(f)]
        else:
            raise ValueError(f"Site({site}) not found in site_fracts({site_fracts}).")

        # Assign symbols and meta
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
    """Allocate atoms according to the site of fractions.

    Args:
        filepath (str): POSCAR file path.
        outdir (str): Output directory.
        factors (tuple[int, int, int]): Supercell factors.
        struct_info (dict[str, dict], optional): Structure information.
        seeds (list[int | None], optional): Seeds for shuffling. Defaults to [None].

    Returns:
        list[str]: Output file paths.
    """
    # Read POSCAR
    atoms = SimplePoscar.read_poscar(filepath)
    logging.debug(f"Atoms: {atoms}")

    # Generate integer site of fractions
    if struct_info is None:
        site_fracts = None
    else:
        struct_info = struct_info.copy()
        struct_info.pop("cell", None)
        site_fracts = defaultdict(dict)
        for site, data in struct_info.items():
            if "sofs" not in data:
                raise ValueError(f"SOFs data not found for site({site})")
            fracts: dict = data["sofs"]
            if abs(sum(fracts.values()) - 1) > 1e-6:
                raise ValueError(f"The sum of fractions for site({site}) not close to 1. Check config.")
            site_fracts[site] = _integer_fracts(fracts=fracts, factors=factors, multi=int(site[0]))
        logging.info(f"Site of fractions: {dict(site_fracts)}")

    outputs = []
    sl = len(seeds)
    ssl = len(str(sl))
    for t, seed in enumerate(seeds, start=1):
        logging.info(f"Allocating {t}/{sl}")
        new_atoms = atoms.copy(clean=True)
        site_subatoms = allocate_atoms(atoms=atoms, site_fracts=site_fracts, seed=seed)
        for site, subatoms in site_subatoms.items():
            new_atoms.extend(subatoms)
            symbol_str = "".join(s for s, c in subatoms.symbol_count)
            output = os.path.join(outdir, f"POSCAR-allocate{t:0{ssl}d}-{site}-{symbol_str}.vasp")
            comment = f"Allocated-seed={seed}-{site}-{symbol_str}"
            SimplePoscar.write_poscar(filepath=output, atoms=subatoms, comment=comment)

        # Save to file
        logging.debug(f"Allocated: {new_atoms}")
        symbol_str = "".join(s for s, c in new_atoms.symbol_count)
        output = os.path.join(outdir, f"POSCAR-allocate{t:0{ssl}d}-{symbol_str}.vasp")
        comment = f"Allocate-seed={seed}-{symbol_str}"
        SimplePoscar.write_poscar(filepath=output, atoms=new_atoms, comment=comment)
        outputs.append(output)
        logging.info(f"POSCAR Saved to {output}")

    return outputs
