"""allocate.py"""

import os
import random
import tomllib
from collections import defaultdict

import numpy as np
from tqdm import tqdm

from SimplePoscar import Atoms, SimplePoscar


def _integer_fractions(fracts: dict, factors: tuple[int, int, int], multi: int) -> dict:
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
        diffs = {s: (abs(f - round(f)), 1 if f - round(f) > 0 else -1) for s, f in super_fracts.items()}
        sorted_diffs = sorted(diffs.items(), key=lambda x: x[1][0], reverse=True)

        # Determine the adjustment needed
        adjustment = target_total - total_rounded

        # Apply adjustments
        for i in range(abs(adjustment)):
            s, (_, direction) = sorted_diffs[i]
            rounded_fracts[s] += direction

    return rounded_fracts

    # ###### Hamilton's Method useless ######
    # # Step 1: Assign integer parts
    # int_fracts = {s: int(f) for s, f in super_fracts.items()}
    # remaining = total - sum(int_fracts.values())

    # # Step 2: Sort symbols by their decimal parts (largest first)
    # sorted_symbols = sorted(super_fracts.keys(),
    #                         key=lambda s: super_fracts[s] - int(super_fracts[s]),
    #                         reverse=True)

    # # Step 3: Assign remaining counts based on decimal parts
    # for symbol in sorted_symbols[:remaining]:
    #     int_fracts[symbol] += 1

    # return int_fracts


def allocate_atoms(atoms: Atoms, vac_sites: dict[str, str],
                   site_fracts: dict[str, dict[str, float]],
                   shuffle: bool = False) -> Atoms:
    """Allocate atoms according to the integer site fractions.

    Args:
        atoms (Atoms): Atoms template.
        vac_sites (dict[str, str]): Dictionary mapping vacancy symbol to sublattice site.
        site_fracts (dict[str, dict[str, float]]): Dictionary mapping symbol to site fractions.
        shuffle (bool, optional): Whether to shuffle atoms. Defaults to False.
    Returns:
        Atoms: Allocated atoms.
    """

    pbar = tqdm(total=len(atoms), ncols=80, desc="Allocating atoms")

    new_atoms = atoms.copy(clean=True)
    for vac, vac_atoms in atoms.group_atoms:
        site = vac_sites[vac]

        # Optional: Shuffle atoms
        if shuffle:
            random.shuffle(vac_atoms)

        # Assign symbols and comments
        l = len(str(len(vac_atoms)))
        symbol_iter = (s for s, c in site_fracts[site].items() for _ in range(c))
        for idx, (symbol, atom) in enumerate(zip(symbol_iter, vac_atoms), start=1):
            atom.symbol = symbol
            atom.comment = f"{site}-{vac}-#{atom.index+1:0{l}d} {idx:0{l}d} {symbol}"
            pbar.update(1)

        new_atoms.extend(vac_atoms)
    pbar.close()
    return new_atoms


def allocate2file(filepath: str, structure: dict[str, dict],
                  factors: tuple[int, int, int], shuffle: bool = False) -> str:
    """Allocate atoms according to the site fractions."""
    # Read POSCAR
    poscar = SimplePoscar()
    atoms = poscar.read_poscar(filepath)

    # Generate vacancy sites and site fractions
    info = structure.copy()
    info.pop("cell")
    vac_sites = {}  # e.g. {'Ag': '1a', 'Cu': '3c'}
    site_fracts = defaultdict(dict)  # template symbol to sofs
    for site, value in info.items():
        vac, _ = value["atoms"]
        vac_sites[vac] = site

        fracts: dict = value["sofs"]
        if abs(sum(fracts.values()) - 1) > 1e-6:
            raise ValueError("The sum of fractions must be 1.")
        site_fracts[site] = _integer_fractions(fracts, factors, int(site[0]))

    print(f"Site fractions: {dict(site_fracts)}")

    # Allocate atoms
    if shuffle:
        print("Would shuffle before allocating atoms.")

    new_atoms = allocate_atoms(atoms.copy(), vac_sites, site_fracts, shuffle)

    # Save to file
    symbol_str = "".join(s for s, _ in new_atoms.symbol_count)
    output = f"{os.path.splitext(filepath)[0]}-{symbol_str}.vasp"
    poscar.write_poscar(output, new_atoms)
    print(f"Allocated atoms saved to {output}")

    return output
