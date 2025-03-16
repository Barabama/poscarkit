"""AtomSeparate.py"""

import logging
import os
import re
from collections import defaultdict
from typing import TypedDict

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import KDTree
from tqdm import tqdm

from .SimplePoscar import Atom, Atoms, SimplePoscar
from .Utils import color_map


class CNData(TypedDict):
    center: Atom
    neighbors: list[Atom]
    cn: int


def group_by_sites(atoms: Atoms):
    """Group atoms by their site names."""
    grouped_atoms = defaultdict(list)
    for atom in atoms:
        match = re.search(r"([^#-]+)", atom.comment)
        if match:
            site = match.group(1)
            grouped_atoms[site].append(atom)

    return grouped_atoms


def separate2files(filepath: str) -> list[str]:
    """Separate atoms by their site names and save to separate files."""
    poscar = SimplePoscar()
    atoms = poscar.read_poscar(filepath)

    logging.info(f"Separating {filepath}...")
    outputs = []
    for site, atom_list in group_by_sites(atoms).items():
        new_atoms = atoms.rebuild(atom_list)
        output = f"{os.path.splitext(filepath)[0]}-{site}.vasp"
        poscar.write_poscar(output, new_atoms)
        logging.info(f"POSCAR saved to {output}")
        outputs.append(output)

    return outputs


def calculate_nearest_neighbors(atoms: Atoms, cut_off: float):
    """
    Calculate the nearest neighbors of each atom and their coordination numbers
    within the given cut-off distance.
    Use KDTree to improve efficiency and reduce memory usage.
    """
    # coords = atoms.direct_coords
    # if res := verify_atoms(coords=coords, atoms=atoms):
    #     raise ValueError(res)

    # Get cartesian coordinates
    coords = atoms.cartesian_coords

    # Make KDTree
    tree = KDTree(coords)

    nn_map = defaultdict(list)
    pair_count = 0

    for i, coord in enumerate(tqdm(coords, desc="Searching for NN", ncols=80)):
        # Search for neighbors in rqstd + tolerance
        neighbor_indices = tree.query_ball_point(coord, r=1.5 * cut_off)

        # Calculate distances
        distances = np.linalg.norm(coords[neighbor_indices] - coord, axis=1)

        # Filter out neighbors
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


def countCN2files(filepath: str):
    cut_off = float(input("Please input the cut-off distance (A): "))

    poscar = SimplePoscar()
    atoms = poscar.read_poscar(filepath)

    output = os.path.splitext(filepath)[0]
    if not os.path.exists(output):
        os.makedirs(output)

    symbol_df = defaultdict(pd.DataFrame)
    symbol_cn_freq = defaultdict(list)
    for symbol, atom_list in atoms.group_atoms:
        symbol_atoms = Atoms(atoms.cell, is_direct=atoms.is_direct, atom_list=atom_list)
        del atom_list
        # Get nearest neighbors
        nn_map, pair_count = calculate_nearest_neighbors(symbol_atoms, cut_off)
        logging.info(f"{symbol}-{symbol} pair count: {pair_count}")

        # Collect CN data map
        cn_data_map = defaultdict(list[CNData])
        for center, neighbors in nn_map.items():
            cn = len(neighbors)
            cn_data_map[cn].append(CNData(center=center, neighbors=neighbors, cn=cn))

        # Save nearest neighbor atoms to POSCAR file
        nn_atoms = set()
        for cn, cn_data_list in cn_data_map.items():
            atom_sets = set()
            for cn_data in cn_data_list:
                sets = {cn_data["center"], *cn_data["neighbors"]}
                atom_sets.update(sets)

            filename = os.path.join(output, f"POSCAR-nn-{symbol}-{cn}.vasp")
            comment = f"Coordination number of {symbol}={cn}"
            poscar.write_poscar(filename, atoms.rebuild(list(atom_sets)), comment)

            nn_atoms.update(atom_sets)

        comment = f"Nearest Neighbors of {symbol} pair count={pair_count}"
        filename = os.path.join(output, f"POSCAR-nn-{symbol}.vasp")
        poscar.write_poscar(filename, atoms.rebuild(list(nn_atoms)), comment)

        # Collect CN data to df
        cn_counts = [[f"{symbol}*-{cn}{symbol}", len(d)] for cn, d in cn_data_map.items()]
        cn_df = pd.DataFrame(data=cn_counts, columns=["CN", "Count"])
        symbol_df[symbol] = cn_df

        # Collect CN data to list
        cn_freq = [d["cn"] for ds in cn_data_map.values() for d in ds]
        symbol_cn_freq[symbol] = cn_freq

    # Write CN data to CSV
    all_df = pd.concat(symbol_df.values(), ignore_index=True)
    filename = os.path.join(output, "CN_Counts.csv")
    all_df.to_csv(filename, index=False)
    logging.info(f"Coordination Number counts saved to {output}")

    # Plot the histogram of Coordinate Numbers
    symbols, cn_freqs = zip(*symbol_cn_freq.items())
    colors = [color_map[s] for s in symbols]
    plt.figure(figsize=(10, 6))
    plt.hist(cn_freqs, bins=list(range(0, 12)),
             alpha=0.5, label=symbols, color=colors, edgecolor="black")
    plt.legend()
    plt.title(f"Histogram of Coordination Numbers")
    plt.xlabel("Coordinate Number")
    plt.ylabel("Frequency")
    
    # Define ticks
    plt.xticks(range(0, 12))
    
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.savefig(os.path.join(output, f"CN-Histogram.png"))
    plt.close()

    return output
