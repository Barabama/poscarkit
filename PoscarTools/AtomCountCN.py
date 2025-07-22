# AtomCountCN.py

import logging
import os
import shutil
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


def calculate_nearest_neighbors(atoms: Atoms, cut_off: float):
    """
    Calculate the nearest neighbors of each atom and their coordination numbers
    within the given cut-off distance.
    Use KDTree to improve efficiency and reduce memory usage.
    """
    nn_map = defaultdict(list)
    pair_count = 0

    # Get cartesian coordinates
    coords = atoms.cartesian_coords
    # Make KDTree
    tree = KDTree(coords)
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


def countCN2files(filepath: str, outdir: str):
    cut_off = float(input("Please input the cut-off distance (A): "))

    # Make output directory
    outdir = os.path.join(outdir, f"{os.path.splitext(os.path.basename(filepath))[0]}-cn")
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)

    # Read POSCAR
    atoms = SimplePoscar.read_poscar(filepath)
    logging.debug(f"Atoms: {atoms}")

    # Count nearest neighbors and coordinate number
    symbol_df = defaultdict(pd.DataFrame)
    symbol_cn_freq = defaultdict(list)
    for symbol, subatoms in atoms.group_atoms():
        # Get nearest neighbors
        nn_map, pair_count = calculate_nearest_neighbors(subatoms, cut_off)
        logging.info(f"{symbol}-{symbol} pair count: {pair_count}")

        # Collect CN data map
        cn_data_map = defaultdict(list[CNData])
        for center, neighbors in nn_map.items():
            cn = len(neighbors)
            cn_data_map[cn].append(CNData(center=center, neighbors=neighbors, cn=cn))

        # Get each coordinate number atoms
        nn_atoms_set = set()
        for cn, cn_data_list in cn_data_map.items():
            cn_atom_sets = set()
            for cn_data in cn_data_list:
                sets = {cn_data["center"], *cn_data["neighbors"]}
                cn_atom_sets.update(sets)

            nn_atoms_set.update(cn_atom_sets)
            cn_atoms = atoms.copy(atom_list=list(cn_atom_sets))
            logging.debug(f"{symbol} {cn} atoms: {cn_atoms}")

            # Save to file
            output = os.path.join(outdir, f"POSCAR-d1nn-{symbol}*-{cn}{symbol}.vasp")
            comment = f"CoordinationNumber-{symbol}={cn}"
            SimplePoscar.write_poscar(filepath=output, atoms=cn_atoms, comment=comment)

        # Get nearest neighbors atoms
        nn_atoms = atoms.copy(atom_list=list(nn_atoms_set))
        logging.debug(f"NN atoms: {nn_atoms}")

        # Save to file
        output = os.path.join(outdir, f"POSCAR-d1nn-{symbol}.vasp")
        comment = f"NearestNeighbors-{symbol}-pair_count={pair_count}"
        SimplePoscar.write_poscar(filepath=output, atoms=nn_atoms, comment=comment)

        # Collect CN data to df
        cn_counts = [[f"{symbol}*-{cn}{symbol}", len(d)] for cn, d in cn_data_map.items()]
        cn_df = pd.DataFrame(data=cn_counts, columns=["CN", "Count"])
        symbol_df[symbol] = cn_df

        # Collect CN data to list
        cn_freq = [d["cn"] for ds in cn_data_map.values() for d in ds]
        symbol_cn_freq[symbol] = cn_freq

    # Write CN data to CSV
    all_df = pd.concat(symbol_df.values(), ignore_index=True)
    output = os.path.join(outdir, "CN_Counts.csv")
    all_df.to_csv(output, index=False)
    logging.info(f"Coordination Number counts saved to {output}")

    # Plot the histogram of Coordinate Numbers
    symbols, cn_freqs = zip(*symbol_cn_freq.items())
    colors = [color_map[s] for s in symbols]
    plt.figure(figsize=(10, 6))
    plt.hist(cn_freqs, bins=list(range(0, 12)), alpha=0.5, label=symbols, color=colors, edgecolor="black")
    plt.legend()
    plt.title("Histogram of Coordination Numbers")
    plt.xlabel("Coordinate Number")
    plt.ylabel("Frequency")
    plt.xticks(range(0, 12))  # Define ticks
    plt.grid(True, linestyle='--', alpha=0.7)
    output = os.path.join(outdir, "CN-Histogram.png")
    plt.savefig(output)
    plt.close()
    logging.info(f"Coordination Number histogram saved to {output}")

    return output
