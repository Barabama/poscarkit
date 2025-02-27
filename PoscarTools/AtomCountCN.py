"""AtomSeparate.py"""

import math
import os
import re
import sys
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree
from tqdm import tqdm

from SimplePoscar import Atoms, SimplePoscar
from AtomSlice import _get_basis, group_by_direction
from Utils import color_map


def group_by_sites(atoms: Atoms):
    grouped_atoms = defaultdict(list)
    for atom in atoms:
        match = re.search(r"([^#-]+)", atom.comment)
        if match:
            site = match.group(1)
            grouped_atoms[site].append(atom)

    return grouped_atoms


def separate2files(filepath: str):

    poscar = SimplePoscar()
    atoms = poscar.read_poscar(filepath)

    print(f"Separating {filepath}...")
    for site, atom_list in group_by_sites(atoms).items():
        new_atoms = atoms.copy(clean=True)
        new_atoms.extend(atom_list)
        output = f"{os.path.splitext(filepath)[0]}-{site}.vasp"
        poscar.write_poscar(output, new_atoms)
        print(f"POSCAR saved to {output}")
        yield output


def verify_atoms(coords: np.ndarray, atoms: Atoms):
    for i, coord in enumerate(tqdm(coords, desc="Verifying atoms", ncols=80)):
        atom = atoms[i]
        if not np.allclose(coord, atom.coord):
            return coord, atom
    return False


def calculate_nearest_neighbors(atoms: Atoms, rqstd: float, tolerance: float = 1e-2):
    """
    Calculate the nearest neighbors of each atom and their coordination numbers
    within the given request distance rqstd and tolerance.
    Use KDTree to improve efficiency and reduce memory usage.
    """
    # coords = atoms.direct_coords
    # if res := verify_atoms(coords=coords, atoms=atoms):
    #     raise ValueError(res)

    # Get cartesian coordinates
    coords = atoms.cartesian_coords

    # Make KDTree
    tree = KDTree(coords)

    nearest_neighbors = defaultdict(list)
    pair_count = 0

    for i, coord in enumerate(tqdm(coords, desc="Searching for NN", ncols=80)):
        # Search for neighbors in rqstd + tolerance
        neighbor_indices = tree.query_ball_point(coord, r=rqstd + tolerance)

        # Calculate distances
        distances = np.linalg.norm(coords[neighbor_indices] - coord, axis=1)

        # Filter out neighbors
        for j, dist in zip(neighbor_indices, distances):
            if dist < tolerance:
                continue
            elif np.abs(dist - rqstd) > tolerance:
                continue

            atom_i = atoms[i]
            atom_j = atoms[j]
            nearest_neighbors[atom_i].append(atom_j)
            if atom_j not in nearest_neighbors:
                pair_count += 1

    return nearest_neighbors, pair_count


def countCN2files(filepath: str, factors: tuple[int, int, int]):
    for sep_output in separate2files(filepath):
        poscar = SimplePoscar()
        atoms = poscar.read_poscar(sep_output)

        rqstds = np.diagonal(atoms.cell) / np.array(factors)
        if len(set(rqstds)) != 1:
            print(f"Failed to determine crystal constant, as {factors}")
            rqstd = int(input("Please input crystal constant:"))
        else:
            rqstd = rqstds[0] * math.sqrt(2) / 2

        if len([s for s, _ in atoms.symbol_count]) <= 1:
            print(f"{atoms} is pure, skipping")
            continue

        output = os.path.splitext(os.path.abspath(sep_output))[0]
        if not os.path.exists(output):
            os.makedirs(output)
            
        all_cn_data = defaultdict(list)

        for symbol, atom_list in atoms.group_atoms:
            symbol_atoms = Atoms(atoms.cell, is_direct=atoms.is_direct, atom_list=atom_list)

            # fn = os.path.join(output, f"POSCAR-{symbol}.vasp")
            # poscar.write_poscar(fn, symbol_atoms)

            neighbor_atoms = atoms.copy(clean=True)
            cn_values = []
            # for i, (proj, layer) in enumerate(group_by_direction(symbol_atoms, basis)):
            #     # fn = os.path.join(output, f"POSCAR-{symbol}-{i}.vasp")
            #     # poscar.write_poscar(fn, layer)

            #     neighbors, count = calculate_nearest_neighbors(layer, rqstd)
            #     if not neighbors:
            #         continue

            #     print(f"{symbol}-{symbol}", count)
            #     neighbor_atoms.extend(neighbors)

            # In Cube
            neighbors_map, pair_count = calculate_nearest_neighbors(symbol_atoms, rqstd)
            print(f"{symbol}-{symbol} pair count: {pair_count}")
            for neighbors in neighbors_map.values():
                cn_values.append(len(neighbors))
                neighbor_atoms.extend(neighbors)
                # print(neighbor_atoms)

            # Count Coordinate Number
            avg_cn = np.mean(cn_values)
            std_cn = np.std(cn_values)
            print(f"{symbol}: Average CN = {avg_cn:.2f}, Standard Deviation CN = {std_cn:.2f}")
            all_cn_data[symbol].extend(cn_values)

            # Save neighbor atoms to POSCAR file
            filname = os.path.join(output, f'POSCAR-neighbors-{symbol}.vasp')
            poscar.write_poscar(filname, neighbor_atoms)

        # Plot the histogram of Coordinate Numbers
        symbols, all_cn_values = zip(*all_cn_data.items())
        colors = [color_map[s] for s in symbols]
        min_cn = min(min(values) for values in all_cn_data.values())
        max_cn = max(max(values) for values in all_cn_data.values())
        plt.hist(all_cn_values, bins=np.arange(min_cn, max_cn + 2),
                 alpha=0.5, label=symbols, color=colors, edgecolor='black')
        plt.legend()
        plt.title(f"Histogram of Coordination Numbers")
        plt.xlabel("Coordinate Number")
        plt.ylabel("Frequency")
        plt.savefig(os.path.join(output, f"CN-Histogram.png"))


def __test():
    filepath = 'data/POSCAR-CoNiV-303030-r1-CoNiV.vasp'
    # direction = (0, 0, 1)
    factor = 30

    # basis = _get_basis(direction)
    countCN2files(filepath, factor)


if __name__ == "__main__":
    __test()
