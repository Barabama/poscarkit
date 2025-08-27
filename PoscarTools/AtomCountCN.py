# AtomCountCN.py

import logging
import os
import shutil
from collections import defaultdict
from typing import TypedDict, Dict, List, Set, Tuple, Any, Optional

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm

from .SimplePoscar import Atom, Atoms, SimplePoscar
from .Utils import color_map


class CNData(TypedDict):
    center: Atom
    neighbors: list[Atom]
    cn: int


def detect_cutoff_distance(atoms: Atoms, symbol: Optional[str] = None, sample_size: int = 1000) -> float:
    """
    Automatically detect the first nearest neighbor distance based on distance distribution analysis.
    
    Args:
        atoms: Atoms object containing atomic positions
        symbol: Atom symbol to analyze (if None, analyze all atoms)
        sample_size: Number of atoms to sample for distance calculation (to avoid memory issues)
        
    Returns:
        Detected cutoff distance for first nearest neighbors
    """
    # If a specific symbol is given, filter atoms
    if symbol:
        filtered_atoms = atoms.copy(atom_list=[atom for atom in atoms if atom.symbol == symbol])
    else:
        filtered_atoms = atoms
    
    coords = filtered_atoms.cartesian_coords
    
    # For very large systems, use sampling to avoid memory issues
    if len(coords) > sample_size:
        logging.info(f"Sampling {sample_size} atoms from {len(coords)} total atoms for cutoff detection")
        indices = np.random.choice(len(coords), size=sample_size, replace=False)
        coords = coords[indices]
    
    # Calculate pairwise distances using pdist
    distances = pdist(coords)
    
    # Filter out very small distances (same atom or numerical error)
    distances = distances[distances > 0.1]
    
    if len(distances) == 0:
        raise ValueError("No valid distances found")
    
    # Sort distances to find the first group of similar distances (first nearest neighbors)
    sorted_distances = np.sort(distances)
    
    # For small systems, use the minimum distance
    if len(sorted_distances) < 100:
        return float(sorted_distances[0])
    
    # For larger systems, use a more sophisticated approach
    # Look at the first 5% of distances
    subset_size = min(len(sorted_distances) // 20, 500)
    subset_distances = sorted_distances[:subset_size]
    
    # Calculate the gradient to find where the rate of increase changes significantly
    diffs = np.diff(subset_distances)
    
    # Look for the first significant increase in the rate of distance increase
    # This indicates the transition from 1st to higher nearest neighbors
    threshold = np.mean(diffs) + np.std(diffs)  # Threshold for significant increase
    
    # Find the first distance where the gap to the next distance is significantly larger
    cutoff_index = 0
    for i in range(len(diffs)):
        if diffs[i] > threshold:
            cutoff_index = i
            break
    
    # If no significant gap found, use a statistical approach
    if cutoff_index == 0:
        # Use the minimum distance or a small percentile
        cutoff_distance = np.percentile(subset_distances, 5)
    else:
        cutoff_distance = subset_distances[cutoff_index]
    
    return float(cutoff_distance)


def calculate_nearest_neighbors(atoms: Atoms, cut_off: float) -> Tuple[Dict[Atom, List[Atom]], int]:
    """
    Calculate the nearest neighbors of each atom and their coordination numbers
    within the given cut-off distance.
    Use KDTree to improve efficiency and reduce memory usage.

    Args:
        atoms: Atoms object containing the atomic positions
        cut_off: Cut-off distance for determining nearest neighbors

    Returns:
        Tuple containing a dictionary mapping atoms to their neighbors and pair count
    """
    nn_map: Dict[Atom, List[Atom]] = defaultdict(list)
    pair_count = 0

    # Get cartesian coordinates
    coords = atoms.cartesian_coords
    # Make KDTree
    tree = KDTree(coords)
    for i, coord in enumerate(tqdm(coords, desc="Searching for NN", ncols=80)):
        # Search for neighbors in rqstd + tolerance
        neighbor_indices = tree.query_ball_point(coord, r=1.2 * cut_off)

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


def generate_poscar_map(atoms: Atoms, cut_off: float, outdir: str) -> Tuple[Dict[Tuple[str, int], str], Dict[str, Any], Dict[str, List[int]]]:
    """
    Generate POSCAR files for each symbol and coordination number combination.

    Args:
        atoms: Atoms object containing atomic positions
        cut_off: Cut-off distance for determining nearest neighbors
        outdir: Output directory for generated files

    Returns:
        Tuple containing:
        - Map of (symbol, cn) to filepath
        - Symbol data dictionary for CSV generation
        - Symbol coordination number frequencies for histogram
    """
    cn_file_map = {}  # map {(symbol, cn): filepath}
    symbol_df = {}
    symbol_cn_freq = {}

    for symbol, subatoms in atoms.group_atoms():
        # Get nearest neighbors
        nn_map, pair_count = calculate_nearest_neighbors(subatoms, cut_off)
        logging.info(f"{symbol}-{symbol} pair count: {pair_count}")

        # Collect CN data map
        cn_data_map = defaultdict(list)
        for center, neighbors in nn_map.items():
            cn = len(neighbors)
            cn_data_map[cn].append(CNData(center=center, neighbors=neighbors, cn=cn))

        # Get each coordinate number atoms and generate POSCAR files
        # nn_atoms_set: Set[Atom] = set()
        for cn, cn_data_list in cn_data_map.items():
            cn_atom_sets: Set[Atom] = set()
            for cn_data in cn_data_list:
                sets = {cn_data["center"], *cn_data["neighbors"]}
                cn_atom_sets.update(sets)

            # nn_atoms_set.update(cn_atom_sets)
            cn_atoms = atoms.copy(atom_list=list(cn_atom_sets))
            logging.debug(f"{symbol}*-{cn}{symbol} atoms: {cn_atoms}")

            # Save to file
            output = os.path.join(outdir, f"POSCAR-d1nn-{symbol}-{cn}.vasp")
            comment = f"CoordinationNumber-{symbol}-{cn}"
            SimplePoscar.write_poscar(filepath=output, atoms=cn_atoms, comment=comment)

            # Add to map
            cn_file_map[(symbol, cn)] = output

        # # Get nearest neighbors atoms
        # nn_atoms = atoms.copy(atom_list=list(nn_atoms_set))
        # logging.debug(f"NN atoms: {nn_atoms}")

        # # Save to file
        # output = os.path.join(outdir, f"POSCAR-d1nn-{symbol}.vasp")
        # comment = f"NearestNeighbors-{symbol}-pair_count={pair_count}"
        # SimplePoscar.write_poscar(filepath=output, atoms=nn_atoms, comment=comment)

        # Collect CN data to df
        cn_counts = [[f"{symbol}*-{cn}{symbol}", len(d)] for cn, d in cn_data_map.items()]
        cn_df = pd.DataFrame(data=cn_counts, columns=["CN", "Count"])
        symbol_df[symbol] = cn_df

        # Collect CN data to list
        cn_freq = [d["cn"] for ds in cn_data_map.values() for d in ds]
        symbol_cn_freq[symbol] = cn_freq

    return cn_file_map, symbol_df, symbol_cn_freq


def merge_by_symbol(cn_file_map: Dict[Tuple[str, int], str], outdir: str):
    """
    Merge POSCAR files by symbol.

    Args:
        cn_file_map: Map of (symbol, cn) to filepath
        outdir: Output directory for generated files
    """
    symbol_files = defaultdict(list)
    for (symbol, cn), filepaths in cn_file_map.items():
        symbol_files[symbol].append(filepaths)

    # Merge CN files for each symbol
    for symbol, filepaths in symbol_files.items():
        if len(filepaths) == 1:
            # If only one file, no need to merge
            continue

        # Merge all CN files for this symbol
        merged_atoms: Optional[Atoms] = None
        for filepath in filepaths:
            add_atoms = SimplePoscar.read_poscar(filepath)
            if merged_atoms is None:
                merged_atoms = add_atoms.copy()
            else:
                merged_atoms.extend(add_atoms)

        # Save merged atoms to file
        output = os.path.join(outdir, f"POSCAR-d1nn-{symbol}.vasp")
        comment = f"All coordination numbers for symbol {symbol}"
        if merged_atoms is not None:
            SimplePoscar.write_poscar(filepath=output, atoms=merged_atoms, comment=comment)


def merge_by_cn(cn_file_map: Dict[Tuple[str, int], str], outdir: str):
    """
    Merge POSCAR files by coordination number.

    Args:
        cn_file_map: Map of (symbol, cn) to filepath
        outdir: Output directory for generated files
    """
    cn_files = defaultdict(list)
    for (symbol, cn), filepaths in cn_file_map.items():
        cn_files[cn].append(filepaths)

    # Merge POSCAR files for each coordination number
    for cn, filepaths in cn_files.items():
        if len(filepaths) == 1:
            # If only one file, no need to merge
            continue

        # Merge CN files for this coordination number
        merged_atoms: Optional[Atoms] = None
        for filpath in filepaths:
            add_atoms = SimplePoscar.read_poscar(filpath)
            if merged_atoms is None:
                merged_atoms = add_atoms.copy()
            else:
                merged_atoms.extend(add_atoms)

        # Save merged atoms to file
        output = os.path.join(outdir, f"POSCAR-d1nn-all-{cn}.vasp")
        comment = f"All atoms with coordination number {cn}"
        if merged_atoms is not None:
            SimplePoscar.write_poscar(filepath=output, atoms=merged_atoms, comment=comment)


def plot_histogram(symbol_cn_freq: Dict[str, List[int]], outdir: str):
    """
    Plot coordination number histogram.

    Args:
        symbol_cn_freq: Dictionary mapping symbols to their coordination number frequencies
        outdir: Output directory for generated files
    """
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
    output = os.path.join(outdir, "cn-histogram.png")
    plt.savefig(output)
    plt.close()
    logging.info(f"Coordination Number histogram saved to {output}")


def countCN2files(filepath: str, outdir: str) -> str:
    """
    Count coordination numbers and generate POSCAR files for each coordination type.

    Args:
        filepath: Path to the input POSCAR file
        outdir: Output directory for generated files

    Returns:
        Path to the directory containing the output files
    """
    # Read POSCAR
    atoms = SimplePoscar.read_poscar(filepath)
    logging.debug(f"Atoms: {atoms}")

    # detect cutoff distance
    cut_off = detect_cutoff_distance(atoms)
    logging.info(f"Automatically detected cutoff distance: {cut_off:.3f} Å")

    # Make output directory
    outdir = os.path.join(outdir, f"{os.path.splitext(os.path.basename(filepath))[0]}-cn")
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir, exist_ok=True)

    # Generate POSCAR map
    cn_file_map, symbol_df, symbol_cn_freq = generate_poscar_map(
        atoms=atoms, cut_off=cut_off, outdir=outdir)

    # Merge files by symbol
    merge_by_symbol(cn_file_map=cn_file_map, outdir=outdir)

    # Merge files by coordination number
    merge_by_cn(cn_file_map=cn_file_map, outdir=outdir)

    # Write CN data to CSV
    all_df = pd.concat(symbol_df.values(), ignore_index=True)
    output = os.path.join(outdir, "cn-counts.csv")
    all_df.to_csv(output, index=False)
    logging.info(f"Coordination Number counts saved to {output}")

    # Plot histogram
    plot_histogram(symbol_cn_freq, outdir)

    return outdir
