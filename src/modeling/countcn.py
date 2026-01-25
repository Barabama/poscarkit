# src/modeling/countcn.py

import os
import logging
import shutil
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import matplotlib

matplotlib.use("Agg")  # Set matplotlib to use non-interactive backend
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist
from ase.neighborlist import NeighborList, natural_cutoffs

from src.modeling import color_map
from src.modeling.base import Atom, Struct, SimplePoscar
from src.utils.progress import progress


@dataclass
class CNData:
    symbols: tuple[str, str]  # center, neighbor
    center: Atom
    neighbors: list[Atom]
    cn: int  # coordination number


class CNCounter:
    def __init__(self, name: str, poscar: Path):
        self._hatch_patterns = ["//", "\\\\", "||", "--", "++", "xx", "oo", "O)", "..", "**"]
        self.name = name
        self.poscar = poscar

    def detect_cutoff(self, sample_size: int = 1000) -> float:
        """
        Detect the cutoff distance for coordination number counting.

        Args:
            sample_size: Number of atoms to sample to avoid large memory usage
        Returns:
            float: Cutoff distance
        """
        struct = self.struct
        coords = struct.get_coords(direct=False)
        # Sample for a large structure
        len_coords = len(coords)
        if len_coords > sample_size:
            logging.info(f"Sampling atoms {sample_size}/{len_coords}")
            indices = np.random.choice(len_coords, size=sample_size, replace=False)
            coords = coords[indices]

        # Calculate pairwise distances and filter out small distances
        distances = pdist(coords)
        distances = distances[distances > 0.1]

        if len(distances) == 0:
            raise ValueError("No distances found")

        # Sort distances to find the first nearest neighbors
        dists_sorted = np.sort(distances)

        # Use the nearest distance for small structure
        len_dists = len(dists_sorted)
        if len_dists < 100:
            return float(dists_sorted[0])

        # Use a more sophisticated approach for larger structure
        subset_size = min(len_dists // 20, 500)
        dists_subset = dists_sorted[:subset_size]

        # Calculate the gradient
        diffs = np.diff(dists_subset)

        # Look for the first significant change in gradient
        # This indicates the transition from 1st to further neighbors
        threshold = np.mean(diffs) + np.std(diffs)

        # Find the first where the gap is larger than the further
        index_cutoff = 0
        for i, diff in enumerate(diffs):
            if diff > threshold:
                index_cutoff = i
                break

        # If no significant change is found, use a statistical approach
        dist_cutoff = (
            np.percentile(dists_subset, 5) if index_cutoff == 0 else dists_subset[index_cutoff]
        )

        return float(dist_cutoff)

    def calculate_cn(
        self, cutoff: float, cutoff_mult: float = 1.1
    ) -> tuple[list[CNData], dict[frozenset, int]]:
        """
        Calculate each atom's coordination number within a given cutoff distance.
        Use KDTree to improve performance and reduce memory usage.

        Args:
            cutoff: Cutoff distance for coordination number counting
            cutoff_mult: Multiplier for cutoff radius
        Returns:
            Tuple[List[CNData], Dict[frozenset, int]]: List of CNData and pair counts
        """

        struct = self.struct
        coords = struct.get_coords(direct=False)

        cndata_list: list[CNData] = []
        pair_counts: dict[frozenset, int] = defaultdict(int)

        # Apply cutoff multiplier
        cutoff = cutoff * cutoff_mult

        # Make a KDTree
        tree = KDTree(coords)
        tolerance = cutoff * 0.1
        for i, coord in progress(enumerate(coords), len(coords), desc="Calculating CN"):
            # Search for neighbors within the cutoff + tolerance
            indices = tree.query_ball_point(coord, cutoff + tolerance)
            if len(indices) <= 0:
                continue

            atom_ct = struct[i]
            s_ct = atom_ct.symbol
            nn_map = defaultdict(list)

            for j in indices:
                if i == j:
                    continue
                dist = np.linalg.norm(coords[j] - coord)
                if (dist < tolerance) or (dist > cutoff + tolerance):
                    continue

                atom_nb = struct[j]
                s_nb = atom_nb.symbol
                nn_map[s_nb].append(atom_nb)

            for s_nb, nn_list in nn_map.items():
                cndata_list.append(
                    CNData(
                        symbols=(s_ct, s_nb),
                        center=atom_ct,
                        neighbors=nn_list,
                        cn=len(nn_list),
                    )
                )
                # Count pairs
                pair_counts[frozenset([s_ct, s_nb])] += 1
        pair_counts = dict(pair_counts)
        return cndata_list, pair_counts

    def calculate_cn_ase(
        self, cutoff_mult: float = 1.2
    ) -> tuple[list[CNData], dict[frozenset, int]]:
        """
        Calculate each atom's coordination number by ASE within a given cutoff multiplier.

        Args:
            cutoff_mult: Cutoff multiplier for coordination number counting
        Returns:
            Tuple[List[CNData], Dict[frozenset, int]]: List of CNData and pair counts
        """
        struct = self.struct
        atoms = SimplePoscar.struct2atoms(struct)
        cutoffs = [cutoff_mult * c for c in natural_cutoffs(atoms)]
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(atoms)

        cndata_list: list[CNData] = []
        pair_counts = defaultdict(int)
        for i in progress(range(len(atoms)), len(atoms), desc="Calculating CN"):
            neighbors, _ = nl.get_neighbors(i)
            if len(neighbors) <= 0:
                continue
            
            atom_ct = struct[i]
            s_ct = atom_ct.symbol
            nn_map = defaultdict(list)

            for j in neighbors:
                atom_nb = struct[j]
                s_nb = atom_nb.symbol
                nn_map[s_nb].append(atom_nb)

            for s_nb, nn_list in nn_map.items():
                cndata_list.append(
                    CNData(
                        symbols=(s_ct, s_nb),
                        center=atom_ct,
                        neighbors=nn_list,
                        cn=len(nn_list),
                    )
                )
                # Count pairs
                pair_counts[frozenset([s_ct, s_nb])] += 1
        pair_counts = dict(pair_counts)
        return cndata_list, pair_counts

    def generate_cn_structs(self) -> dict[tuple[str, str, int], Struct]:
        """
        Generate Struct for each (center, neighbor, coordination number).

        Returns:
            Dict[Tuple[str, str, int], Struct]: {cn_info: Struct}
        """
        struct = self.struct
        cndata_list = self.cndata_list

        # Collect cndata
        cn_structs = {}
        for cndata in progress(cndata_list, desc="Generating CN Structs"):
            s_ct, s_nb = cndata.symbols
            cn = cndata.cn
            substruct = struct.copy(atom_list=[cndata.center] + cndata.neighbors)
            cn_structs[(s_ct, s_nb, cn)] = substruct

        return cn_structs

    def generate_poscar(self, cn_structs: dict[tuple[str, str, int], Struct], outdir: Path):
        """
        Generate POSCAR files for each (center, neighbor, coordination number).

        Args:
            cn_structs: Dict{cn_info: Struct}
            outdir: Output directory
        """
        name = self.name
        for (s_ct, s_nb, cn), substruct in cn_structs.items():
            output = outdir.joinpath(f"{name}-d1nn-{s_ct}-{s_nb}-{cn}.vasp")
            SimplePoscar.write_poscar(poscar=output, struct=substruct, comment=str(output.stem))

    def generate_poscar_by_cn(
        self, cn_structs: dict[tuple[str, str, int], Struct], outdir: Path
    ):
        """
        Generate POSCAR files by coordination number.

        Args:
            cn_structs: Dict{cn_info: Struct}
            outdir: Output directory
        """
        name = self.name
        cn_substruct_list = defaultdict(list[Struct])
        for (s_ct, s_nb, cn), struct in cn_structs.items():
            cn_substruct_list[cn].append(struct)

        for cn, struct_list in cn_substruct_list.items():
            if len(struct_list) <= 1:
                continue
            struct_all = struct_list[0]
            for struct in struct_list[1:]:
                struct_all.extend(struct)
            if not struct_all:
                continue

            output = outdir.joinpath(f"{name}-d1nn-cn{cn}.vasp")
            SimplePoscar.write_poscar(
                poscar=output, struct=struct_all, comment=str(output.stem)
            )

    def save_dataframe(self, outdir: Path):
        """
        Save dataframe to CSV file.

        Args:
            outdir: Output directory
        """
        cndata_list = self.cndata_list
        # group cndata
        cndata_dict = defaultdict(list)
        for cndata in cndata_list:
            s_ct, s_nb = cndata.symbols
            cn = cndata.cn
            cndata_dict[(s_ct, s_nb, cn)].append(cndata)

        # Save to CSV
        data = defaultdict(list)
        for (s_ct, s_nb, cn), cndata_list in cndata_dict.items():
            data["CN"].append(f"{s_ct}*-{cn}{s_nb}")
            data["Count"].append(len(cndata_list))

        output = outdir.joinpath(f"{self.name}-d1nn-cn-count.csv")
        df = pd.DataFrame(data)
        df.to_csv(output, index=False)
        logging.info(f"Dataframe saved to {output}")

    def plot_histogram_faceted(self, outdir: Path):
        """
        Plot histogram faceted.

        Args:
            outdir: Output directory
        """
        name = self.name
        cndata_list = self.cndata_list
        pair_counts = self.pair_counts
        hatche_patterns = self._hatch_patterns
        symbols = sorted(list(set(d.symbols[0] for d in cndata_list)))
        nums = len(symbols)

        # Group data
        cn_stats = defaultdict(list)
        for cndata in cndata_list:
            cn_stats[cndata.symbols].append(cndata.cn)

        symbol2hatch = {
            s: hatche_patterns[i % len(hatche_patterns)] for i, s in enumerate(symbols)  #
        }

        # n * n
        fig, axes = plt.subplots(
            nums, nums, figsize=(nums * 4, nums * 4), sharex=True, sharey=True
        )
        if nums == 1:
            axes = np.array([[axes]])  # Make it 2D for consistent indexing

        for i, s_ct in enumerate(symbols):
            for j, s_nb in enumerate(symbols):
                ax = axes[i, j]
                key = (s_ct, s_nb)
                data = cn_stats.get(key, [])
                if len(data) <= 0:
                    continue

                ax.hist(
                    data,
                    bins=range(0, max(data) + 1),
                    alpha=0.7,
                    edgecolor="black",
                    hatch=symbol2hatch[s_nb],
                    align="left",
                )

                pairs = pair_counts.get(frozenset([s_ct, s_nb]), 0)
                ax.text(
                    0.95,
                    0.95,
                    f"{s_ct}-{s_nb} pairs: {pairs}",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
                )

                if i == 0:  # Set column titles
                    ax.set_title(f"Neighbor: {s_nb}")
                if j == 0:  # Set row titles
                    ax.set_ylabel(f"Center: {s_ct}")

                ax.set_xlim(0, 12)  # Max CN 12
                ax.grid(True, alpha=0.3)

        fig.suptitle("Coordination Number Distribution", fontsize=16)
        plt.tight_layout()
        output = outdir.joinpath(f"{name}-cn-histogram-faceted.png")
        plt.savefig(output)
        plt.close()
        logging.info(f"Histogram faceted saved to {output}")

    def plot_histogram_stacked(self, outdir: Path):
        """
        Plot histogram stacked.

        Args:
            outdir: Output directory
        """
        name = self.name
        cndata_list = self.cndata_list
        pair_counts = self.pair_counts
        hatch_patterns = self._hatch_patterns
        # Group data
        cn_stats = defaultdict(lambda: defaultdict(list))
        for cndata in cndata_list:
            s_ct, s_nb = cndata.symbols
            cn_stats[s_ct][s_nb].append(cndata.cn)

        # Plot stacked histogram for each center
        for s_ct, data_nb in cn_stats.items():
            s_nb_list = list(data_nb.keys())
            data = [data_nb[s_nb] for s_nb in s_nb_list]
            cn_max = max(max(d) for d in data) if data else 12
            pairs = [
                f"{s_ct}-{s_nb} pairs: {pair_counts.get(frozenset([s_ct, s_nb]), 0)}"
                for s_nb in s_nb_list
            ]

            hatches = [hatch_patterns[i % len(hatch_patterns)] for i in range(len(s_nb_list))]

            plt.figure(figsize=(10, 6))
            plt.hist(
                data,
                bins=range(0, cn_max + 1),
                label=pairs,
                alpha=0.7,
                edgecolor="black",
                hatch=hatches,
                stacked=True,
                align="left",
            )
            plt.title(f"Coordination Number of Center {s_ct}")
            plt.xlabel("Coordination Number")
            plt.ylabel("Frequency")
            plt.legend(title="Neighbor")
            plt.grid(True, alpha=0.3)
            plt.xticks(range(0, cn_max + 1))

            output = outdir.joinpath(f"{name}-cn-histogram-{s_ct}.png")
            plt.savefig(output, dpi=300, bbox_inches="tight")
            plt.close()
            logging.info(f"Coordination Number of Center {s_ct} saved to {output}")

    def plot_heatmap(self, outdir: Path):
        """
        Plot heatmap.

        Args:
            outdir: Output directory
        """
        name = self.name
        cndata_list = self.cndata_list
        # Group data
        cndata_dict = defaultdict(list)
        for cndata in cndata_list:
            cndata_dict[cndata.symbols].append(cndata)

        cn_avg = defaultdict(dict)
        for (s_ct, s_nb), cndata_list in cndata_dict.items():
            cn_avg[s_ct][s_nb] = np.mean([d.cn for d in cndata_list])

        df = pd.DataFrame(cn_avg).fillna(0)  # center 为行, neighbor 为列

        plt.figure(figsize=(8, 8))
        sns.heatmap(df, annot=True, fmt=".1f", cmap="YlGnBu", cbar_kws={"label": "Average CN"})
        plt.title("Average Coordination Number Heatmap")
        plt.xlabel("Center Atom")
        plt.ylabel("Neighbor Atom")
        plt.tight_layout()
        output = outdir.joinpath(f"{name}-cn-heatmap.png")
        plt.savefig(output, dpi=300, bbox_inches="tight")
        plt.close()
        logging.info(f"Coordination Number Heatmap saved to {output}")

    def countCN2files(
        self,
        outdir: Path | str,
        cutoff_mult: float = 1.1,
        parallel: int = 2,
        by_ase: bool = False,
    ) -> Path:
        """
        Calculate coordination number and save to files.

        Args:
            outdir: Output directory
            cutoff_mult: Multiplier for cutoff radius
            parallel: Number of parallel processes
            by_ase: Whether to use ASE for CN calculation
        Returns:
            Path: Path to output directory
        """
        name = self.name
        poscar = self.poscar
        # Create output directory
        outdir = Path(outdir) if isinstance(outdir, str) else outdir
        outdir = outdir.joinpath(f"{name}-cn-count")
        if outdir.exists():
            shutil.rmtree(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # Read POSCAR
        self.struct = SimplePoscar.read_poscar(poscar)
        logging.debug(f"Struct: {self.struct}")

        # Search CN
        logging.info(f"Applying cutoff multiplier: {cutoff_mult}")
        if by_ase:
            logging.info("Searching CN by ase")
            self.cndata_list, self.pair_counts = self.calculate_cn_ase(cutoff_mult)
        else:
            cutoff = self.detect_cutoff()
            logging.info(f"Detected cutoff: {cutoff:.2f}")
            self.cndata_list, self.pair_counts = self.calculate_cn(cutoff, cutoff_mult)

        def write_single(args: tuple[tuple[str, str, int], Struct]):
            (s_ct, s_nb, cn), substruct = args
            output = outdir.joinpath(f"{name}-d1nn-{s_ct}-{s_nb}-{cn}.vasp")
            SimplePoscar.write_poscar(output, substruct, comment=output.stem)

        def write_by_cn(args: tuple[int, list[Struct]]):
            cn, structs = args
            if len(structs) <= 1:
                return
            merged = structs[0]
            for s in structs[1:]:
                merged = merged.copy()
                merged.extend(s)
            output = outdir.joinpath(f"{name}-d1nn-cn{cn}.vasp")
            SimplePoscar.write_poscar(output, merged, comment=output.stem)

        # Generate POSCAR
        cn_structs = self.generate_cn_structs()
        cn_groups = defaultdict(list)
        for k, s in cn_structs.items():
            cn_groups[k[2]].append(s)

        avail_parallel = min(parallel, (os.cpu_count() or 1))
        with ThreadPoolExecutor(max_workers=avail_parallel) as exe:
            list(exe.map(write_single, cn_structs.items()))
            list(exe.map(write_by_cn, cn_groups.items()))

        # Save to CSV
        self.save_dataframe(outdir=outdir)

        # Plot histogram
        self.plot_histogram_faceted(outdir=outdir)
        self.plot_histogram_stacked(outdir=outdir)
        self.plot_heatmap(outdir=outdir)

        return outdir
