# src/modeling/slice.py

import logging
import shutil
from collections import defaultdict
from itertools import groupby
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")  # Set matplotlib to use non-interactive backend
import matplotlib.pyplot as plt
from ase.build.tools import cut

from src.modeling import color_map
from src.modeling.base import Struct, SimplePoscar
from src.utils.progress import progress


class Slicer:

    def __init__(self, name: str, poscar: Path, miller_index: tuple[int, int, int]):
        self.basis_map = {
            (0, 0, 1): [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
            (1, 1, 0): [(0, 0, -1), (-1, 1, 0), (1, 1, 0)],
            (1, 1, 1): [(1, 1, -2), (-1, 1, 0), (1, 1, 1)],
        }
        self.name = name
        self.poscar = poscar
        self.miller_index = miller_index
        self.basis = self.get_basis()
        self.basis_norm = np.array([self.normalize(v) for v in self.basis])
        self.struct = SimplePoscar.read_poscar(poscar)
        self.transformed = self.transform_struct()

    def normalize(self, vector: np.ndarray) -> np.ndarray:
        return vector / np.linalg.norm(vector)

    def get_basis(self) -> np.ndarray:
        """
        Get basis vectors by miller index.
        Returns:
            np.ndarray: Basis vectors
        """
        miller_index = self.miller_index
        if miller_index in self.basis_map:
            basis = [np.array(v) for v in self.basis_map[self.miller_index]]
        else:
            n = np.array(miller_index)
            # Find other two vectors vertical to n
            t0 = np.array([1, 0, 0]) if abs(n[0]) < abs(n[1]) else np.array([0, 1, 0])
            b1 = np.cross(n, t0)
            b2 = np.cross(n, b1)
            basis = [b1, b2, n]

        basis = np.array(basis)
        return basis

    def transform_struct(self) -> Struct:
        """
        Transform a Struct with basis.

        Returns:
            Struct: Transformed Struct
        """
        struct = self.struct
        basis_norm = self.basis_norm
        atoms = SimplePoscar.struct2atoms(struct)
        a, b, c = basis_norm
        transformed = cut(atoms, a, b, c)
        newstruct = SimplePoscar.atoms2struct(transformed)
        return newstruct

    def group_by_normal(self, precision: int = 2):
        """
        Group atoms by distance to the normal vector.

        Args:
            precision: Precision of distance
        Yields:
            Tuple: (Projection, Struct)
        """
        transfd = self.transformed
        basis_norm = self.basis_norm
        # Calculate projections rounded
        coords = transfd.get_coords(direct=False)
        projs = np.dot(coords, basis_norm[2])  # Projection onto normal
        projs_rounded = np.round(projs, precision)

        # Group atoms
        indics_sorted = np.argsort(projs_rounded)
        for proj, group in groupby(indics_sorted, key=lambda x: projs_rounded[x]):
            layer = transfd.copy(atom_list=[transfd[i] for i in group])
            yield proj, layer

    def plot_layer(
        self,
        imgpath: Path,
        title: str,
        layer: Struct,
        pair_counts: dict[str, int] = {},
    ):
        """
        Plot layer structure based on basis.

        Args:
            imgpath: Path to image file
            title: Title of plot
            layer: Struct object
            pair_counts: Pair counts
        """
        basis = self.basis
        # Calculate projections
        layer.sort()
        coords = layer.get_coords(direct=False)
        b1, b2, n = [self.normalize(v) for v in layer.cell]
        projs_normal = np.dot(coords, n)  # projection onto normal
        projs_plane = coords - np.outer(projs_normal, n)  # projection onto plane
        xs = np.dot(projs_plane, b1)  # portion of x-axis
        ys = np.dot(projs_plane, b2)  # portion of y-axis
        coords_proj = np.column_stack((xs, ys))

        # Group by projections
        symbol_counts = defaultdict(list)
        for atom, coord in zip(layer, coords_proj):
            symbol_counts[atom.symbol].append(coord)

        # Get range of basis vectors
        x_min, x_max = 0.0, np.linalg.norm(layer.cell[0])
        y_min, y_max = 0.0, np.linalg.norm(layer.cell[1])
        x_margin = (x_max - x_min) * 0.1
        y_margin = (y_max - y_min) * 0.1

        # Plot
        plt.figure(figsize=(8, 8))
        for symbol, coords in symbol_counts.items():
            color = color_map.get(symbol, "#FF00FF")
            x, y = zip(*coords)
            label = (
                symbol
                if not (pair_counts and symbol in pair_counts)
                else f"{symbol}-{symbol} pairs: {pair_counts[symbol]}"
            )
            plt.scatter(x, y, marker="o", s=10, color=color, alpha=1.0, label=label)

        plt.title(title)
        plt.xlabel(f"[{' '.join(str(v) for v in basis[0])}] Coordinate (Å)")
        plt.ylabel(f"[{' '.join(str(v) for v in basis[1])}] Coordinate (Å)")
        # plt.axis("equal")
        plt.grid()
        plt.legend(title="Symbols", bbox_to_anchor=(1, 1), loc="upper left")
        plt.xlim(-x_margin, x_max + x_margin)
        plt.ylim(-y_margin, y_max + y_margin)
        # plt.tight_layout(rect=(0, 0, 1, 1))
        plt.savefig(imgpath, bbox_inches="tight")
        plt.close()

    def slice2files(self, outdir: Path | str) -> list[Path]:
        """
        Slice a structure and save to files.

        Args:
            outdir: Output directory
        Returns:
            List[Path]: Path to layer files
        """
        name = self.name
        basis = self.basis
        miller_index = self.miller_index
        transfd = self.transformed
        # Get basis
        logging.info(f"Basis: {basis}")

        # Output directory
        miller_index_str = "".join(str(d) for d in miller_index)
        outdir = Path(outdir) if isinstance(outdir, str) else outdir
        outdir = outdir.joinpath(f"{name}-sliced({miller_index_str})")
        if outdir.exists():
            shutil.rmtree(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        # Transform
        output = outdir.joinpath(f"Transformed({miller_index_str}).vasp")
        SimplePoscar.write_poscar(poscar=output, struct=transfd, comment=str(output.stem))

        # Group by normal
        layers = [ls for ls in self.group_by_normal()]
        num_layers = len(layers)
        logging.info(f"Number of layers: {num_layers}")
        ll = len(str(num_layers))
        results = []
        for i, (proj, layer) in progress(
            enumerate(layers, start=1),
            total=num_layers,
            desc="Slicing layers",
        ):
            logging.info(f"Layer {i}: {proj:.2f}")
            logging.info(f"Layer: {layer}")

            # Save layer
            output = outdir.joinpath(f"Transformed({miller_index_str})-layer{i:0{ll}d}.vasp")
            SimplePoscar.write_poscar(poscar=output, struct=layer, comment=str(output.stem))
            results.append(output)
            # Plot layer
            imgpath = output.with_suffix(".png")
            self.plot_layer(imgpath=imgpath, title=str(imgpath.stem), layer=layer)

        logging.info(f"Sliced Saved to {outdir}")
        return results
