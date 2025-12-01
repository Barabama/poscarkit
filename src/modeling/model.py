# src/modeling/model.py

import re
import random
import logging
from collections import defaultdict
from pathlib import Path
from typing import Any, Generator

import numpy as np
from sqsgenerator import StructureFormat, parse_config, optimize
from sqsgenerator.core import SqsConfiguration

from src.modeling.base import SimplePoscar, Struct
from src.modeling.supercell import make_supercell


def _integer_fractions(
    fracts: dict[str, float], factors: tuple[int, int, int], multipl: int
) -> dict[str, int]:
    """
    Convert decimal fractions to integer fractions.

    Args:
        fracts: Dictionary {symbol: decimal fraction}
        factors: Factors for supercell (nx, ny, nz)
        multip: Multiplicity of sublattice in unitcell

    Returns:
        dict: Dictionary {symbol: integer fraction}.
    """
    factor = np.prod(factors)
    # Super_fracts is every fract * multi * factor
    super_fracts = {s: f * multipl * factor for s, f in fracts.items()}
    # Round the values to get integers
    fracts_rounded = {s: round(f) for s, f in super_fracts.items()}
    rounded_total = sum(fracts_rounded.values())
    target_total = multipl * factor

    # Adjusting rounding errors
    if rounded_total != target_total:
        # Calculate differences and adjust based on closeness to the next integer
        diffs = [
            (s, float(abs(f - round(f))), 1 if f - round(f) > 0 else -1)
            for s, f in super_fracts.items()
        ]
        diffs.sort(key=lambda x: x[1], reverse=True)
        # Determine the adjustment needed
        adjustment = target_total - rounded_total

        # Apply adjustments
        for i in range(abs(adjustment)):
            symbol, decimal, direction = diffs[i]
            fracts_rounded[symbol] += direction

    return fracts_rounded


def _ask_normalize_fractions(site: str, fractions: dict[str, float]):
    """Ask user to normalize the fractions."""
    total = sum(fractions.values())
    if abs(total - 1) > 1e-6:
        logging.warning(f"The Sum(fractions of site {site}) not close to 1.")
        logging.warning(f"Fractions: {fractions}")

        # Ask user to normalize the fractions
        user_answer = input("Normalize the fractions? [Y/N]|Y\n>>> ")
        if user_answer.lower().startswith("y"):
            fractions = {s: f / total for s, f in fractions.items()}
        else:
            raise ValueError(f"Sum(fractions of site {site}) not close to 1, stop!")
    return fractions


class ModelStruct:

    def __init__(
        self,
        name: str,
        poscar: Path,
        factors: tuple[int, int, int],
        struct_info: dict[str, Any],
        batch_size: int,
    ):
        """
        Args:
            name: Name of the model
            factors: Supercell factors (nx, ny, nz)
            struct: class Struct to be modeled
            struct_info: Dict including call and site fractions
                e.g. {
                    "cell": [2.7, 2.7, 2.7] for BCC or 3x3 matrix for HEX,
                    "1a": {
                        "atoms": ["Al", [[0, 0, 0], ... ]],
                        "sofs": {"Co": 0.5, "Ni": 0.5, ... },
                    },
                    "1b": {
                        "atoms": ["Ni", [[0.5, 0.5, 0.5], ... ]],
                        "sofs": {"Co": 0.5, "Ni": 0.5, ... },
                    },
                    ... for other sublattices
        }
        """
        self.name = name
        self.factors = factors
        self.poscar = poscar
        self.structure_info = struct_info
        self.batch_size = batch_size
        self.unitcell = SimplePoscar.read_poscar(poscar)
        self.site_integers = self._gen_site_integers()

    def _gen_site_integers(self):
        """
        Generate site integers.

        Returns:
            site_integers: Dict{(site, elem): {symbol: integer fraction}}
        """
        factors = self.factors
        structure_info = self.structure_info
        unitcell = self.unitcell

        site_integers = defaultdict(dict)
        struct_info = structure_info.copy()
        if not struct_info:
            logging.warning("No structure info provided.")
            return site_integers

        symbol_count = {s: c for s, c in unitcell.symbol_count}
        for site, data in struct_info.items():
            # Continue cell
            if site == "cell":
                continue
            # Get site fractions
            if "sofs" not in data:
                raise ValueError(f"SOFs_info not found in stie {site}.")
            fractions: dict = data["sofs"]
            fractions = _ask_normalize_fractions(site=site, fractions=fractions)

            # # Get multiplicity from site identifier
            # site_integers[site] = _integer_fractions(
            #     fracts=fractions, factors=factors, multipl=int(site[0])
            # )

            # Get multiplicity from length of atoms
            if "atoms" not in data:
                raise ValueError(f"Atoms_info not found in stie {site}.")
            elem, coords = data["atoms"]
            if elem not in symbol_count:
                raise ValueError(f"Element {elem} not found in struct {unitcell}.")
            site_integers[(site, elem)] = _integer_fractions(
                fracts=fractions,
                factors=factors,
                multipl=symbol_count[elem],
            )
        logging.info(f"SOFs_info: {dict(site_integers)}")
        return site_integers

    def _allocate_atoms(
        self,
        supercell: Struct,
        site_integers: dict[tuple[str, str], dict[str, int]] = {},
        seed: int | None = None,
    ) -> dict[str, Struct]:
        """
        Allocate atoms according to site of intergers.

        Args:
            supercell: Supercell Struct
            site_integers: Dict {site: {symbol: intergers}
            seed: Random seed to shuffle for reproducibility

        Returns:
            dict: Dict {site: Allocated substruct}.
        """
        if seed is not None:
            random.seed(seed)

        site_substruct = {}
        for note, substruct in supercell.group_structs(key="note"):
            match = re.search(r"(\d+[a-z])-([A-Za-z]+)", note)
            if not match:
                raise ValueError(
                    f"Unknwon note {note} to recognize the sublattice in {substruct}."
                )
            site, elem = match.groups()
            site = (site, elem)

            # Re-index atoms
            for i, atom in enumerate(substruct):
                atom.index = i

            # Shuffle atoms within the same sublattice
            sub_list = substruct.atom_list
            random.shuffle(sub_list)

            if not site_integers:
                symbols = [atom.symbol for atom in sub_list]
            elif site in site_integers:
                symbols = [s for s, f in site_integers[site].items() for i in range(f)]
            else:
                raise ValueError(
                    f"Site {site} not found in site_integers {site_integers}."
                )

            # Assign symbols and meta
            slsl = len(str(len(sub_list)))
            for idx, (symbol, atom) in enumerate(zip(symbols, sub_list), start=1):
                atom.symbol = symbol
                atom.meta = f"{idx:0{slsl}d}"

            new_substruct = substruct.copy(atom_list=sub_list)
            site_substruct[site] = new_substruct

        return site_substruct

    def model_by_shuffle(
        self, seeds: list[int | None] = [None]
    ) -> Generator[tuple[Path, Struct], None, None]:
        """
        Model the structure by shuffling sublattices.

        Args:
            batch_size: Batch size
            seed: Random seed

        Yields:
            Tuple(filename, Struct)
        """
        name = self.name
        factors = self.factors
        unitcell = self.unitcell
        site_integers = self.site_integers
        supercell = make_supercell(unitcell, factors)
        batch_size = self.batch_size
        if len(seeds) < batch_size:
            seeds += [None] * (batch_size - len(seeds))
        ssl = len(str(len(seeds)))
        for t, seed in enumerate(seeds, start=1):
            logging.info(f"Modeling {name} by shuffle, batch {t}/{batch_size}")

            shuffled = supercell.copy(clean=True)
            site_substruct = self._allocate_atoms(
                supercell=supercell, site_integers=site_integers, seed=seed
            )
            for (site, elem), substruct in site_substruct.items():
                # Get every site of substruct
                filename = Path(f"{name}-shuffle{t:0{ssl}d}-{site}.vasp")
                yield filename, substruct
                shuffled.extend(substruct)

            # Get all sites of newstruct
            filename = Path(f"{name}-shuffle{t:0{ssl}d}.vasp")
            yield filename, shuffled

    def model_by_sqsgen(
        self, iterations: int = 10000000
    ) -> Generator[tuple[Path, Struct], None, None]:
        """
        Model the structure by sqsgenerator.

        Args:
            Iterations
        Yields:
            Tuple(filename, Struct)

        """
        name = self.name
        factors = self.factors
        unitcell = self.unitcell
        batch_size = self.batch_size
        site_integers = self.site_integers

        if not site_integers:
            logging.error("SQS modeling requires site_integers from structure_info")
            return

        seeds = range(batch_size)
        ssl = len(str(len(seeds)))
        for t, seed in enumerate(seeds, start=1):
            logging.info(f"Modeling {name} by sqsgen, batch {t}/{batch_size}")

            sqsgen_list: list[Struct] = []
            for (site, elem), comps in site_integers.items():
                struct_site = unitcell.classify(elem)
                json_info = {
                    "iterations": int(iterations),
                    "sublattice_mode": "split",
                    "structure": {
                        "lattice": struct_site.cell.tolist(),
                        "coords": struct_site.get_coords().tolist(),
                        "species": struct_site.symbols,
                        "supercell": list(factors),
                    },
                    "composition": [{**comps, "sites": elem}],
                }
                sqscfg = parse_config(dict(json_info))
                if not isinstance(sqscfg, SqsConfiguration):
                    raise ValueError(f"Error parsed sqscfg {sqscfg}.")

                logging.info(f"Sqs generating for site {site}...")
                pack = optimize(sqscfg)
                if len(pack) < 1:
                    raise ValueError("Sqs result empty.")

                best = pack.best()
                sqs_struct = best.structure()
                subsqsgen = SimplePoscar.read_poscar(
                    sqs_struct.dump(StructureFormat.poscar)
                )
                filename = Path(f"{name}-sqsgen{t:0{ssl}d}-{site}.vasp")
                sqsgen_list.append(subsqsgen)
                yield filename, subsqsgen

            if len(sqsgen_list) <= 1:
                continue
            merged = sqsgen_list[0]
            for subsqsgen in sqsgen_list[1:]:
                merged.extend(subsqsgen)
            filename = Path(f"{name}-sqsgen{t:0{ssl}d}.vasp")
            yield filename, merged
