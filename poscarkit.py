# poscarkit.py

import argparse
import logging
import os
import sys
import tomllib
import traceback
from copy import deepcopy
from dataclasses import dataclass, field

import PoscarTools.Utils
from PoscarTools.AtomSupercell import unitcell2file, supercell2file
from PoscarTools.AtomAllocate import allocate2files
from PoscarTools.AtomCountCN import countCN2files
from PoscarTools.AtomSlice import slice2file


INFO_EXEC = f"""
============================= POSCARKIT (v0.8.0) ==============================
This toolkit has many Functions such as Making Supercell, Allocating atoms
based on SOFs, Counting Coordinate Numbers between same type atoms, Slicng,
developed by FZU-MCMF, main developers: Gao Min-Liang, Wu Bo*, Qiao Yang,
Yang Su-wen, et al..
Please contact via wubo@fzu.edu.cn, 654489521@qq.com in case need further guide
and information.
"""
INFO_CHOICES = """
========================== Options with instructions ==========================
  1) Help            : Instructions of this toolkit.
  2) Read config     : from <config.toml> as changing <Structure> or <SOFs>.
  3) Run Workflow    : First [Make Supercell] and then [Allocate Atoms].
  4) Make Supercell  : based on <SupercellFactors> along the basis vectors.
  5) Allocate Atoms  : based on <SupercellFactors> <Structure> <ShuffleSeeds>.
  6) Count CN        : CN(Coordinate Numbers) and NN(Nearest Neighbors).
  7) Slice to layers : normal to <SliceDirection> also known as Miller Index.
"""
INFO_HELP = f"""
============================ How to use POSCARKIT =============================
  Note: No need to configure everything. Just write down SOFs data and run the
        program and you will be asked to input the information.

  2) Read config    : Read the configuration <config.toml> again as changing:
    <Filepath> ---------- File path of POSCAR to be processed.
    <outdir>   ---------- Output directory of results.
    <SupercellFactors> -- 3-int factors alone basis to make supercell.
    <Structure> --------- Specify the structure to be loaded.
    <ShuffleSeeds> ------ Seeds for shuffle when allocating atoms.
    <SliceDirection> ---- Direction of slicing.
    <StructureInfo> ----- include Structure information and SOFs data.

  3) Run Workflow    : Workflow is set to do First [Make Supercell] and then
    [Allocate Atoms]. See 4) and 5) for more details.

  4) Make Supercell  : Choose a POSCAR or the unitcell of prescribed prototype
    built by config, and then expand it to a supercell based on
    <SupercellFactors> (such as (3x3x3), (30x30x30) for Cube; (2x2x3) for Hex).

  5) Allocate Atoms  : Shuffle based on <ShuffleSeeds> and Allocate atoms based
    on <SupercellFactors> and SOFs data of <Structure>. These works will be
    grouped by sublattices, output each sublattice POSCARs and a whole POSCAR.
    If not <Structure> provided, Shuffle atoms only.

  6) Count CN        : Count CN(Coordinate Numbers) of all the same type atoms
    in the supercell and calculate the NN(Nearest Neighbors). Mi*-Mi as d1NN
    (the Distance of the First Nearest Neighbors of the same type atom Mi),
    the CN of an atom Mi coordinated with its own type.

  7) Slice to layers : Slice atoms to get all the single layers normal to
    <SliceDirections> (such as [001], [110], [111]) and plot the layers.
"""
CHOICES = (1, 2, 3, 4, 5, 6, 7)

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
workdir = os.path.dirname(sys.argv[0])


@dataclass
class Config:
    # Filepath of POSCAR to operate
    Filepath: str = field(default="")

    # Output directory
    Outdir: str = field(default=os.path.join(workdir, "out"))

    # Used for 'Supercell', 'Allocation'
    SupercellFactors: list[int] = field(default_factory=lambda: [])

    # Used for 'Allocation'
    Structure: str = field(default="")

    # Used for 'Allocation'
    ShuffleSeeds: list = field(default_factory=lambda: [None])

    # Used for 'Slice
    SliceDirection: list[int] = field(default_factory=lambda: [])

    FCC: dict[str, dict[str, list | dict]] = field(default_factory=dict)
    BCC: dict[str, dict[str, list | dict]] = field(default_factory=dict)
    HCP: dict[str, dict[str, list | dict]] = field(default_factory=dict)


class PoscarKit:
    def __init__(self):
        self.config = self.read_config()
        self._func_map = {
            1: self.show_help,
            2: self.read_config,
            3: self.handle_workflow,
            4: self.handle_supercell,
            5: self.handle_allocate,
            6: self.handle_countCN,
            7: self.handle_slice,
        }

    def show_help(self, **kwargs) -> None:
        print(INFO_HELP)

    def read_config(self, **kwargs) -> Config:
        cfg_path = os.path.join(workdir, "config.toml")
        if not os.path.isfile(cfg_path):
            logging.warning("config.toml not found! Using default configuration.")
            with open(cfg_path, "w", encoding="utf-8") as tf:
                tf.write(PoscarTools.Utils.default_config)
            self.config = self.read_config(**kwargs)
        else:
            with open(cfg_path, "rb") as tf:
                self.config = Config(**tomllib.load(tf))
                logging.info(f"Loaded configuration from {cfg_path}")
        return self.config

    def _handle_option(self, option: int = 0):
        """Check option, passed option==2."""
        while not option:
            try:
                option = int(input(f"{INFO_CHOICES}Enter choice >>> "))
                if option not in CHOICES:
                    raise ValueError("Invalid option")
            except ValueError:
                logging.error(f"Invalid option({option}). Please try again.")
                option = 0
        return option

    def _handle_filepath(self, filepath: str = "", force: bool = True) -> str:
        """Check filepath."""
        if force:
            prompt = "Enter Filepath >>> "
            filepath = filepath or input(prompt)
            while not os.path.isfile(os.path.abspath(filepath)):
                logging.error(f"No such file: {filepath}")
                filepath = input(prompt)
        else:
            prompt = "Enter Filepath or NONE to use StructureInfo >>> "
            filepath = filepath or input(prompt)
            if filepath:
                filepath = self._handle_filepath(filepath=filepath, force=True)
        return filepath

    def _handle_outdir(self) -> str:
        """Check output directory."""
        outdir = self.config.Outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        return outdir

    def _handle_factors(self) -> tuple[int, int, int]:
        """Check supercell factors."""
        prompt = "Enter SupercellFactors (x y z, default=3 3 3) >>> "
        factors = tuple(self.config.SupercellFactors)
        while True:
            try:
                user_input = input(prompt).split()[:3]
                factors = factors or tuple(map(int, user_input)) or (3, 3, 3)
                if len(factors) != 3 or any(f <= 0 for f in factors):
                    raise ValueError(f"Needs 3 positive integers. Got {factors}")
                return factors
            except ValueError as e:
                logging.warning(f"{e}. Please try again.")
                factors = ()

    def _handle_miller(self) -> tuple[int, int, int]:
        """Check slice direction."""
        prompt = "Enter SliceDirection (x y z, default=0 0 1) >>> "
        miller = tuple(self.config.SliceDirection)
        while True:
            try:
                user_input = input(prompt).split()[:3]
                miller = miller or tuple(map(int, user_input)) or (0, 0, 1)
                if len(miller) != 3 or any(not isinstance(i, int) for i in miller):
                    raise ValueError(f"Needs 3 integers. Got {miller}.")
                return miller
            except ValueError as e:
                logging.warning(f"{e}. Please try again.")
                miller = ()

    def _handle_structure(self, force: bool = True) -> dict[str, dict] | None:
        """return a structure (FCC, BCC, HCP) information from the configuration."""
        structure = self.config.Structure.upper()
        if force:
            prompt = "To build Unitcell, Enter Structure (fcc, bcc, hcp, ...) >>> "
            structure = structure or input(prompt).upper()
            while structure not in self.config.__dict__:
                logging.warning(f"Structure({structure}) not found in config. Please try again.")
                structure = input(prompt).upper()
        else:
            prompt = "To load SOFs data, Enter Structure (fcc, bcc, hcp, ...) or NONE to shuffle only >>> "
            structure = structure or input(prompt).upper()
            if structure not in self.config.__dict__:
                logging.warning(f"Structure({structure}) not found in config. Shuffle only.")
        return deepcopy(self.config.__dict__.get(structure))

    def _handle_seeds(self) -> list[int | None]:
        """Check shuffle seeds."""
        seeds = self.config.ShuffleSeeds
        for seed in seeds:
            if not (isinstance(seed, int) or seed is None):
                seed = ord(seed)
        return seeds

    def handle_supercell(self, filepath: str) -> str:
        filepath = self._handle_filepath(filepath=filepath, force=False)
        outdir = self._handle_outdir()
        if not filepath:
            logging.warning("No filepath provided. Using struture info instead.")
            struct_info = self._handle_structure(force=True)
            filepath = unitcell2file(struct_info=struct_info, outdir=outdir)
        else:
            filepath = self._handle_filepath(filepath)
        factors = self._handle_factors()
        supercell = supercell2file(filepath=filepath, outdir=outdir, factors=factors)
        return supercell

    def handle_allocate(self, filepath: str) -> list[str]:
        filepath = self._handle_filepath(filepath=filepath)
        outdir = self._handle_outdir()
        factors = self._handle_factors()
        struct_info = self._handle_structure(force=False)
        seeds = self._handle_seeds()
        allocateds = allocate2files(filepath=filepath, outdir=outdir, factors=factors,
                                    struct_info=struct_info, seeds=seeds)
        return allocateds

    def handle_countCN(self, filepath: str) -> str:
        filepath = self._handle_filepath(filepath)
        outdir = self._handle_outdir()
        return countCN2files(filepath=filepath, outdir=outdir)

    def handle_slice(self, filepath: str) -> str:
        filepath = self._handle_filepath(filepath)
        outdir = self._handle_outdir()
        miller_index = self._handle_miller()
        return slice2file(filepath=filepath, outdir=outdir, miller_index=miller_index)

    def handle_workflow(self, filepath: str):
        supercell = self.handle_supercell(filepath)
        allocateds = self.handle_allocate(supercell)
        return allocateds

    def run(self, filepath: str = "", option: int = 0):
        """Main function."""
        print(INFO_EXEC)
        filepath = filepath or self.config.Filepath
        while True:
            try:
                option = self._handle_option(option)
                if option not in self._func_map:
                    continue
                func = self._func_map[option]
                result = func(**{"filepath": filepath})
                # Reset
                filepath = ""
                option = 0

            except ValueError as e:
                traceback.print_exc()
                logging.error(f"Invalid input: {e}")
                filepath = ""
                option = 0
            except KeyboardInterrupt:
                print("\n")
                print("Thank you for using POSCAR tool! Exiting...")
                sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Poscar Tool")
    parser.add_argument("filepath", type=str, nargs="?", default="", help="Poscar file path")
    parser.add_argument("-f", "--file", type=str, default="", help="Poscar file path")
    parser.add_argument("-c", "--choice", type=int, default=0, choices=CHOICES,
                        help="Choice of operation")
    args = parser.parse_args()

    try:
        poscarkit = PoscarKit()
        poscarkit.run(filepath=args.filepath or args.file, option=args.choice)

    except Exception as e:
        logging.critical(e)
        traceback.print_exc()
        input("Press Enter to exit...")
        sys.exit(1)
