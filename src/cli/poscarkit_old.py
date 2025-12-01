# src/cli/poscarkit.py

import os
import sys
import argparse
import logging
import tomllib
import traceback
from copy import deepcopy
from pathlib import Path
from typing import TypedDict

from src.modeling.base import SimplePoscar
from src.modeling.countcn import CNCounter
from src.modeling.slice import Slicer
from src.modeling.supercell import unitcell2file, supercell2file
from src.workflow.modeling import run_modeling
from src.workflow.slice_to_countcn import slice2files_with_countcn
from src.config import VERSION, CONTACT, DEVELOPER, DEFAULT_CONFIG

INFO_CLI = f"""
================================== POSCARKIT ==================================
This toolkit has many Functions such as Making Supercell, Allocating atoms
based on SOFs, Counting Coordinate Numbers between same type atoms, Slicing.
{DEVELOPER}
{VERSION}
{CONTACT}
"""
INFO_OPTIONS = """
========================== Options with instructions ==========================
  1) Help            : Instructions of this toolkit.
  2) Read config     : from <config.toml> as changing <phase> or <sofs>.
  3) Run Modeling    : First [Make Supercell] and then [Allocate Atoms].
  5) Make Supercell  : based on <supercell_factors> along the basis vectors.
  6) Count CN        : CN(Coordinate Numbers) and NN(Nearest Neighbors).
  7) Slice to layers : normal to <slice_direction> also known as Miller Index.
  8) Slice to CountCN: Slice to layers and then count CN for each layer.
"""
INFO_HELP = f"""
============================= How to use POSCARKIT ============================
  Note: No need to configure everything. Just write down SOFs data and run the
        program and you will be asked to input the information.

  2) Read config     : Read the configuration <config.toml> again as changing:
    <name> -------------   Name of this work.
    <poscar> ------------- POSCAR file path to be processed.
    <outdir>   ----------- Output directory of results.
    <phase> -------------- Specify the phase to be loaded.
    <supercell_factors> -- 3-int factors alone basis to make supercell.
    <slice_direction> ---- Direction of slicing.
    <structure_info> ----- Including Structure information and SOFs data.
    <shuffle_seeds> ------ Seeds for shuffle when allocating atoms.
    <batch_size> --------- Batch size for modeling.
    <enable_sqs> --------- Enable SQS generation.

  3) Run Modeling    : Modeling is set to do First [Make Supercell] and then
    shuffle Atoms. Shuffle based on <shuffle_seeds> and Allocate atoms based
    on <supercell_factors> and SOFs data of <phase>. These works will be
    grouped by sublattices, output each sublattice POSCARs and a whole POSCAR.
    If not <phase> provided, Shuffle atoms only.
    
  5) Make Supercell  : Choose a POSCAR or the unitcell of prescribed prototype
    built by config, and then expand it to a supercell based on
    <supercell_factors> (such as (3x3x3), (30x30x30) for Cube; (2x2x2) for Hex).

  6) Count CN        : Count CN(Coordinate Numbers) of all the pairs of atoms
    in the supercell and calculate the NN(Nearest Neighbors). A*-B as d1NN
    (the Distance of the First Nearest Neighbors of the same type atom),
    the CN of symbol A coordinated with symbol B.

  7) Slice to layers : Slice atoms to get all the single layers normal to
    <slice_direction> (such as [001], [110], [111]) and plot the layers.
    
  8) Slice to CountCN: First slice atoms to layers, then count coordination
    numbers for each layer separately. This combines functions 7 and 8.
"""


OPTIONS = (1, 2, 3, 5, 6, 7, 8)

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
workdir = Path(sys.argv[0]).absolute().parent


class Config(TypedDict):
    name: str
    poscar: str
    outdir: str
    phase: str
    supercell_factors: list[int]
    slie_direction: list[int]
    shuffle_seeds: list[int | None]
    batch_size: int
    enable_sqs: bool


class PoscarKit:

    def __init__(self, cfg_path: str = "", cfg: dict = {}):
        self._func_map = {
            1: self.show_help,
            2: self.read_config,
            3: self.handle_modeling,
            5: self.handle_supercell,
            6: self.handle_countcn,
            7: self.handle_slice,
            8: self.handle_slice_to_countcn,
        }
        self.config = self.read_config(cfg_path, cfg)

    def show_help(self, **kwargs):
        print(INFO_HELP)

    def read_config(self, cfg_path: str = "", cfg: dict = {}, **kwargs) -> Config:
        """
        Read config file
        Args:
            cfg_path: config.toml path
            cfg: Dict to parse
        Returns:
            Config: config
        """
        if cfg:
            return Config(**cfg)
        path = workdir.joinpath("config.toml")
        # path = Path(cfg_path).absolute()
        if not path.is_file():
            logging.warning(
                f"Config file {str(path)} does not exist, using default one."
            )
            with open(path, "w", encoding="utf-8") as tf:
                tf.write(DEFAULT_CONFIG)
            config = self.read_config(cfg_path=str(path), **kwargs)
        elif path.is_file():
            with open(path, "r", encoding="utf-8") as tf:
                config = Config(**tomllib.loads(tf.read()))
                logging.info(f"Read config from {str(path)}")
        return config

    def handle_option(self, option: int = 0) -> int:
        """
        Check if option in CHOICES

        Args:
            option: Option
        Returns:
            int: Option
        """
        try:
            option = option or int(input(f"{INFO_OPTIONS}Enter option\n> "))
            if option not in OPTIONS:
                raise ValueError
        except ValueError as e:
            logging.error(f"{e}, Invalid option {option}. Please try again.")
            option = self.handle_option()
        return option

    def handle_name(self, name: str = "") -> str:
        """
        Check name

        Args:
            name: Name
        Returns:
            str: Name
        """
        name = name or self.config.get("name", "") or input("Enter work name\n> ")
        return self.handle_name(name) if not name else name

    def handle_filepath(self, poscar: str = "", force: bool = True) -> Path:
        """
        Check POSCAR file path

        Args:
            poscar: POSCAR file path
            force: Force to get a POSCAR file
        Returns:
            Path: POSCAR file path
        """
        poscar = poscar or self.config.get("poscar", "")
        if force:
            prompt = "Enter POSCAR file path\n> "
            path = Path(poscar or input(prompt).strip()).absolute()
            if not path.is_file():
                logging.error(f"Invalid file path {str(path)}. Please try again.")
                path = self.handle_filepath(force=force)
        else:
            prompt = (
                "Enter POSCAR file path or leave blank to use structure_info[ ]\n> "
            )
            poscar = poscar or input(prompt).strip()
            path = self.handle_filepath(poscar, force=True) if poscar else Path(poscar)
        return path

    def handle_outdir(self, outdir: str = "") -> Path:
        """
        Check output directory

        Args:
            outdir: Output directory
        Returns:
            Path: Output directory
        """
        outdir = (
            outdir
            or self.config.get("outdir", "")
            or input("Enter output directory[output/]\n>")
        )
        outdir = outdir or "output"
        out = workdir.joinpath(outdir)
        if not out.exists():
            out.mkdir(parents=True, exist_ok=True)
        return out

    def handle_structure(self, phase: str = "", force: bool = True) -> dict[str, dict]:
        """
        Check structure info

        Args:
            phase: Phase of structure
            force: Force to get a phase
        Returns:
            Dict: Structure info of phase
        """
        phase = phase or self.config.get("phase", "")
        if force:
            prompt = "Enter phase to build unitcell (fcc, bcc, hcp, ...)\n> "
            phase = phase or input(prompt)
            while phase not in self.config:
                logging.error(f"Phase {phase} not found. Please try again.")
                phase = input(prompt)
        else:
            prompt = "Enter phase to load sofs data or leave blank to shuffle atoms only[ ]\n> "
            phase = phase or input(prompt).strip()
            if phase not in self.config:
                logging.warning(f"Phase {phase} not found. Shuffle atoms only.")
        return deepcopy(self.config.get(phase, {}))

    def handle_factors(self, factors: list[int] = []) -> tuple[int, int, int]:
        """
        Check supercell factors
        Args:
            factors: Supercell factors
        Returns:
            Tuple: Supercell factors
        """
        prompt = "Enter supercell factors (x y z)[3 3 3]\n> "
        factors = factors or self.config.get("supercell_factors", [])
        try:
            factors = factors or list(map(int, input(prompt).split()[:3])) or [3, 3, 3]
            if len(factors) != 3 or any(f <= 0 for f in factors):
                raise ValueError(f"Needs 3 positive integers, got {factors}.")
            return tuple(factors)  # type: ignore
        except ValueError as e:
            logging.error(f"{e}. Please try again.")
            return self.handle_factors()

    def handle_miller(self, miller: list[int] = []) -> tuple[int, int, int]:
        """
        Check Miller index.
        Args:
            miller: Miller index
        Return:
            Tuple: Miller index
        """
        prompt = "Enter slice direction (x y z)[0 0 1]\n> "
        miller = miller or self.config.get("slie_direction", [])
        try:
            miller = miller or list(map(int, input(prompt).split()[:3])) or [0, 0, 1]
            if len(miller) != 3 or any(not isinstance(m, int) for m in miller):
                raise ValueError(f"Needs 3 integers, got {miller}.")
            return tuple(miller)  # type: ignore
        except ValueError as e:
            logging.error(f"{e}. Please try again.")
            return self.handle_miller()

    def handle_seeds(self, seeds: list[int | None] = [None]) -> list[int | None]:
        """
        Check shuffle seeds.
        Args:
            seeds: Shuffle seeds
        Returns:
            List: Shuffle seeds
        """

        seeds = seeds or self.config.get("shuffle_seeds", [None])
        seeds = [ord(s) if not (isinstance(s, int) or s is None) else s for s in seeds]
        return seeds

    def handle_supercell(
        self,
        poscar: str = "",
        outdir: str = "",
        factors: list[int] = [],
        by_ase: bool = False,
    ) -> Path:
        """
        Handle supercell generation

        Args:
            poscar: POSCAR file path
            outdir: Output directory
            factors: Supercell factors
            by_ase: Whether to use ase.make_supercell()
        Returns:
            Path: Supercell POSCAR file path
        """

        path = self.handle_filepath(poscar, force=False)
        out = self.handle_outdir(outdir)
        if path.is_file():
            logging.info(f"Using POSCAR file {path}")
            path = self.handle_filepath(poscar, force=True)
        else:
            logging.warning(f"POSCAR file {path} not found, using structure_info")
            structure_info = self.handle_structure(force=True)
            path = unitcell2file(structure_info=structure_info, outdir=out)
        factrs = self.handle_factors(factors)
        supercell = supercell2file(
            poscar=path, outdir=out, factors=factrs, by_ase=by_ase
        )
        return supercell

    def handle_modeling(
        self,
        name: str = "",
        poscar: str = "",
        outdir: str = "",
        factors: list[int] = [],
        seeds: list[int | None] = [None],
        batch_size: int = 1,
        enable_sqs: bool = False,
    ) -> list[Path]:
        """
        Handle modeling workflow

        Args:
            name: Name of modeling
            poscar: POSCAR file path
            outdir: Output directory
            factors: Supercell factors
            seeds: Shuffle seeds
            batch_size: Batch size for modeling
            enable_sqs: Whether to enable SQS
        Returns:
            List: List of output files
        """
        name = self.handle_name(name)
        path = self.handle_filepath(poscar, force=False)
        out = self.handle_outdir(outdir)
        factrs = self.handle_factors(factors)
        structure_info = self.handle_structure(force=False)
        seeds = self.handle_seeds(seeds)
        batch_size = batch_size or self.config.get("batch_size", 1)
        enable_sqs = enable_sqs or self.config.get("enable_sqs", False)

        results = run_modeling(
            name=name,
            poscar=path,
            outdir=out,
            supercell_factors=factrs,
            structure_info=structure_info,
            shuffle_seeds=seeds,
            batch_size=batch_size,
            enable_sqs=enable_sqs,
        )
        return results

    def handle_countcn(
        self, name: str = "", poscar: str = "", outdir: str = ""
    ) -> Path:
        """
        Handle CN counting

        Args:
            name: Name of CN counting
            poscar: POSCAR file path
            outdir: Output directory
        Returns:
            Path: Output directory
        """
        name = self.handle_name(name)
        path = self.handle_filepath(poscar)
        out = self.handle_outdir(outdir)
        counter = CNCounter(name, path)
        result = counter.countCN2files_acc(outdir=out)
        return result

    def handle_slice(
        self, name: str = "", poscar: str = "", outdir: str = "", miller: list[int] = []
    ) -> list[Path]:
        """
        Handle slicing

        Args:
            name: Name of slicing
            poscar: POSCAR file path
            outdir: Output directory
            miller: Miller index
        Returns:
            List: List of output files
        """
        name = self.handle_name(name)
        path = self.handle_filepath(poscar)
        out = self.handle_outdir(outdir)
        milr = self.handle_miller(miller)
        slicer = Slicer(name, path, milr)
        results = slicer.slice2files(outdir=out)
        return results

    def handle_slice_to_countcn(
        self, name: str = "", poscar: str = "", outdir: str = "", miller: list[int] = []
    ) -> list[Path]:
        """
        Handle slicing and counting CNs

        Args:
            name: Name of slicing and counting
            poscar: POSCAR file path
            outdir: Output directory
            miller: Miller index
        Returns:
            List: List of output dirs
        """
        name = self.handle_name(name)
        path = self.handle_filepath(poscar)
        out = self.handle_outdir(outdir)
        milr = self.handle_miller(miller)
        results = slice2files_with_countcn(
            name=name, poscar=path, outdir=out, miller_index=milr
        )
        return results

    def run(self, poscar: str = "", option: int = 0):
        print(INFO_CLI)
        poscar = poscar or self.config.get("poscar", "")
        while True:
            try:
                option = self.handle_option(option)
                if option not in self._func_map:
                    continue
                func = self._func_map[option]
                try:
                    result = func(**{"poscar": poscar})
                except KeyboardInterrupt:
                    print("\n")
                except Exception:
                    traceback.print_exc()
                finally:
                    input("Press Enter to continue...")
                    poscar = ""
                    option = 0
            except KeyboardInterrupt:
                break

        print("\n")
        print("Thank you for using PoscarKit")


def main():
    parser = argparse.ArgumentParser(description="PoscarKit")
    parser.add_argument(
        "poscar",
        type=str,
        nargs="?",
        default="",
        help="POSCAR file",
    )
    parser.add_argument(
        "-p",
        "--poscar",
        type=str,
        default="",
        help="POSCAR file",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default="",
        help="Config file",
    )
    parser.add_argument(
        "-o",
        "--option",
        type=int,
        default=0,
        choices=OPTIONS,
        help="Option",
    )
    args = parser.parse_args()
    try:
        poscarkit = PoscarKit(cfg_path=args.config)
        poscarkit.run(poscar=args.poscar, option=args.option)
    except Exception:
        traceback.print_exc()
        input("Press Enter to exit...")
        sys.exit(1)


if __name__ == "__main__":
    main()

# 编译指令
# nuitka --standalone --onefile --output-dir=dist --jobs=4 --lto=yes `
# --enable-plugin=tk-inter --enable-plugin=no-qt --windows-console-mode=disable `
# --windows-icon-from-ico="icon.ico" --onefile-no-compression `
# --enable-plugin=upx --upx-binary="D:\\Programs\\upx-5.0.2-win64\\upx.exe" `
# --nofollow-import-to=matplotlib.tests --nofollow-import-to=pandas.tests `
# --nofollow-import-to=pytest --nofollow-import-to=setuptools.tests `
# --output-filename=poscarkit-0.8.3.exe `
# --file-version=0.8.3 `
# --copyright="(C) 2025 MCMF, Fuzhou University" `
# poscarkit.py
