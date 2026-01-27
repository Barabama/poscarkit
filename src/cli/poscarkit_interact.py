# src/cli/poscarkit_interact.py

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
  4) Count CN        : CN(Coordinate Numbers) and NN(Nearest Neighbors).
  5) Slice to layers : normal to <slice_direction> also known as Miller Index.
  6) Slice to CountCN: Slice to layers and then count CN for each layer.
  7) Make Supercell  : based on <supercell_factors> along the basis vectors.
  8) Compare         : Compare two POSCAR files.
  9) Merge           : Merge two POSCAR files.
 10) Separate        : Separate a POSCAR file by groups.
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
    <cutoff_mult> -------- Multiplier for cutoff distance.
    <by_ase> ------------- Use ASE to CountCN and make supercell.
    <separate_key> ------- Key for separating POSCAR file.

  3) Run Modeling    : Modeling is set to do First [Make Supercell] and then
    shuffle Atoms. Shuffle based on <shuffle_seeds> and Allocate atoms based
    on <supercell_factors> and SOFs data of <phase>. These works will be
    grouped by sublattices, output each sublattice POSCARs and a whole POSCAR.
    If not <phase> provided, Shuffle atoms only.
    
  4) Count CN        : Count CN(Coordinate Numbers) of all the pairs of atoms
    in the supercell and calculate the NN(Nearest Neighbors). A*-B as d1NN
    (the Distance of the First Nearest Neighbors of the same type atom),
    the CN of symbol A coordinated with symbol B.

  5) Slice to layers : Slice atoms to get all the single layers normal to
    <slice_direction> (such as [001], [110], [111]) and plot the layers.
    
  6) Slice to CountCN: First slice atoms to layers, then count coordination
    numbers for each layer separately. This combines functions Slice and CountCN.

  7) Make Supercell  : Choose a POSCAR or the unitcell of prescribed prototype
    built by config, and then expand it to a supercell based on
    <supercell_factors> (such as (3x3x3), (30x30x30) for Cube; (2x2x2) for Hex).

  8) Compare         : Compare two POSCAR files.

  9) Merge           : Merge two POSCAR files.

 10) Separate        : Separate a POSCAR file by groups according to symbol,
   coords, x, y, z, note. 
"""


OPTIONS = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
WORKDIR = Path(sys.argv[0]).absolute().parent
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s[%(levelname)s]%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


class Config(TypedDict, total=False):
    name: str
    poscar: str
    outdir: str
    phase: str
    supercell_factors: list[int]
    slie_direction: list[int]
    shuffle_seeds: list[int | None]
    batch_size: int
    enable_sqs: bool
    cutoff_mult: float
    by_ase: bool
    separate_key: str


class PoscarkitInteract:

    def __init__(self, cfg_path: str = "", cfg: dict = {}):
        self._func_map = {
            1: self.show_help,
            2: self.read_config,
            3: self.run_modeling,
            4: self.run_countcn,
            5: self.run_slice,
            6: self.run_slice_to_countcn,
            7: self.run_supercell,
            8: self.run_compare,
            9: self.run_merge,
            10: self.run_separate,
        }
        self.config = self.read_config(cfg_path, cfg)

    def show_help(self, **kwargs):
        print(INFO_HELP)

    def read_config(self, cfg_path: Path | str = "", cfg: dict = {}, **kwargs) -> Config:
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
        cfg_path = Path(cfg_path) if isinstance(cfg_path, str) else cfg_path
        if cfg_path.is_file():
            with open(cfg_path, "r", encoding="utf-8") as tf:
                config = Config(**tomllib.loads(tf.read()))
                logging.info(f"Read config from {str(cfg_path)}")
        else:
            logging.warning(f"Config file {str(cfg_path)} does not exist, using default one.")
            cfg_path = WORKDIR.joinpath("config.toml")
            with open(cfg_path, "w", encoding="utf-8") as tf:
                tf.write(DEFAULT_CONFIG)
            config = self.read_config(cfg_path=cfg_path, **kwargs)
        return config

    def _handle_option(self, option: int = 0) -> int:
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
            option = self._handle_option()
        return option

    def _handle_name(self, name: str = "") -> str:
        """
        Check name

        Args:
            name: Name
        Returns:
            str: Name
        """
        name = name or self.config.get("name", "") or input("Enter work name\n> ")
        return self._handle_name(name) if not name else name

    def _handle_poscar(self, poscar: str = "", force: bool = True) -> Path:
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
                path = self._handle_poscar(force=force)
        else:
            prompt = "Enter POSCAR file path or leave blank to use structure_info[ ]\n> "
            poscar = poscar or input(prompt).strip()
            path = self._handle_poscar(poscar, force=True) if poscar else Path(poscar)
        return path

    def _handle_outdir(self, outdir: str = "") -> Path:
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
        out = WORKDIR.joinpath(outdir)
        if not out.exists():
            out.mkdir(parents=True, exist_ok=True)
        return out

    def _handle_structure(self, phase: str = "", force: bool = True) -> dict[str, dict]:
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

    def _handle_factors(self, factors: list[int] = []) -> tuple[int, int, int]:
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
            return self._handle_factors()

    def _handle_miller(self, miller: list[int] = []) -> tuple[int, int, int]:
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
            return self._handle_miller()

    def _handle_seeds(self, seeds: list[int | None] = [None]) -> list[int | None]:
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

    def run_modeling(self) -> list[Path]:
        """
        Run modeling workflow

        Returns:
            List: List of output files
        """
        name = self._handle_name()
        poscar = self._handle_poscar(force=False)
        outdir = self._handle_outdir()
        factrs = self._handle_factors()
        info = self._handle_structure(force=False)
        seeds = self._handle_seeds()
        batch_size = self.config.get("batch_size", 1)
        enable_sqs = self.config.get("enable_sqs", False)
        interations = self.config.get("iterations", 1e7)

        results = run_modeling(
            name=name,
            poscar=poscar,
            outdir=outdir,
            supercell_factors=factrs,
            structure_info=info,
            shuffle_seeds=seeds,
            batch_size=batch_size,
            enable_sqs=enable_sqs,
            iterations=interations,
        )
        return results

    def run_countcn(self) -> Path:
        """
        Run CN counting

        Returns:
            Path: Output directory
        """
        name = self._handle_name()
        posc = self._handle_poscar()
        out = self._handle_outdir()
        mult = self.config.get("cutoff_mult", 1.1)
        parallel = self.config.get("parallel", 2)
        by_ase = self.config.get("by_ase", False)

        counter = CNCounter(name, posc)
        result = counter.countCN2files(
            outdir=out,
            cutoff_mult=mult,
            parallel=parallel,
            by_ase=by_ase,
        )
        return result

    def run_slice(self) -> list[Path]:
        """
        Run slicing

        Args:

        Returns:
            List: List of output files
        """
        name = self._handle_name()
        poscar = self._handle_poscar()
        outdir = self._handle_outdir()
        miller = self._handle_miller()
        slicer = Slicer(name, poscar, miller)
        results = slicer.slice2files(outdir=outdir)
        return results

    def run_slice_to_countcn(self) -> list[Path]:
        """
        Run slicing and counting CNs

        Returns:
            List: List of output dirs
        """
        name = self._handle_name()
        poscar = self._handle_poscar()
        outdir = self._handle_outdir()
        miller = self._handle_miller()
        results = slice2files_with_countcn(
            name=name,
            poscar=poscar,
            outdir=outdir,
            miller_index=miller,
        )
        return results

    def run_supercell(
        self,
    ) -> Path:
        """
        Run supercell generation

        Args:

        Returns:
            Path: Supercell POSCAR file path
        """

        poscar = self._handle_poscar(force=False)
        outdir = self._handle_outdir()
        by_ase = self.config.get("by_ase", False)
        if poscar.is_file():
            logging.info(f"Using POSCAR file {poscar}")
            poscar = self._handle_poscar(force=True)
        else:
            logging.warning(f"POSCAR file {poscar} not found, using structure_info")
            structure_info = self._handle_structure(force=True)
            poscar = unitcell2file(structure_info=structure_info, outdir=outdir)
        factrs = self._handle_factors(factors)
        supercell = supercell2file(poscar=poscar, outdir=outdir, factors=factrs, by_ase=by_ase)
        return supercell

    def run_compare(self):
        """Run comparing"""
        poscar1 = self._handle_poscar(force=True)
        poscar2 = self._handle_poscar(force=True)
        SimplePoscar.compare_poscar(poscar1, poscar2)

    def run_merge(self) -> Path:
        """Run merging

        Returns:
            Path: Merged POSCAR file path
        """
        poscar1 = self._handle_poscar(force=True)
        poscar2 = self._handle_poscar(force=True)
        outdir = self._handle_outdir()
        merged = SimplePoscar.merge_poscar(poscar1, poscar2, outdir)
        return merged

    def run_separate(self) -> list[Path]:
        """Run separating

        Returns:
            List: List of output files
        """
        poscar = self._handle_poscar(force=True)
        outdir = self._handle_outdir()
        key = self.config.get("separate_key", "note")
        outputs = SimplePoscar.separate_poscar(poscar, outdir, key)
        return outputs

    def run(self):
        print(INFO_CLI)
        while True:
            try:
                option = self._handle_option()
                if option not in self._func_map:
                    continue
                func = self._func_map[option]
                try:
                    result = func()
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
    try:
        poscarkit = PoscarkitInteract()
        poscarkit.run()
    except Exception:
        traceback.print_exc()
        input("Press Enter to exit...")
        sys.exit(1)


if __name__ == "__main__":
    main()
