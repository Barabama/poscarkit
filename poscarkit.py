# poscarkit.py

import argparse
import glob
import logging
import os
import sys
import tomllib
import traceback
from copy import deepcopy
from dataclasses import dataclass, field

from PoscarTools.Supercell import supercell2file
from PoscarTools.AtomSlice import slice2file
from PoscarTools.AtomShuffle import shuffle2files
from PoscarTools.AtomAllocate import allocate2file
from PoscarTools.AtomCountCN import countCN2files


VERSION = "0.7.5"
INFO_EXEC = f"""
--- POSCAR tool (ver{VERSION}) ---
This tool has many uses of supercell, slice, shuffle, allocation, etc.
developed by FZU-MCMF, main developers: Gao Min-Liang, Wu Bo*, et al.,
please contact via wubo@fzu.edu.cn, 654489521@qq.com in case.
Ctrl+C to exit.
------------------------------
"""
INFO_CHOICES = """
=======================================
  1. Read config
  2. Supercell   3. Slice    4. Shuffle
  5. Allocation  6. WorkFlow 7. CountCN
  Exit. Ctrl+C
=======================================
"""
CHOICES = (2, 3, 4, 5, 6, 7)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


@dataclass
class Config:
    Filepath: str = field(default="")

    # Used for 'Supercell', 'Allocation', 'ShuffleAllocation'
    SupercellFactors: list[int] = field(default_factory=lambda: [])

    # Used for 'Shuffle', 'Allocation', 'ShuffleAllocation'
    Structure: str = field(default="")

    # Used for 'Shuffle'
    ShuffleSeeds: list[int | None] = field(default_factory=lambda: [None])

    # Used for 'Allocation', 'ShuffleAllocation'
    Shuffle: bool = field(default=False)

    # Direction of slicing Used for 'Slice'
    SliceDirection: list[int] = field(default_factory=lambda: [])

    FCC: dict[str, dict[str, list | dict]] = field(default_factory=dict)
    BCC: dict[str, dict[str, list | dict]] = field(default_factory=dict)
    HCP: dict[str, dict[str, list | dict]] = field(default_factory=dict)


class PoscarKit:

    def __init__(self):
        self.config = self.read_config()
        self._params_map = {
            # "filepath": self._handle_filepath,
            "factors": self._handle_supercell_factors,
            "miller_index": self._handle_miller_index,
            "structure": self._handle_structure,
            "seeds": self._handle_seeds,
            "shuffle": self._handle_shuffle,
        }
        self._func_map = {
            2: {"func": supercell2file, "args": ["factors"]},
            3: {"func": slice2file, "args": ["miller_index"]},
            4: {"func": shuffle2files, "args": ["structure", "seeds"]},
            5: {"func": allocate2file, "args": ["structure", "factors", "shuffle"]},
            6: {"func": self._workflow, "args": ["factors", "structure", "seeds", "shuffle"]},
            7: {"func": countCN2files, "args": []},
        }

    def read_config(self) -> Config:
        cfg_path = os.path.join(os.path.dirname(sys.argv[0]), "config.toml")
        if not os.path.isfile(cfg_path):
            logging.warning("config.toml not found! Using default configuration.")
            self.config = Config()
        else:
            with open(cfg_path, "rb") as tf:
                self.config = Config(**tomllib.load(tf))
        return self.config

    def _handle_option(self, option: int = 0) -> int:
        """Check option, passed option==1."""
        option = option or int(input(f"{INFO_CHOICES}Enter choice >>> "))
        if option == 1:
            self.read_config()
            logging.info("Configurations reloaded.")
            return 0
        return option

    def _handle_filepath(self, filepath: str = "") -> list[str]:
        """Check filepath."""
        files = []
        filepath = filepath or input("Enter filepath >>> ")
        filepath = os.path.abspath(filepath)
        if os.path.isfile(filepath):
            files.append(filepath)
        elif os.path.isdir(filepath):
            patterns = ["*.vasp", "*POSCAR*", "*poscar*"]
            for pattern in patterns:
                files.extend(glob.glob(os.path.join(filepath, "**", pattern), recursive=True))
        else:
            logging.warning(f"No (*.vasp, *POSCAR*, *poscar*) files found in {filepath}")
        return files

    def _handle_supercell_factors(self) -> tuple[int, int, int]:
        """Check supercell factors."""
        factors = tuple(self.config.SupercellFactors)
        while True:
            try:
                prompt = "Enter Supercell factors (x y z, default=3 3 3) >>> "
                factors = factors or tuple(map(int, input(prompt).split()[:3])) or (3, 3, 3)
                if len(factors) != 3 or any(f <= 0 for f in factors):
                    raise ValueError("Input must be 3 positive integers")
                return factors
            except ValueError as e:
                logging.warning(f"{e}. Please try again.")
                factors = ()

    def _handle_miller_index(self) -> tuple[int, int, int]:
        """Check slice direction."""
        miller_index = tuple(self.config.SliceDirection)
        while True:
            try:
                prompt = "Enter slice direction (x y z, default=0 0 1) >>> "
                miller_index = miller_index or tuple(map(int, input(prompt).split()[:3])) or (0, 0, 1)

                if len(miller_index) != 3 or any(not isinstance(i, int) for i in miller_index):
                    raise ValueError("Input must be a tuple of three integers.")
                return miller_index
            except ValueError as e:
                logging.warning(f"{e}. Please try again.")
                miller_index = ()

    def _handle_structure(self) -> dict[str, dict]:
        """return a structure (FCC, BCC, HCP) information from the configuration."""
        structure = self.config.Structure.upper()
        while True:
            structure = structure or input("Enter Structure (fcc, bcc, hcp) >>> ").upper()
            if structure in self.config.__dict__:
                return deepcopy(self.config.__dict__[structure])
            else:
                logging.warning(f"{structure} not found. Please try again.")

    def _handle_seeds(self) -> list[int | None]:
        return self.config.ShuffleSeeds

    def _handle_shuffle(self) -> bool:
        return self.config.Shuffle

    def _workflow(self, filepath: str, factors: tuple[int, int, int],
                  structure: dict[str, dict], seeds: list[int | None], shuffle: bool):
        supercell = supercell2file(filepath, factors)
        shuffleds = shuffle2files(supercell, structure, seeds)
        fs = [allocate2file(f, structure, factors, shuffle) for f in shuffleds]
        return fs

    def run(self, filepath: str = "", option: int = 0):
        """Main function."""
        print(INFO_EXEC)
        filepath = filepath or self.config.Filepath
        while True:
            try:
                option = self._handle_option(option)
                if option not in CHOICES:
                    continue

                func: function = self._func_map[option]["func"]
                args: list[str] = self._func_map[option]["args"]
                kwargs = {k: v() for k, v in self._params_map.items() if k in args}
                filepath = filepath or self.config.Filepath
                outputs = [func(file, **kwargs) for file in self._handle_filepath(filepath)]

                # Reset
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
