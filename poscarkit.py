"""main.py"""

import argparse
import os
import sys
import tomllib
from copy import deepcopy
from dataclasses import dataclass, field

sys.path.append(os.path.join(os.path.dirname(__file__), "PoscarTools"))
from Supercell import supercell2file
from AtomSlice import slice2file
from AtomShuffle import shuffle2files
from AtomAllocate import allocate2file
from AtomCountCN import countCN2files


VERSION = "0.7"
INFO_EXEC = f"""
--- POSCAR tool (ver{VERSION}) ---
This tool has many uses of supercell, slice, shuffle, allocation, etc.
developed by FZU-MCMF, main developers: Gao Min-Liang, Wu Bo*, et al.,
please contact via wubo@fzu.edu.cn, 654489521@qq.com in case.
Ctrl+C to exit.
----------------------------
"""

INFO_CHOICES = """
=======================================
  1. Read config
  2. Supercell   3. Slice    4. Shuffle
  5. Allocation  6. FlowWork 7. CountCN
  Exit. Ctrl+C
=======================================
"""
CHOICES = (1, 2, 3, 4, 5, 6 ,7)


@dataclass
class Config:
    Filepath: str = field(default="")

    # Used for 'Supercell', 'Allocation', 'ShuffleAllocation'
    SupercellFactors: list[int] = field(default_factory=list)

    # Used for 'Shuffle', 'Allocation', 'ShuffleAllocation'
    Structure: str = field(default="")

    # Used for 'Shuffle'
    ShuffleSeeds: list[int | None] = field(default_factory=lambda: [None])

    # Used for 'Allocation', 'ShuffleAllocation'
    Shuffle: bool = field(default=False)

    # Direction of slicing Used for 'Slice'
    SliceDirection: list[int] = field(default_factory=list)

    FCC: dict[str, dict[str, list | dict]] = field(default_factory=dict)
    BCC: dict[str, dict[str, list | dict]] = field(default_factory=dict)
    HCP: dict[str, dict[str, list | dict]] = field(default_factory=dict)


def read_config() -> Config:
    with open(os.path.join(os.path.dirname(sys.argv[0]), "config.toml"), "rb") as tf:
        config = Config(**tomllib.load(tf))
    return config


def handle_supercell_factors(factors: tuple[int, int, int] = ()) -> tuple[int, int, int]:
    """Check supercell factors."""
    while True:
        try:
            prompt = "Enter Supercell factors (x y z, default 3 3 3) >>> "
            factors = factors or tuple(map(int, input(prompt).split()[:3])) or (3, 3, 3)

            if len(factors) != 3 or any(f <= 0 for f in factors):
                raise ValueError("Input must be 3 positive integers")
            return factors
        except ValueError as e:
            print(f"Invalid input: {e}. Please try again.")
            factors = ()


def handle_structure(config: Config, structure: str = "") -> dict[str, list | dict]:
    """return a structure (FCC, BCC, HCP) information from the configuration."""
    while True:
        structure = structure or input("Enter Structure (fcc, bcc, hcp) >>> ").upper()
        if structure in config.__dict__:
            return deepcopy(config.__dict__[structure])
        else:
            print(f"{structure} not found in configuration. Please try again.")


def handle_slice_direction(direction: tuple[int, int, int] = ()) -> tuple[int, int, int]:
    """Check slice direction."""
    while True:
        try:
            prompt = "Enter slice direction (x y z, default=0 0 1) >>> "
            direction = direction or tuple(map(int, input(prompt).split()[:3])) or (0, 0, 1)

            if len(direction) != 3 or any(not isinstance(i, int) for i in direction):
                raise ValueError("Input must be a tuple of three integers.")
            return direction
        except ValueError as e:
            print(f"Invalid input: {e}. Please try again.")
            direction = ()


def main(config: Config, filepath: str = "", option: int = 0):
    """Main function."""
    print(INFO_EXEC)
    while True:
        try:
            # Check option and filepath, passed option==1
            option = option or int(input(f"{INFO_CHOICES}Enter choice >>> "))
            if option != 1:
                filepath = filepath or config.Filepath or input("Enter filepath >>> ")
                if not os.path.isfile(filepath):
                    print("[ERROR]File not found! Please try again.")
                    filepath = ""
                    config.Filepath = ""
                    option = 0
                    continue

            # Execute option
            # Flag = False
            match option:
                case 1:
                    print("Reading config...", end="")
                    config = read_config()
                    filepath = ""
                    option = 0
                    print("Done.")

                case 2:
                    factors = tuple(config.SupercellFactors)
                    if factors := handle_supercell_factors(factors):
                        supercell2file(filepath, factors)
                    # Flag = True
                case 3:
                    direction = tuple(config.SliceDirection)
                    if direction := handle_slice_direction(direction):
                        slice2file(filepath, direction)
                    # Flag = True

                case 4:
                    structure = handle_structure(config, config.Structure.upper())
                    seeds = config.ShuffleSeeds or [None]
                    print
                    if structure and seeds:
                        fs = [f for f in shuffle2files(filepath, structure, seeds)]
                    # Flag = True

                case 5:
                    structure = handle_structure(config, config.Structure.upper())
                    factors = handle_supercell_factors(tuple(config.SupercellFactors))
                    shuffle = config.Shuffle
                    if structure and factors:
                        allocate2file(filepath, structure, factors, shuffle)
                    # Flag = True

                case 6:
                    factors = handle_supercell_factors(tuple(config.SupercellFactors))
                    structure = handle_structure(config, config.Structure.upper())
                    seeds = config.ShuffleSeeds or [None]
                    shuffle = config.Shuffle
                    if structure and factors:
                        supercell_file = supercell2file(filepath, factors)
                        shuffled_files = [f for f in shuffle2files(supercell_file, structure, seeds)]
                        for f in shuffled_files:
                            _ = allocate2file(f, structure, factors, shuffle)
                    # Flag = True

                case 7:
                    factors = handle_supercell_factors(tuple(config.SupercellFactors))
                    if factors:
                        countCN2files(filepath, factors)

                
                case _:
                    raise ValueError(f"Invalid option {option}")

            # # Reset
            # if Flag:
            filepath = ""
            option = 0

        except ValueError as e:
            print(f"[ERROR]Invalid input: {e}", end=" ")
            filepath = ""
            option = 0
        except KeyboardInterrupt:
            print("\nExiting... Thank you for using POSCAR tool!")
            sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Poscar Tool")
    parser.add_argument("filepath", type=str, nargs="?", default="", help="Poscar file path")
    parser.add_argument("-f", "--file", type=str, default="", help="Poscar file path")
    parser.add_argument("-c", "--choice", type=int, default=0, choices=CHOICES,
                        help="Choice of operation")
    args = parser.parse_args()
    file = args.filepath or args.file
    opt = args.choice
    cfg = read_config()
    main(config=cfg, filepath=file, option=opt)
