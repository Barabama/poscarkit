import sys
import argparse
import logging
import traceback
from pathlib import Path
from typing import Any


from src.config import VERSION, CONTACT, DEVELOPER, DEFAULT_CONFIG
from src.cli.poscarkit import (
    cmd_help,
    cmd_modeling,
    cmd_supercell,
    cmd_compare,
    cmd_merge,
    cmd_separate,
    cmd_countcn,
    cmd_slice,
    cmd_slice_to_countcn,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")

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
    numbers for each layer separately. This combines functions 7 and 8.

  7) Make Supercell  : Choose a POSCAR or the unitcell of prescribed prototype
    built by config, and then expand it to a supercell based on
    <supercell_factors> (such as (3x3x3), (30x30x30) for Cube; (2x2x2) for Hex).

  8) Compare         : Compare two POSCAR files.

  9) Merge           : Merge two POSCAR files.

 10) Separate        : Separate a POSCAR file by groups.

"""

OPTIONS = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)


class PoscarkitInteract:
    """Poscarkit interactive mode"""

    def __init__(self):
        self._func_map = {
            1: self.show_help,
            2: self.read_config,
            3: self.handle_modeling,
            4: self.handle_countcn,
            5: self.handle_slice,
            6: self.handle_slice_to_countcn,
            7: self.handle_supercell,
            8: self.handle_compare,
            9: self.handle_merge,
            10: self.handle_separate,
        }
        self.config = self.read_config()

    def show_help(self, **kwargs):
        """Show help information"""
        print(INFO_HELP)

    def read_config(self, cfg_path: str = "", cfg: dict = {}, **kwargs):
        """
        Read config file
        Args:
            cfg_path: config.toml path
            cfg: Dict to parse
        Returns:
            dict: config
        """
        import tomllib
        from pathlib import Path
        
        if cfg:
            self.config = cfg
            return self.config
        
        # Get work directory
        workdir = Path(sys.argv[0]).parent
        path = workdir.joinpath("config.toml")
        
        if not path.is_file():
            logging.warning(
                f"Config file {str(path)} does not exist, using default one."
            )
            with open(path, "w", encoding="utf-8") as tf:
                tf.write(DEFAULT_CONFIG)
            self.config = self.read_config(cfg_path=str(path), **kwargs)
        elif path.is_file():
            with open(path, "r", encoding="utf-8") as tf:
                self.config = tomllib.loads(tf.read())
                logging.info(f"Read config from {str(path)}")
        
        return self.config

    def handle_option(self, option: int = 0) -> int:
        """Handle menu option selection"""
        try:
            option = option or int(input(f"{INFO_OPTIONS}Enter option\n> "))
            if option not in OPTIONS:
                raise ValueError
        except ValueError as e:
            logging.error(f"{e}, Invalid option {option}. Please try again.")
            option = self.handle_option()
        return option

    def handle_name(self, name: str = "") -> str:
        """Handle name input"""
        try:
            name = name or input("Enter work name\n> ")
            return self.handle_name() if not name else name
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            raise

    def handle_filepath(self, poscar: str = "", force: bool = True) -> str:
        """Handle filepath input"""
        try:
            if force:
                prompt = "Enter POSCAR file path\n> "
                path = poscar or input(prompt).strip()
                if not path:
                    print("Invalid file path. Please try again.")
                    path = self.handle_filepath(force=force)
            else:
                prompt = "Enter POSCAR file path or leave blank to use structure_info\n> "
                path = poscar or input(prompt).strip()
            return path
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            raise

    def handle_outdir(self, outdir: str = "") -> str:
        """Handle output directory input"""
        try:
            outdir = outdir or input("Enter output directory[output/]\n>")
            outdir = outdir or "output"
            return outdir
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            raise

    def handle_factors(self, factors: list = []) -> tuple:
        """Handle supercell factors input"""
        try:
            prompt = "Enter supercell factors (x y z)[3 3 3]\n> "
            factors = factors or input(prompt).strip()
            factors = factors or "3 3 3"
            try:
                factors = list(map(int, factors.split()[:3]))
                if len(factors) != 3 or any(f <= 0 for f in factors):
                    raise ValueError
                return tuple(factors)
            except ValueError:
                print("Invalid factors. Please try again.")
                return self.handle_factors()
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            raise

    def handle_miller(self, miller: list = []) -> tuple:
        """Handle Miller index input"""
        try:
            prompt = "Enter slice direction (x y z)[0 0 1]\n> "
            miller = miller or input(prompt).strip()
            miller = miller or "0 0 1"
            try:
                miller = list(map(int, miller.split()[:3]))
                if len(miller) != 3:
                    raise ValueError
                return tuple(miller)
            except ValueError:
                print("Invalid Miller index. Please try again.")
                return self.handle_miller()
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            raise

    def handle_modeling(self, **kwargs):
        """Handle modeling menu"""
        try:
            name = self.handle_name()
            poscar = self.handle_filepath(kwargs.get("poscar", ""), force=False)
            outdir = self.handle_outdir()
            factors = self.handle_factors()

            # Create namespace for cmd_modeling
            class Namespace:
                def __init__(self):
                    self.name = name
                    self.poscar = poscar
                    self.factors = factors
                    self.outdir = outdir
                    self.config = None
                    self.phase = None
                    self.seeds = None
                    self.batch_size = 1
                    self.enable_sqs = False
                    self.iterations = 1e7

            args = Namespace()
            cmd_modeling(args)
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            return
        except Exception as e:
            print(f"Error: {e}")
            traceback.print_exc()

    def handle_countcn(self, **kwargs):
        """Handle countcn menu"""
        try:
            name = self.handle_name()
            poscar = self.handle_filepath(kwargs.get("poscar", ""))
            outdir = self.handle_outdir()

            # Create namespace for cmd_countcn
            class Namespace:
                def __init__(self):
                    self.name = name
                    self.poscar = poscar
                    self.outdir = outdir
                    self.cutoff_mult = 1.1
                    self.parallel = 2
                    self.by_ase = False

            args = Namespace()
            cmd_countcn(args)
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            return
        except Exception as e:
            print(f"Error: {e}")
            traceback.print_exc()

    def handle_slice(self, **kwargs):
        """Handle slice menu"""
        try:
            name = self.handle_name()
            poscar = self.handle_filepath(kwargs.get("poscar", ""))
            outdir = self.handle_outdir()
            miller = self.handle_miller()

            # Create namespace for cmd_slice
            class Namespace:
                def __init__(self):
                    self.name = name
                    self.poscar = poscar
                    self.outdir = outdir
                    self.miller_index = miller

            args = Namespace()
            cmd_slice(args)
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            return
        except Exception as e:
            print(f"Error: {e}")
            traceback.print_exc()

    def handle_slice_to_countcn(self, **kwargs):
        """Handle slice-to-countcn menu"""
        try:
            name = self.handle_name()
            poscar = self.handle_filepath(kwargs.get("poscar", ""))
            outdir = self.handle_outdir()
            miller = self.handle_miller()

            # Create namespace for cmd_slice_to_countcn
            class Namespace:
                def __init__(self):
                    self.name = name
                    self.poscar = poscar
                    self.outdir = outdir
                    self.miller_index = miller

            args = Namespace()
            cmd_slice_to_countcn(args)
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            return
        except Exception as e:
            print(f"Error: {e}")
            traceback.print_exc()

    def handle_supercell(self, **kwargs):
        """Handle supercell menu"""
        try:
            poscar = self.handle_filepath(kwargs.get("poscar", ""))
            outdir = self.handle_outdir()
            factors = self.handle_factors()

            # Create namespace for cmd_supercell
            class Namespace:
                def __init__(self):
                    self.poscar = poscar
                    self.outdir = outdir
                    self.factors = factors
                    self.by_ase = False

            args = Namespace()
            cmd_supercell(args)
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            return
        except Exception as e:
            print(f"Error: {e}")
            traceback.print_exc()

    def handle_compare(self, **kwargs):
        """Handle compare menu"""
        try:
            poscar1 = self.handle_filepath(kwargs.get("poscar", ""))
            poscar2 = self.handle_filepath("", force=True)

            # Create namespace for cmd_compare
            class Namespace:
                def __init__(self):
                    self.poscar1 = poscar1
                    self.poscar2 = poscar2

            args = Namespace()
            cmd_compare(args)
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            return
        except Exception as e:
            print(f"Error: {e}")
            traceback.print_exc()

    def handle_merge(self, **kwargs):
        """Handle merge menu"""
        try:
            poscar1 = self.handle_filepath(kwargs.get("poscar", ""))
            poscar2 = self.handle_filepath("", force=True)
            outdir = self.handle_outdir()

            # Create namespace for cmd_merge
            class Namespace:
                def __init__(self):
                    self.poscar1 = poscar1
                    self.poscar2 = poscar2
                    self.outdir = outdir

            args = Namespace()
            cmd_merge(args)
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            return
        except Exception as e:
            print(f"Error: {e}")
            traceback.print_exc()

    def handle_separate(self, **kwargs):
        """Handle separate menu"""
        try:
            poscar = self.handle_filepath(kwargs.get("poscar", ""))
            outdir = self.handle_outdir()
            key = input("Enter key to separate by [note]\n>") or "note"

            # Create namespace for cmd_separate
            class Namespace:
                def __init__(self):
                    self.poscar = poscar
                    self.outdir = outdir
                    self.key = key

            args = Namespace()
            cmd_separate(args)
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            return
        except Exception as e:
            print(f"Error: {e}")
            traceback.print_exc()

    def run(self, poscar: str = "", option: int = 0):
        """Run interactive mode"""
        print(INFO_CLI)
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
    """Main function"""
    try:
        app = PoscarkitInteract()
        app.run()
    except Exception:
        traceback.print_exc()
        input("Press Enter to exit...")
        sys.exit(1)


if __name__ == "__main__":
    main()
