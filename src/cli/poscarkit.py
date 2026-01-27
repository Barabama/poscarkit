# src/cli/poscarkit.py

import sys
import argparse
import logging
import tomllib
from pathlib import Path
from typing import Any

from src.modeling.base import SimplePoscar
from src.modeling.countcn import CNCounter
from src.modeling.slice import Slicer
from src.workflow.modeling import run_modeling
from src.modeling.supercell import unitcell2file, supercell2file
from src.workflow.slice_to_countcn import slice2files_with_countcn
from src.config import VERSION, CONTACT, DEVELOPER


WORKDIR = Path(sys.argv[0]).absolute().parent
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s[%(levelname)s]%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def cmd_help(args: argparse.Namespace) -> int:
    print(
        f"""
================================== POSCARKIT ==================================
{DEVELOPER}
{VERSION}
{CONTACT}

A tool for modeling structure POSCAR files, based on Sublattice Occupied 
Fractions (SOFs).

COMMANDS:
  help          Show this help message
  modeling      Run modeling workflow (make supercell and allocate atoms)
  countcn       Calculate coordination numbers and save to files
  slice         Slice a structure by Miller index
  slice-to-countcn Slice a structure and count CNs in each layer
  supercell     Generate supercell from POSCAR file or structure info
  compare       Compare two POSCAR files
  merge         Merge two POSCAR files
  separate      Separate a POSCAR file by groups

EXAMPLES:
  poscarkit help
  poscarkit modeling --poscar POSCAR.vasp --factors 3 3 3 --name my_model
  poscarkit countcn --poscar POSCAR.vasp --name my_cn --outdir output/
  poscarkit slice --poscar POSCAR.vasp --miller-index 1 1 1 --outdir output/
  poscarkit slice-to-countcn --poscar POSCAR.vasp --miller-index 1 1 1 --outdir output/
  poscarkit supercell --poscar POSCAR.vasp --factors 2 2 2 --outdir output/
  poscarkit compare --poscar1 POSCAR1.vasp --poscar2 POSCAR2.vasp
  poscarkit merge --poscar1 POSCAR1.vasp --poscar2 POSCAR2.vasp --outdir output/
  poscarkit separate --poscar POSCAR.vasp --key note --outdir output/
"""
    )
    return 0


def cmd_modeling(args: argparse.Namespace) -> int:
    name = args.name or "modeling"
    poscar = Path(args.poscar) if args.poscar else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    factors = tuple(args.factors)
    phase = args.phase
    config = args.config
    seeds = args.seeds if args.seeds else [None]
    batch_size = args.batch_size
    enable_sqs = args.enable_sqs
    iterations = args.iterations

    structure_info: dict[str, Any] = {}

    if config:
        with open(config, "rb") as f:
            cfg = tomllib.load(f)
    else:
        cfg = {}
    structure_info = cfg.get(phase, {})

    outdir.mkdir(parents=True, exist_ok=True)

    if poscar and not poscar.is_file():
        logging.error(f"POSCAR file not found: {poscar}")
        return 1

    if not poscar and not structure_info:
        logging.error("Either --poscar or (--config and --phase) must be provided")
        return 1

    if poscar:
        poscar_path = poscar
    elif structure_info:
        poscar_path = unitcell2file(structure_info=structure_info, outdir=outdir)
    else:
        poscar_path = Path("")

    results = run_modeling(
        name=name,
        poscar=poscar_path,
        outdir=outdir,
        supercell_factors=factors,
        structure_info=structure_info,
        shuffle_seeds=seeds,
        batch_size=batch_size,
        enable_sqs=enable_sqs,
        iterations=iterations,
    )

    logging.info(f"Modeling completed. Generated {len(results)} structures:")
    for result in results:
        logging.info(f"  - {result}")
    return 0


def cmd_supercell(args: argparse.Namespace) -> int:
    poscar = Path(args.poscar) if args.poscar else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    factors = tuple(args.factors)
    by_ase = args.by_ase

    if not poscar or not poscar.is_file():
        logging.error(f"POSCAR file not found: {poscar}")
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    result = supercell2file(poscar=poscar, outdir=outdir, factors=factors, by_ase=by_ase)

    logging.info(f"Supercell generated: {result}")
    return 0


def cmd_compare(args: argparse.Namespace) -> int:
    poscar1 = Path(args.poscar1) if args.poscar1 else None
    poscar2 = Path(args.poscar2) if args.poscar2 else None

    if not poscar1 or not poscar1.is_file():
        logging.error(f"POSCAR file not found: {poscar1}")
        return 1

    if not poscar2 or not poscar2.is_file():
        logging.error(f"POSCAR file not found: {poscar2}")
        return 1

    SimplePoscar.compare_poscar(poscar1, poscar2)
    return 0


def cmd_merge(args: argparse.Namespace) -> int:
    poscar1 = Path(args.poscar1) if args.poscar1 else None
    poscar2 = Path(args.poscar2) if args.poscar2 else None
    outdir = Path(args.outdir) if args.outdir else Path("output")

    if not poscar1 or not poscar1.is_file():
        logging.error(f"POSCAR file not found: {poscar1}")
        return 1

    if not poscar2 or not poscar2.is_file():
        logging.error(f"POSCAR file not found: {poscar2}")
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    merged_file = SimplePoscar.merge_poscar(poscar1, poscar2, outdir)
    logging.info(f"Merged POSCAR files saved to {merged_file}")
    return 0


def cmd_separate(args: argparse.Namespace) -> int:
    poscar = Path(args.poscar) if args.poscar else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    key = args.key if args.key else "note"

    if not poscar or not poscar.is_file():
        logging.error(f"POSCAR file not found: {poscar}")
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    separated_files = SimplePoscar.separate_poscar(poscar, outdir, key)
    logging.info(f"Separated POSCAR files saved to {outdir}")
    for file in separated_files:
        logging.info(f"- {file}")
    return 0


def cmd_countcn(args: argparse.Namespace) -> int:
    name = args.name or "countcn"
    poscar = Path(args.poscar) if args.poscar else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    cutoff_mult = args.cutoff_mult
    parallel = args.parallel
    by_ase = args.by_ase

    if not poscar or not poscar.is_file():
        logging.error(f"POSCAR file not found: {poscar}")
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    counter = CNCounter(name=name, poscar=poscar)
    result = counter.countCN2files(
        outdir=outdir,
        cutoff_mult=cutoff_mult,
        parallel=parallel,
        by_ase=by_ase,
    )

    logging.info(f"Coordination number analysis completed. Results saved to: {result}")
    return 0


def cmd_slice(args: argparse.Namespace) -> int:
    name = args.name or "slice"
    poscar = Path(args.poscar) if args.poscar else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    miller_index = tuple(args.miller_index)

    if not poscar or not poscar.is_file():
        logging.error(f"POSCAR file not found: {poscar}")
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    slicer = Slicer(name=name, poscar=poscar, miller_index=miller_index)
    results = slicer.slice2files(outdir=outdir)

    logging.info(f"Slicing completed. Generated {len(results)} layers:")
    for result in results:
        logging.info(f"  - {result}")
    return 0


def cmd_slice_to_countcn(args: argparse.Namespace) -> int:
    name = args.name or "slice-to-countcn"
    poscar = Path(args.poscar) if args.poscar else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    miller_index = tuple(args.miller_index)

    if not poscar or not poscar.is_file():
        logging.error(f"POSCAR file not found: {poscar}")
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    results = slice2files_with_countcn(
        name=name, poscar=poscar, outdir=outdir, miller_index=miller_index
    )

    logging.info(f"Slice to CN count completed. Results saved to:")
    for result in results:
        logging.info(f"  - {result}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="poscarkit",
        description="A tool for modeling structure POSCAR files based on SOFs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Help command
    parser_help = subparsers.add_parser("help", help="Show help message")
    parser_help.set_defaults(func=cmd_help)

    # Modeling command
    parser_modeling = subparsers.add_parser("modeling", help="Run modeling workflow")
    parser_modeling.add_argument(
        "--name",
        "-n",
        type=str,
        default="modeling",
        help="Name of the modeling task (default: modeling)",
    )
    parser_modeling.add_argument(
        "--poscar",
        "-p",
        type=str,
        help="Path to POSCAR file",
    )
    parser_modeling.add_argument(
        "--factors",
        "-f",
        type=int,
        nargs=3,
        default=(3, 3, 3),
        metavar=("X", "Y", "Z"),
        help="Supercell factors (e.g., 3 3 3)",
    )
    parser_modeling.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser_modeling.add_argument(
        "--config",
        "-c",
        type=str,
        help="Path to config.toml file",
    )
    parser_modeling.add_argument(
        "--phase",
        type=str,
        help="Phase name from config file",
    )
    parser_modeling.add_argument(
        "--seeds",
        "-s",
        type=int,
        nargs="+",
        help="Shuffle seeds (e.g., 1 2 3)",
    )
    parser_modeling.add_argument(
        "--batch-size",
        "-b",
        type=int,
        default=1,
        help="Batch size for modeling (default: 1)",
    )
    parser_modeling.add_argument(
        "--enable-sqs",
        "-q",
        action="store_true",
        help="Enable SQS generation (default: False)",
    )
    parser_modeling.add_argument(
        "--iterations",
        "-i",
        type=int,
        default=1e7,
        help="Number of iterations for sqsgenerator (default: 1e7)",
    )
    parser_modeling.set_defaults(func=cmd_modeling)

    # Supercell command
    parser_supercell = subparsers.add_parser("supercell", help="Generate supercell")
    parser_supercell.add_argument(
        "--poscar",
        "-p",
        type=str,
        required=True,
        help="Path to POSCAR file",
    )
    parser_supercell.add_argument(
        "--factors",
        "-f",
        type=int,
        nargs=3,
        default=(3, 3, 3),
        metavar=("X", "Y", "Z"),
        help="Supercell factors (e.g., 3 3 3)",
    )
    parser_supercell.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser_supercell.add_argument(
        "--by-ase",
        "-a",
        action="store_true",
        help="Use ASE to make supercell (default: False)",
    )
    parser_supercell.set_defaults(func=cmd_supercell)

    # Compare command
    parser_compare = subparsers.add_parser("compare", help="Compare two POSCAR files")
    parser_compare.add_argument(
        "--poscar1",
        "-p1",
        type=str,
        required=True,
        help="Path to first POSCAR file",
    )
    parser_compare.add_argument(
        "--poscar2",
        "-p2",
        type=str,
        required=True,
        help="Path to second POSCAR file",
    )
    parser_compare.set_defaults(func=cmd_compare)

    # Merge command
    parser_merge = subparsers.add_parser("merge", help="Merge two POSCAR files")
    parser_merge.add_argument(
        "--poscar1",
        "-p1",
        type=str,
        required=True,
        help="Path to first POSCAR file",
    )
    parser_merge.add_argument(
        "--poscar2",
        "-p2",
        type=str,
        required=True,
        help="Path to second POSCAR file",
    )
    parser_merge.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser_merge.set_defaults(func=cmd_merge)

    # Separate command
    parser_separate = subparsers.add_parser("separate", help="Separate a POSCAR file by groups")
    parser_separate.add_argument(
        "--poscar",
        "-p",
        type=str,
        required=True,
        help="Path to POSCAR file",
    )
    parser_separate.add_argument(
        "--key", "-k", type=str, default="note", help="Key to group atoms (default: note)"
    )
    parser_separate.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser_separate.set_defaults(func=cmd_separate)

    # CountCN command
    parser_countcn = subparsers.add_parser("countcn", help="Calculate coordination numbers")
    parser_countcn.add_argument(
        "--name",
        "-n",
        type=str,
        default="countcn",
        help="Name for the coordination number analysis (default: countcn)",
    )
    parser_countcn.add_argument(
        "--poscar",
        "-p",
        type=str,
        required=True,
        help="Path to POSCAR file",
    )
    parser_countcn.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser_countcn.add_argument(
        "--cutoff-mult",
        "-m",
        type=float,
        default=1.1,
        help="Multiplier for cutoff radius (default: 1.1)",
    )
    parser_countcn.add_argument(
        "--parallel",
        "-j",
        type=int,
        default=2,
        help="Number of parallel processes (default: 2)",
    )
    parser_countcn.add_argument(
        "--by-ase",
        "-a",
        action="store_true",
        help="Whether to use ASE for CN calculation (default: False)",
    )
    parser_countcn.set_defaults(func=cmd_countcn)

    # Slice command
    parser_slice = subparsers.add_parser("slice", help="Slice a structure by Miller index")
    parser_slice.add_argument(
        "--name",
        "-n",
        type=str,
        default="slice",
        help="Name for the slicing analysis (default: slice)",
    )
    parser_slice.add_argument(
        "--poscar",
        "-p",
        type=str,
        required=True,
        help="Path to POSCAR file",
    )
    parser_slice.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser_slice.add_argument(
        "--miller-index",
        "-i",
        type=int,
        nargs=3,
        required=True,
        metavar=("H", "K", "L"),
        help="Miller index of the slice (e.g., 1 1 1)",
    )
    parser_slice.set_defaults(func=cmd_slice)

    # Slice to CountCN command
    parser_slice_to_countcn = subparsers.add_parser(
        "slice-to-countcn", help="Slice a structure and count CNs in each layer"
    )
    parser_slice_to_countcn.add_argument(
        "--name",
        "-n",
        type=str,
        default="slice-to-countcn",
        help="Name for the slice to CN count analysis (default: slice-to-countcn)",
    )
    parser_slice_to_countcn.add_argument(
        "--poscar",
        "-p",
        type=str,
        required=True,
        help="Path to POSCAR file",
    )
    parser_slice_to_countcn.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser_slice_to_countcn.add_argument(
        "--miller-index",
        "-i",
        type=int,
        nargs=3,
        required=True,
        metavar=("H", "K", "L"),
        help="Miller index of the slice (e.g., 1 1 1)",
    )
    parser_slice_to_countcn.set_defaults(func=cmd_slice_to_countcn)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 0

    try:
        exit_code = args.func(args)
        return exit_code
    except KeyboardInterrupt:
        logging.info("Interrupted by user")
        return 130
    except Exception as e:
        logging.error(f"Error: {e}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
