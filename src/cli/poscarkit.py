# src/cli/poscarkit.py

import sys
import argparse
import logging
import tomllib
import traceback
from pathlib import Path
from typing import Any

from src.modeling.base import SimplePoscar
from src.modeling.countcn import CNCounter
from src.modeling.slice import Slicer
from src.workflow.modeling import run_modeling
from src.modeling.supercell import unitcell2file, supercell2file
from src.workflow.slice_to_countcn import slice2files_with_countcn
from src.workflow.thermo import run_thermo
from src.config import VERSION, CONTACT, DEVELOPER, normalize_config_keys


WORKDIR = Path.cwd()


def cmd_help(args: argparse.Namespace) -> int:
    print(
        f"""
================================== POSCARKIT ==================================
{DEVELOPER}
{VERSION}
{CONTACT}

A tool for modeling structure POSCAR files, based on Sublattice Occupying 
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
  import-to-model  Import SOFs from CSV/XLSX and run modeling
  surface       Generate surface slabs from bulk POSCAR
  thermo        Calculate Sconf and DeltaG from SOF data + TDB

EXAMPLES:
  poscarkit help
  poscarkit modeling --poscar POSCAR.vasp --factors 3 3 3 --name my_model
  poscarkit countcn --poscar POSCAR.vasp --name my_cn --outdir output/
  poscarkit slice --poscar POSCAR.vasp --miller-index 1 1 1 --outdir output/
  poscarkit slice-to-countcn --poscar POSCAR.vasp --miller-index 1 1 1 --outdir output/
  poscarkit supercell --poscar POSCAR.vasp --factors 2 2 2 --outdir output/
  poscarkit compare --poscar1 POSCAR1.vasp --poscar2 POSCAR2.vasp
  poscarkit merge --poscars POSCAR1.vasp POSCAR2.vasp POSCAR3.vasp --outdir output/
  poscarkit separate --poscar POSCAR.vasp --key note --outdir output/
  poscarkit import-to-model --csv tc_exps.csv --phase fcc -t 473 873 1273
  poscarkit import-to-model --csv pandat-data.csv --phase fcc -t 873 -o print
  poscarkit thermo --data sof_data.xlsx --tdb database.TDB --outdir output/
"""
    )
    return 0


def cmd_modeling(args: argparse.Namespace) -> int:
    name = args.name or "modeling"
    poscar = Path(args.poscar) if args.poscar else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    factors = tuple(args.factors)
    phase = args.phase.upper() if args.phase else args.phase
    config = args.config
    seeds = args.seeds if args.seeds else [None]
    batch_size = args.batch_size
    enable_sqs = args.enable_sqs
    iterations = args.iterations

    structure_info: dict[str, Any] = {}

    if config:
        with open(config, "rb") as f:
            cfg = normalize_config_keys(tomllib.load(f))
    else:
        cfg = {}
    structure_info = cfg.get(phase.upper() if phase else "", {})

    outdir.mkdir(parents=True, exist_ok=True)

    if poscar and not poscar.is_file():
        logging.error(f"POSCAR file not found: {poscar}")
        return 1

    if not poscar and not structure_info:
        logging.error("Either --poscar or (--config and --phase) must be provided")
        return 1

    results = run_modeling(
        name=name,
        poscar=Path(poscar or ""),
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


def cmd_countcn(args: argparse.Namespace) -> int:
    name = args.name or "countcn"
    poscar = Path(args.poscar) if args.poscar else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    cutoff_mult = args.cutoff_mult
    parallel = args.parallel
    by_ase = args.by_ase
    pbc = getattr(args, "pbc", False)

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
        pbc=pbc,
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
    pbc = getattr(args, "pbc", False)
    by_ase = getattr(args, "by_ase", False)

    if not poscar or not poscar.is_file():
        logging.error(f"POSCAR file not found: {poscar}")
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    results = slice2files_with_countcn(
        name=name, poscar=poscar, outdir=outdir, miller_index=miller_index,
        pbc=pbc, by_ase=by_ase,
    )

    logging.info(f"Slice to CountCN completed. Results saved to:")
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
    poscars = [Path(p) for p in args.poscars] if args.poscars else []
    outdir = Path(args.outdir) if args.outdir else Path("output")

    if len(poscars) < 2:
        logging.error("At least 2 POSCAR files are required for merging")
        return 1

    for poscar in poscars:
        if not poscar.is_file():
            logging.error(f"POSCAR file not found: {poscar}")
            return 1

    outdir.mkdir(parents=True, exist_ok=True)

    merged_file = SimplePoscar.merge_poscar(poscars, outdir)
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


def cmd_thermo(args: argparse.Namespace) -> int:
    data_path = Path(args.data) if args.data else None
    tdb_path = Path(args.tdb) if args.tdb else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    name = args.name or "thermo"

    if not data_path or not data_path.is_file():
        logging.error(f"Data file not found: {data_path}")
        return 1
    if not tdb_path or not tdb_path.is_file():
        logging.error(f"TDB file not found: {tdb_path}")
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    csv_path, png_path = run_thermo(
        data_path=str(data_path),
        tdb_path=str(tdb_path),
        name=name,
        outdir=str(outdir),
    )

    logging.info(f"Results saved: {csv_path}, {png_path}")
    return 0


def cmd_import_to_model(args: argparse.Namespace) -> int:
    csv_path = Path(args.csv)
    if not csv_path.is_file():
        logging.error(f"CSV/XLSX file not found: {csv_path}")
        return 1

    config_path = args.config or str(WORKDIR / "config.toml")
    if not Path(config_path).is_file():
        logging.error(f"Config file not found: {config_path}")
        return 1

    temperatures = args.temperatures or []
    if not temperatures:
        logging.error("No temperatures specified. Use --temperatures / -t")
        return 1

    from workflow.import_to_model import run_import_to_model

    results = run_import_to_model(
        csv_path=str(csv_path),
        config_path=config_path,
        phase=args.phase,
        temperatures=list(temperatures),
        sublattice_map=args.sublattice_map,
        outdir=args.outdir or str(WORKDIR / "output"),
        name=args.name or "modeling",
        factors=tuple(args.factors) if args.factors else (3, 3, 3),
        seeds=args.seeds,
        batch_size=args.batch_size or 1,
        enable_sqs=args.enable_sqs or False,
        iterations=args.iterations or int(1e7),
        output_mode=args.output or "run",
    )
    logging.info(f"Import-to-model completed. Files: {len(results)}")
    return 0


def cmd_surface(args: argparse.Namespace) -> int:
    poscar = Path(args.poscar) if args.poscar else None
    outdir = Path(args.outdir) if args.outdir else Path("output")
    miller = tuple(args.miller)
    layers = args.layers
    vacuum = args.vacuum
    fix_layers = args.fix_layers
    fix_z_only = args.fix_z_only

    if not poscar or not poscar.is_file():
        logging.error(f"POSCAR file not found: {poscar}")
        return 1

    outdir.mkdir(parents=True, exist_ok=True)

    from src.modeling.surface import SurfaceBuilder

    builder = SurfaceBuilder(
        poscar=poscar,
        miller=miller,
        layers=layers,
        vacuum=vacuum,
        fix_layers=fix_layers,
        fix_z_only=fix_z_only,
        outdir=outdir,
        precision=getattr(args, "precision", 2),
    )

    try:
        builder._transform_cell()
        builder._identify_layers()
        builder._find_gaps()
        results = builder.build_all(outdir)

        logging.info(f"Surface slab generation completed. Generated {len(results)} slabs:")
        for r in results:
            logging.info(f"  - {r}")
        return 0
    except ValueError as e:
        logging.error(str(e))
        return 1


def main() -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s[%(levelname)s]%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

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
    parser_countcn.add_argument(
        "--pbc",
        action="store_true",
        help="Enable periodic boundary conditions for CN counting",
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
        help="Name for the slice to CountCN analysis (default: slice-to-countcn)",
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
    parser_slice_to_countcn.add_argument(
        "--pbc",
        action="store_true",
        help="Enable periodic boundary conditions for CN counting in layers",
    )
    parser_slice_to_countcn.add_argument(
        "--by-ase",
        "-a",
        action="store_true",
        help="Use ASE backend for CN counting (default: False)",
    )
    parser_slice_to_countcn.set_defaults(func=cmd_slice_to_countcn)

    # Surface command
    parser_surface = subparsers.add_parser(
        "surface",
        help="Generate surface slabs from bulk POSCAR",
    )
    parser_surface.add_argument(
        "poscar",
        type=str,
        help="Path to bulk POSCAR file",
    )
    parser_surface.add_argument(
        "--miller",
        "-m",
        type=int,
        nargs=3,
        default=(0, 0, 1),
        metavar=("H", "K", "L"),
        help="Miller indices (default: 0 0 1)",
    )
    parser_surface.add_argument(
        "--layers",
        "-l",
        type=int,
        default=3,
        help="Number of slab layers, outputs N and N+1 (default: 3)",
    )
    parser_surface.add_argument(
        "--vacuum",
        "-v",
        type=float,
        default=15.0,
        help="Total vacuum thickness in Angstrom (default: 15.0)",
    )
    parser_surface.add_argument(
        "--fix-layers",
        type=int,
        default=None,
        help="Manual override for fixed bottom layers (default: auto)",
    )
    parser_surface.add_argument(
        "--fix-z-only",
        action="store_true",
        help="Fix only z-direction instead of all (default: FFF)",
    )
    parser_surface.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser_surface.add_argument(
        "--precision",
        type=int,
        default=2,
        help="Decimal precision for z-coordinate grouping (default: 2)",
    )
    parser_surface.set_defaults(func=cmd_surface)

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
    parser_merge = subparsers.add_parser("merge", help="Merge multiple POSCAR files")
    parser_merge.add_argument(
        "--poscars",
        "-p",
        type=str,
        nargs="+",
        required=True,
        help="Paths to POSCAR files to merge (at least 2)",
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

    # Import SOFs command
    parser_import = subparsers.add_parser(
        "import-to-model",
        help="Import SOFs from ThermoCalc/Pandat CSV/XLSX files",
    )
    parser_import.add_argument(
        "--csv", "-c", type=str, required=True, help="Path to CSV or XLSX file"
    )
    parser_import.add_argument(
        "--config",
        type=str,
        help="Path to config.toml (default: ./config.toml)",
    )
    parser_import.add_argument(
        "--phase",
        "-p",
        type=str,
        required=True,
        help="Phase name (fcc, bcc, hcp)",
    )
    parser_import.add_argument(
        "--temperatures",
        "-t",
        type=float,
        nargs="+",
        help="Temperatures to import (space-separated, e.g. 473 873 1273)",
    )
    parser_import.add_argument(
        "--sublattice-map",
        type=str,
        help="Override CALPHAD to Wyckoff mapping (e.g. '1:1a,2:3c')",
    )
    parser_import.add_argument(
        "--output",
        "-o",
        type=str,
        choices=["run", "save-config", "print"],
        default="run",
        help="Output mode: run modeling, save config, or print to stdout (default: run)",
    )
    parser_import.add_argument(
        "--name",
        "-n",
        type=str,
        default="modeling",
        help="Work name prefix (default: modeling)",
    )
    parser_import.add_argument(
        "--factors",
        type=int,
        nargs=3,
        default=(3, 3, 3),
        metavar=("X", "Y", "Z"),
        help="Supercell factors (default: 3 3 3)",
    )
    parser_import.add_argument(
        "--outdir", type=str, default="output", help="Output directory"
    )
    parser_import.add_argument(
        "--seeds", type=int, nargs="+", help="Shuffle seeds"
    )
    parser_import.add_argument(
        "--batch-size",
        "-b",
        type=int,
        default=1,
        help="Batch size for modeling (default: 1)",
    )
    parser_import.add_argument(
        "--enable-sqs",
        "-q",
        action="store_true",
        help="Enable SQS generation",
    )
    parser_import.add_argument(
        "--iterations",
        "-i",
        type=int,
        default=int(1e7),
        help="SQS iterations (default: 1e7)",
    )
    parser_import.set_defaults(func=cmd_import_to_model)

    # Thermo command
    parser_thermo = subparsers.add_parser(
        "thermo", help="Calculate Sconf and DeltaG from SOF data + TDB"
    )
    parser_thermo.add_argument(
        "--data", "-d", type=str, required=True, help="Path to data file (xlsx/csv)"
    )
    parser_thermo.add_argument(
        "--tdb", "-t", type=str, required=True, help="Path to TDB thermodynamic database"
    )
    parser_thermo.add_argument(
        "--name",
        "-n",
        type=str,
        default="thermo",
        help="Work name for output subdirectory (default: thermo)",
    )
    parser_thermo.add_argument(
        "--outdir",
        "-o",
        type=str,
        default="output",
        help="Output directory (default: output)",
    )
    parser_thermo.set_defaults(func=cmd_thermo)

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
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
