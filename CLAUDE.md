# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & test

```bash
pip install -e .                        # editable install (needed for import paths)
python -m unittest discover tests       # run all tests
python -m unittest tests.modeling.test_simple_poscar   # single test module
```

Tests use `unittest`, not pytest (despite the pytest dependency). The `[tool.pytest.ini_options]` in pyproject.toml only sets `pythonpath`.

## Entry points

- `main.py` — dispatches to CLI mode (when `sys.argv` has args) or interactive terminal menu (no args)
- CLI: `poscarkit <subcommand>` (8 subcommands: help, modeling, countcn, slice, slice-to-countcn, supercell, compare, merge, separate)
- Interactive: `poscarkit_interact` or `python main.py` with no args, presents numbered menu

## Architecture

### Core data layer (`src/modeling/base.py`)

`Atom` is a slots-based dataclass with `index`, `symbol`, `coord` (np.ndarray), `constr` (selective dynamics flags), `note` (sublattice identifier like `"1a-Au"`), and `meta` (batch index).

`Struct` wraps a list of `Atom` + cell matrix. It provides grouping (`group_structs`, `classify`), sorting, coordinate switching (direct/cartesian), duplicate removal, and constraint management. The `compare` method compares two structs cell, symbols, and atom positions.

`SimplePoscar` is a static-method-only namespace class for POSCAR file I/O (`read_poscar`, `write_poscar`) and `Struct <-> Atoms` conversion (`struct2atoms`, `atoms2struct`). POSCAR operations (compare, merge, separate) are also static methods here.

### Modeling engines

Supercell generation (`src/modeling/supercell.py`):
- `make_supercell()` — manual supercell via broadcasting direct coordinates across grid indices
- `unitcell2file()` — build a unitcell POSCAR from `structure_info` dict (from config.toml)
- `supercell2file()` — CLI-facing wrapper with optional ASE backend (`by_ase=True`)

Atomic allocation (`src/modeling/model.py`):
- `ModelStruct` orchestrates allocation. `_gen_site_integers()` converts fractional SOFs from config into integer atom counts per sublattice.
- `model_by_shuffle()` — randomly shuffles atoms within each sublattice, assigns new symbols according to integer fractions. Uses `multiprocessing.Pool` for batch parallelism.
- `model_by_sqsgen()` — delegates to `sqsgenerator` package for Special Quasirandom Structures optimization. Each sublattice is processed independently then merged.

### Analysis modules

`CNCounter` (`src/modeling/countcn.py`):
- Cutoff detection via pairwise distance gradient analysis (`detect_cutoff`)
- Two CN backends: custom KDTree-based (`calculate_cn`) or ASE `NeighborList` (`calculate_cn_ase`)
- Outputs: per-pair POSCAR files, CSV tables, faceted/stacked histograms, heatmap (uses matplotlib Agg backend)

`Slicer` (`src/modeling/slice.py`):
- Transforms structure to new basis defined by Miller index (`cut` from `ase.build.tools`)
- Groups atoms by projection onto surface normal, exports layers as POSCAR + PNG scatter plot + Excel

### Workflow layer

`src/workflow/modeling.py::run_modeling()` — chains supercell generation + atomic allocation, saves each batch result as POSCAR.

`src/workflow/slice_to_countcn.py::slice2files_with_countcn()` — slices then runs CN counting on each layer.

### CLI layer

`src/cli/poscarkit.py` — argparse-based with subcommands; each `cmd_*` function handles one operation. `main()` dispatches via `args.func(args)`.

`src/cli/poscarkit_interact.py` — `PoscarkitInteract` class runs an interactive REPL menu. Each option is a method (`run_modeling`, `run_countcn`, etc.) that prompts for missing inputs.

### Configuration

`config.toml` defines crystal prototypes (bcc, fcc, hcp) with cell vectors and sublattice atom positions, plus SOFs for each sublattice site. Custom phases can be added. The `src/config/__init__.py` module holds version/contact metadata and the `DEFAULT_CONFIG` template string.

## Key conventions

- Coordinates are always stored as **direct** (fractional) in `Struct` after reading. `get_coords()` switches between direct and cartesian.
- The `note` field on `Atom` encodes sublattice membership as `"<site>-<element>"` (e.g. `"1a-Au"`, `"3c-Cu"`). This is parsed by regex in allocation and grouping operations.
- Matplotlib uses the `Agg` non-interactive backend — don't import `plt` directly without setting the backend first.
- Imports use the `src.` prefix (e.g. `from src.modeling.base import Struct`) because `pip install -e .` makes `src` a package root.
