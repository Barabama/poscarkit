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

- `main.py` — dispatches to **CLI mode** (when `sys.argv` has args) or **GUI mode** (no args, launches `PoscaKitGUI`)
- CLI: `poscarkit <subcommand>` (12 subcommands: help, modeling, countcn, slice, slice-to-countcn, supercell, compare, merge, separate, import-to-model, thermo, surface)
- Interactive terminal: `poscarkit_interact` runs the numbered-menu REPL (`PoscarkitInteract`), now with 13 options including Import to Model, Thermo, and Surface

## Architecture

### Core data layer (`src/modeling/base.py`)

`Atom` is a slots-based dataclass with `index`, `symbol`, `coord` (np.ndarray), `constr` (selective dynamics flags as `["T","T","T"]`), `note` (sublattice identifier like `"1a-Au"`), and `meta` (batch index).

`Struct` wraps a list of `Atom` + cell matrix. It provides grouping (`group_structs`, `classify`), sorting (by symbol, coord, x/y/z, note, meta), coordinate switching (direct/cartesian), duplicate removal, and constraint management. The `compare` method compares two structs' cell, symbols, and atom positions. `constr`, `note`, and `meta` are preserved through `struct2atoms`/`atoms2struct` round-trips.

`SimplePoscar` is a static-method-only namespace class for POSCAR file I/O (`read_poscar`, `write_poscar`) and `Struct <-> Atoms` conversion (`struct2atoms`, `atoms2struct`). POSCAR operations (compare, merge, separate) are also static methods here. Sorting now prioritizes atom notes.

### Modeling engines

Supercell generation (`src/modeling/supercell.py`):
- `make_supercell()` — manual supercell via broadcasting direct coordinates across grid indices
- `unitcell2file()` — build a unitcell POSCAR from `structure_info` dict (from config.toml)
- `supercell2file()` — CLI-facing wrapper with optional ASE backend (`by_ase=True`)

Atomic allocation (`src/modeling/model.py`):
- `ModelStruct` orchestrates allocation. `_gen_site_integers()` converts fractional SOFs from config into integer atom counts per sublattice.
- `model_by_shuffle()` — randomly shuffles atoms within each sublattice, assigns new symbols according to integer fractions. Uses `multiprocessing.Pool` for batch parallelism.
- `model_by_sqsgen()` — delegates to `sqsgenerator` package for Special Quasirandom Structures optimization. Each sublattice is processed independently then merged.
- SOF→integer conversion log is saved as CSV (`<name>-sof-integers.csv`) via `run_modeling()`.

Surface slab builder (`src/modeling/surface.py`):
- `SurfaceBuilder` — generates asymmetric surface slabs from bulk POSCAR. Uses ASE `cut()` for cell transformation, rounded-projection z-grouping for layer identification, gap-based cutting plane detection, coordinate unwrapping for slab assembly, vacuum distribution (default 15 Å), and selective dynamics (FFF or TTF). Outputs N and N+1 layer slabs with metrics summary CSV (composition deviation, dipole moment, surface energy estimate, layer compositions, etc.).
- `Layer` dataclass holds per-layer atoms, z-centroid, composition string.

### Analysis modules

`CNCounter` (`src/modeling/countcn.py`):
- Cutoff detection via pairwise distance gradient analysis (`detect_cutoff`), with optional PBC (minimum-image-convention distances via `_pdist_mic`) and `seed` parameter for reproducible sampling
- Two CN backends: custom KDTree-based (`calculate_cn`) or ASE `NeighborList` (`calculate_cn_ase`)
- Cumulative CN distribution (`CN >= n`) via `generate_cn_cumulative_structs()`
- Outputs: per-pair POSCAR files, CSV tables, faceted/stacked histograms, heatmap (uses matplotlib Agg backend)

`Slicer` (`src/modeling/slice.py`):
- Transforms structure to new basis defined by Miller index (`cut` from `ase.build.tools`)
- Groups atoms by projection onto surface normal using rounded fractional z-coordinates
- Exports layers as POSCAR + PNG scatter plot (with 10 Å scale bar and semi-transparent backdrop) + Excel

### IO layer (`src/io/`)

Format auto-detection and reader dispatch (`src/io/readers/__init__.py`):
- `detect_format(path)` — detects `tc_exps`, `pandat`, or `teacher` from column name patterns
- `read(path, phase_hint)` — dispatches to the appropriate reader, returns `SOFData`

Intermediate representation (`src/io/ir.py`):
- `SOFData` — unified dataclass: `phase`, `site_ratios`, `T`, `Y_subl` (per-sublattice site fractions), `composition`, `elements`, optional `G_real`
- `get_sofs_at(ir, T, sublattice_map)` — extracts SOFs at a temperature, mapping CALPHAD sublattice indices to Wyckoff site names
- `build_structure_info(config_cfg, sofs_by_site)` — merges config geometry with imported SOF data for `run_modeling()`
- `parse_sublattice_map("1:1a,2:3c")` — parses CALPHAD→Wyckoff mapping strings
- `PHASE_SITE_RATIOS` and `PHASE_SUBLATTICE_MAP` — fallback maps for FCC/BCC/HCP

Readers (`src/io/readers/`):
- `tc_exps.py` — ThermoCalc experimental format
- `pandat.py` — Pandat format (multi-sheet xlsx with T-column merge)
- `teacher.py` — Educational/simple format

### Thermodynamics (`src/thermo/`)

- `tdb.py` — TDB file parser: extracts FUNCTION, PARAMETER (G type, order=0), PHASE records. `resolve_expression()` replaces SERXX#/GHSERXX# references with resolved function bodies.
- `sconf.py` — `calc_sconf()` for real SOFs (N-sublattice formula), `calc_sconf_random()` for random mixing
- `deltaG.py` — `calc_deltaG_random()` evaluates weighted end-member expressions + `R*T*Sconf` term. Uses restricted numpy `eval` (no sympy dependency).
- `plot.py` — Dual-panel Sconf-T + DeltaG-T chart (matplotlib Agg backend)

### Workflow layer

`src/workflow/modeling.py::run_modeling()` — chains supercell generation + atomic allocation, saves each batch result as POSCAR, writes SOF→integer CSV log.

`src/workflow/slice_to_countcn.py::slice2files_with_countcn()` — slices then runs CN counting on each layer. Supports `pbc` and `by_ase` parameters.

`src/workflow/import_to_model.py::run_import_to_model()` — reads CSV/XLSX SOF data, maps to Wyckoff sites, optionally patches config, runs modeling for each temperature. Three output modes: `run` (generate POSCARs), `save-config` (write config snippet), `print` (stdout JSON).

`src/workflow/thermo.py::run_thermo()` — reads SOF data + TDB, computes Sconf and DeltaG, saves CSV + PNG plot.

### GUI layer (`src/gui/`)

`app.py::PoscaKitGUI` — Tkinter main window:
- Sidebar with 13 color-coded function buttons (Config, Modeling, Import to Model, Count CN, Slice, Slice to CountCN, Surface, Thermo, Supercell, Compare, Merge, Separate, About)
- Dynamic form area with scrollable canvas; forms built by `forms.py` builders, each returning a `get_args()` callable
- Thread-safe log panel (`_LogHandler` redirects Python logging to tk.Text, max 5000 lines)
- Tasks run in background daemon threads, dispatching to `src.cli.poscarkit.cmd_*` handlers
- Config I/O: load/save TOML key-value pairs; About panel reads metadata via `importlib.metadata`

`forms.py` — Form builders with shared widget helpers (`_entry`, `_file_row`, `_dir_row`, `_combo`, `_int_entries`, `_checkbox`). Key widgets:
- `_sof_editor` — editable SOF table with per-site sum validation (green = 1.0, red = mismatch), add/remove rows, write-back to TOML
- Phase list is dynamically discovered from config dict (standard phases FCC/BCC/HCP always listed first)
- `_link_ase_pbc` — auto-checks PBC when ASE backend is selected (ASE always uses PBC)

### CLI layer

`src/cli/poscarkit.py` — argparse-based with 12 subcommands; each `cmd_*` function handles one operation. `main()` dispatches via `args.func(args)`.

`src/cli/poscarkit_interact.py` — `PoscarkitInteract` class runs an interactive REPL menu (13 options). Each option is a method that prompts for missing inputs. Uses a `Config` TypedDict for config state.

### Configuration

`config.toml` defines crystal prototypes (BCC, FCC, HCP) with cell vectors and sublattice atom positions, plus SOFs for each sublattice site. Custom phases can be added. The `src/config/__init__.py` module holds version/contact metadata, `normalize_config_keys()` (merges case-variant phase keys to uppercase), and the `DEFAULT_CONFIG` template string. Version is read via `importlib.metadata.version("poscarkit")` at runtime.

## Key conventions

- Coordinates are always stored as **direct** (fractional) in `Struct` after reading. `get_coords()` switches between direct and cartesian.
- The `note` field on `Atom` encodes sublattice membership as `"<site>-<element>"` (e.g. `"1a-Au"`, `"3c-Cu"`). This is parsed by regex in allocation and grouping operations.
- Phase keys in config are normalized to uppercase via `normalize_config_keys()`; config lookup always uses uppercase keys.
- Matplotlib uses the `Agg` non-interactive backend — don't import `plt` directly without setting the backend first.
- Imports use the `src.` prefix (e.g. `from src.modeling.base import Struct`) because `pip install -e .` makes `src` a package root.
- Heavy libraries (ASE, matplotlib, pandas) are lazy-imported at point of use, not at module top-level.
- The `constr` field on `Atom` is a list of `"T"`/`"F"` flags for selective dynamics (`["T","T","T"]` = free, `["F","F","F"]` = fixed, `["T","T","F"]` = fix z-only).
