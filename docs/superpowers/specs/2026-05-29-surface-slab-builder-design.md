# Surface Slab Builder Design

**Date:** 2026-05-29
**Status:** Approved

## Overview

Add a surface slab model builder to poscarkit that generates 3-4 layer asymmetric
surface slabs from bulk POSCAR files. Targets high-entropy alloy (HEA) research
on BCC, FCC, HCP, and custom phases. Replaces the manual Materials Studio
workflow of importing bulk structure → cutting surfaces → selecting layers →
adding vacuum → exporting POSCAR.

## Scope

- New `poscarkit surface` CLI subcommand
- New `SurfaceBuilder` class in `src/modeling/surface.py`
- Fix `struct2atoms`/`atoms2struct` bridge to preserve `note`, `constr`, `meta`
- Fix periodic boundary bug in `Slicer.group_by_normal`
- Unit tests for all new and modified code

## Architecture

### Approach: Hybrid (ASE geometry + custom business logic)

ASE's `cut()` handles the hard Miller-index cell transformation (proven to
preserve custom per-atom arrays). All HEA-specific logic — layer identification,
gap detection, slab assembly, vacuum, selective dynamics, metrics — runs on
native `Struct` objects with full `note`/`constr`/`meta` preservation.

### Data Flow

```
POSCAR → Struct
  → struct2atoms (preserves note/constr/meta in atoms.arrays)
  → ASE cut (transforms cell, z-axis = surface normal)
  → atoms2struct (restores note/constr/meta from atoms.arrays)
  → Layer identification (fractional z + circular gap detection)
  → For each gap:
      → Select N / N+1 layers → sub-Struct
      → Unwrap periodic z-coordinates → continuous Cartesian positions
      → Build new cell (a,b in-plane; c = slab_thickness + vacuum)
      → Add vacuum (bottom 2Å + top remaining)
      → Apply Selective Dynamics (fix bottom layers)
      → Write POSCAR
  → Write summary CSV
```

## Component Design

### 1. Bridge Fix: `struct2atoms` / `atoms2struct` (base.py)

**struct2atoms:** Store `note`, `meta`, `constr` in `atoms.new_array()`.

**atoms2struct:** Restore from `atoms.arrays` with fallback to defaults for
backward compatibility with bare ASE Atoms objects.

### 2. Layer Identification Algorithm

1. Get fractional z-coordinates of all atoms after cell transformation
2. Sort by z: `z_0, z_1, ..., z_{n-1}`
3. Compute gaps: `gaps[i] = z_{i+1} - z_i`, plus wrap gap `1.0 - z_{n-1} + z_0`
4. Threshold = `median(gaps) * 3`
5. Large gaps → layer boundaries; small gaps → intra-layer variation
6. If wrap gap is small → merge first and last groups (same layer crossing boundary)
7. Layers ordered by fractional z ascending

### 3. Slab Assembly (per gap)

For gap_i between layer_i and layer_{i+1}:

- Select N layers from layer_i downward (cyclic indexing): `layer_i, layer_{i-1}, ..., layer_{i-N+1}`
- Unwrap z: adjust fractional z by ±1 as needed for continuity from gap_i downward
- `slab_thickness = gap_i_z_cartesian - min_bottom_z_cartesian`
- New cell: `a, b` unchanged; `c = (0, 0, slab_thickness + vacuum_bottom + vacuum_top)`
- Translate atoms: `z += vacuum_bottom - z_min`
- Apply constraints to bottom `fix_layers` layers

### 4. Vacuum Distribution

- `vacuum_bottom = 2.0 Å` (thin gap behind fixed layers)
- `vacuum_top = user_vacuum - 2.0 Å`
- User parameter `--vacuum` defaults to 15.0 Å

### 5. Selective Dynamics

| `--fix-layers` | `--fix-z-only` | Constraint |
|---|---|---|
| auto (default) | off | Bottom layers: F F F, rest: T T T |
| auto (default) | on | Bottom layers: T T F, rest: T T T |
| explicit N | off | Bottom N layers: F F F |
| explicit N | on | Bottom N layers: T T F |

Auto-fix rule: floor((layers+1)/2) layers fixed (3→2, 4→2, 5→3, ...).

### 6. Output

**Directory:** `<outdir>/<name>-slab-<miller>/`

**Files (flat in directory):**
- `<name>-slab-<miller>-term<XX>-layers<N>.vasp` — slab POSCAR
- `<name>-slab-<miller>-term<XX>-layers<N+1>.vasp` — slab POSCAR (N+1 variant)
- `summary.csv` — all slabs with metrics

**termination_id:** term01, term02, ... sorted by gap fractional z ascending.

### 7. Summary CSV Columns

| Column | Description |
|---|---|
| `filename` | Output filename |
| `termination_id` | term01, term02, ... |
| `gap_position_z` | Cutting position in fractional z |
| `n_layers` | Number of layers in slab |
| `n_atoms` | Total atom count |
| `composition` | Overall element composition |
| `layer_compositions` | Per-layer composition (pipe-separated) |
| `composition_deviation` | sum(abs(surface%-bulk%)) / n_elements |
| `surface_energy_est` | Bond-energy model estimate |
| `dipole_moment` | Macroscopic dipole moment (D) |
| `top_layer_z_std` | Surface roughness: std of top-layer z |
| `fix_layers` | Number of fixed bottom layers |
| `fix_mode` | FFF or TTF |
| `vacuum_top` | Top vacuum thickness (Å) |
| `vacuum_bottom` | Bottom vacuum thickness (Å) |

### 8. SurfaceBuilder Class

File: `src/modeling/surface.py`

```
SurfaceBuilder
├── __init__(poscar, miller, layers, vacuum, fix_layers, fix_z_only, orthogonal, precision)
├── _transform_cell()       → Struct
├── _identify_layers()      → list[Layer]
├── _find_gaps()            → list[float]
├── _unwrap_zs(atom_list, gap_z) → np.ndarray
├── _build_slab(gap_z, n_layers) → Struct
├── _add_vacuum(slab)       → Struct
├── _add_constraints(slab)  → Struct
├── _compute_metrics(slab, bulk) → dict
├── build_all(outdir)       → list[Path]
└── _write_summary(slabs, metrics) → Path
```

`Layer` dataclass:
```python
@dataclass
class Layer:
    index: int
    z_centroid: float
    atoms: list[Atom]
    composition: str
```

### 9. CLI Interface

```bash
poscarkit surface <poscar> [options]

Options:
  --miller H K L      Miller indices (default: 0 0 1)
  --layers N          Number of layers, outputs N and N+1 (default: 3)
  --vacuum V          Total vacuum in Å (default: 15.0)
  --fix-layers N      Manual override for fixed bottom layers (default: auto)
  --fix-z-only        Fix only z-direction instead of all (default: FFF)
  --orthogonal        Force orthogonal cell for FCC(110)/BCC(110)
  --outdir DIR        Output directory (default: ./output/)
  --precision P       Decimal precision for z-coordinate grouping (default: 2)
```

### 10. Fix: Slicer.group_by_normal Periodic Bug

Current code uses Cartesian coordinates for layer grouping, which misidentifies
atoms near cell boundaries as separate layers. Fix: use fractional coordinates
along the surface normal with circular gap detection, same algorithm as the new
SurfaceBuilder layer identification.

## Files Changed

| File | Action | Description |
|---|---|---|
| `src/modeling/base.py` | Modify | Fix struct2atoms/atoms2struct bridge |
| `src/modeling/slice.py` | Modify | Fix group_by_normal periodic bug |
| `src/modeling/surface.py` | New | SurfaceBuilder class |
| `src/cli/poscarkit.py` | Modify | Register `surface` subcommand |

## Testing

### New: `tests/modeling/test_surface.py`

| Test | What It Verifies |
|---|---|
| `test_transform_cell_fcc111` | FCC → (111) correctly orients z-axis |
| `test_transform_cell_bcc100` | BCC → (100) transformation |
| `test_transform_cell_custom` | Arbitrary Miller index |
| `test_identify_layers` | Correct layer count and z-centroids |
| `test_find_gaps` | Gap count = layer count, positions between layers |
| `test_build_slab_3layers` | Correct atom count and symbols |
| `test_build_slab_4layers` | N+1 auto-output |
| `test_vacuum_distribution` | bottom=2Å, top=total-2Å |
| `test_constraints_fff` | Default: bottom layers FFF |
| `test_constraints_ttf` | --fix-z-only: bottom layers TTF |
| `test_note_preserved` | note/field survives round-trip |
| `test_layers_periodic_crossing` | Layer crossing cell boundary handled correctly |
| `test_summary_csv` | All columns present, values in valid ranges |

### Modified: `tests/modeling/test_base.py`

- `test_struct_atoms_roundtrip`: Struct → Atoms → Struct preserves note, constr, meta

### Modified: `tests/modeling/test_slice.py`

- `test_group_by_normal_periodic`: Same-layer atoms across cell boundary grouped together

## Non-Goals

- Symmetric slab construction (HEA surfaces use asymmetric slabs)
- Adsorbate placement (separate future feature)
- DFT input file generation (INCAR, KPOINTS, etc.)
- Integration with `run_modeling` workflow (deferred to future PR)
- GUI or interactive mode for surface selection
