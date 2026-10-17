# src/modeling/layers.py

"""Layer identification for Miller-index-transformed structures.

Both :class:`src.modeling.slice.Slicer` and
:class:`src.modeling.surface.SurfaceBuilder` cut a bulk structure along a
Miller plane and then have to split the resulting cell into atomic layers.

Rounding the projection onto the surface normal to a fixed number of decimal
places is not a valid way to do that.  A thermally distorted layer is a
continuous distribution roughly 0.07 A wide, so a fixed 0.01 A bin shreds it
into ~7 fragments, while a 0.7 A bin merges neighbouring layers and still
splits any layer whose centre happens to straddle a bin boundary.  The layer
count then depends on where the crystal sits relative to an arbitrary decimal
grid, which means a rigid translation of an identical structure changes the
result.

The helpers here cluster the sorted projections instead and split at real
gaps.  That is invariant under rigid translations and is expressed in
Angstroms rather than in decimal places.
"""

import logging

import numpy as np
from ase.build.tools import cut

#: Default gap (Angstrom) used to separate atomic layers.  Interlayer
#: spacings are typically > 1 A while thermal/static broadening is < 0.1 A,
#: so any value in between works; 0.3 A sits comfortably inside that window.
DEFAULT_LAYER_GAP = 0.3

#: Tolerances tried, in order, when cutting a commensurate cell.  ASE's
#: ``cut`` masks atoms with a scaled tolerance of ``0.1 * tolerance``.  The
#: smallest value conserves every atom for most cells; a few orientations need
#: a looser value so that atoms sitting exactly on the upper boundary are
#: picked up by a periodic image instead of being dropped.
CUT_TOLERANCES = (1e-6, 1e-5, 1e-4, 1e-3, 1e-2)


def scaled_basis(miller_index: tuple[int, int, int]) -> np.ndarray:
    """Return commensurate basis vectors in scaled (fractional) coordinates.

    The third vector is the surface normal.  These vectors are passed straight
    to :func:`ase.build.tools.cut`, which expects *scaled* coordinates -- it
    computes ``newcell = np.dot([a, b, c], cell)``.  Passing normalised
    Cartesian unit vectors there yields a non-commensurate cell whose third
    edge is only one lattice constant long, which truncates the structure
    mid-plane.

    Args:
        miller_index: Miller index ``(h, k, l)``.

    Returns:
        ``(3, 3)`` array of basis vectors, one per row.
    """
    miller_index = tuple(int(v) for v in miller_index)

    basis_map = {
        (0, 0, 1): [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        (1, 1, 0): [(0, 0, 1), (0.5, -0.5, 0), (1, 1, 0)],
        (1, 1, 1): [(0.5, -0.5, 0), (0.5, 0.5, -1), (1, 1, 1)],
    }
    if miller_index in basis_map:
        return np.array(basis_map[miller_index], dtype=float)

    # General Miller index: build two vectors perpendicular to n, then n.
    n = np.array(miller_index, dtype=float)
    t0 = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < abs(n[1]) else np.array([0.0, 1.0, 0.0])
    b1 = np.cross(n, t0)
    b2 = np.cross(n, b1)
    return np.array([b1, b2, n], dtype=float)


def cut_commensurate(
    atoms,
    basis: np.ndarray,
    tolerances: tuple[float, ...] = CUT_TOLERANCES,
):
    """Cut ``atoms`` with ``basis`` while conserving the atom count.

    The expected count is ``len(atoms) * |det(basis)|``, i.e. the number of
    atoms in the volume spanned by the new cell.  The first tolerance in
    ``tolerances`` that reproduces it is used.

    Args:
        atoms: ASE ``Atoms`` object.
        basis: ``(3, 3)`` scaled basis vectors, one per row.
        tolerances: Cut tolerances to try, smallest first.

    Returns:
        ``(atoms, tolerance)``.  ``tolerance`` is ``None`` when no value
        conserved the atom count; the closest result is returned and a
        warning is logged.
    """
    basis = np.asarray(basis, dtype=float)
    expected = int(round(len(atoms) * abs(float(np.linalg.det(basis)))))

    best = None
    best_tol = None
    for tol in tolerances:
        candidate = cut(atoms, *basis, tolerance=tol)
        if best is None or abs(len(candidate) - expected) < abs(len(best) - expected):
            best, best_tol = candidate, tol
        if len(candidate) == expected:
            return candidate, tol

    logging.warning(
        "cut() could not conserve the atom count for this Miller basis: "
        f"expected {expected} atoms, got {len(best)} (tolerance {best_tol}). "
        "The basis may not be commensurate with the cell."
    )
    return best, None


def layer_indices(
    coords: np.ndarray,
    cell: np.ndarray,
    gap: float = DEFAULT_LAYER_GAP,
) -> list[tuple[float, np.ndarray]]:
    """Split Cartesian coordinates into layers along the c-axis normal.

    Atoms are projected onto the unit normal of the ``a``-``b`` plane, the
    projections are wrapped into one interplanar period, and the sorted
    sequence is split wherever the gap exceeds ``gap``.  The group that
    straddles the periodic boundary is merged.

    Args:
        coords: ``(N, 3)`` Cartesian positions.
        cell: ``(3, 3)`` cell matrix whose rows are the lattice vectors.
        gap: Minimum inter-layer gap in Angstrom.

    Returns:
        List of ``(fractional_z, indices)`` sorted by ``fractional_z``.  The
        centre is the mean fractional coordinate along the c axis (the
        stacking direction), wrapped into ``[0, 1)``.

    Raises:
        ValueError: If the cell is degenerate or ``gap`` is not positive.
    """
    coords = np.asarray(coords, dtype=float)
    cell = np.asarray(cell, dtype=float)
    if gap <= 0:
        raise ValueError(f"Layer gap must be positive, got {gap}")
    if len(coords) == 0:
        return []

    a, b, c = cell
    cross = np.cross(a, b)
    area = float(np.linalg.norm(cross))
    if area < 1e-12:
        raise ValueError("Cannot identify layers: the a-b plane is degenerate")

    normal = cross / area
    period = abs(float(np.dot(c, normal)))
    if period < 1e-12:
        raise ValueError("Cannot identify layers: the c vector is parallel to the a-b plane")

    # Cluster along the physical stacking direction.
    proj = (coords @ normal) % period
    order = np.argsort(proj, kind="stable")
    sorted_proj = proj[order]

    splits = np.flatnonzero(np.diff(sorted_proj) > gap)
    groups = np.split(order, splits + 1)

    # Merge the groups that belong to the same layer across the periodic
    # boundary (first and last in sorted order).
    if len(groups) > 1:
        wrap_gap = period - sorted_proj[-1] + sorted_proj[0]
        if wrap_gap <= gap:
            groups[0] = np.concatenate([groups[-1], groups[0]])
            groups = groups[:-1]

    # Report the centre in fractional coordinates along c, which is the
    # coordinate system the callers (and the written POSCARs) work in.
    frac_z = coords @ np.linalg.inv(cell)
    frac_z = frac_z[:, 2]

    layers = []
    for indices in groups:
        values = frac_z[indices]
        # Unwrap a layer that straddles the periodic boundary before averaging.
        if np.ptp(values) > 0.5:
            values = np.where(values > 0.5, values - 1.0, values)
        centre = float(values.mean()) % 1.0
        # A layer centred exactly on the boundary can land on 1.0 because of
        # floating point; report it as 0.0 so centres stay in [0, 1).
        if centre > 1.0 - 1e-9:
            centre = 0.0
        layers.append((centre, indices))

    layers.sort(key=lambda item: item[0])
    return layers
