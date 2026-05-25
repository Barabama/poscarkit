"""DeltaG_random calculation — TDB end-member parameters + nominal composition.

G_random(T) = sum_{combo in elements^N} (prod(x_i) * G_combo(T)) + R * T * Sconf_random

where N = number of sublattices, and G_combo(T) is the end-member expression
G(phase, elem_1:elem_2:...:elem_N; 0) from the TDB.
"""

from itertools import product

import numpy as np

from src.thermo.tdb import resolve_expression

R = 8.314  # J/(mol*K)


def calc_deltaG_random(
    composition: dict[str, float],
    functions: dict[str, str],
    params: dict[str, str],
    T_values: np.ndarray,
    sconf_random: float,
    phase: str,
    site_ratios: list[float],
) -> np.ndarray:
    """Calculate DeltaG_random(T) for a phase with N sublattices.

    Args:
        composition: nominal composition {ELEMENT: fraction}
        functions: TDB FUNCTION dict {name: expression}
        params: end-member parameter dict {"PHASE,ELEM1:ELEM2:...": expression}
        T_values: temperature points (K)
        sconf_random: random configurational entropy (J/(mol*K))
        phase: phase name (e.g. "FCC")
        site_ratios: sublattice ratios (determines N = len(site_ratios))

    Returns:
        DeltaG_random array (J/mol)
    """
    elements = sorted(composition.keys())
    n_subl = len(site_ratios)
    terms: list[str] = []

    for combo in product(elements, repeat=n_subl):
        weight = 1.0
        for elem in combo:
            weight *= composition[elem]
        if weight < 1e-30:
            continue

        key = f"{phase},{':'.join(combo)}"
        if key not in params:
            raise KeyError(
                f"End-member parameter not found in TDB: {key}. "
                f"Ensure the TDB file contains all end-member definitions "
                f"for the {phase} phase with {n_subl} sublattices."
            )

        raw_expr = params[key]
        resolved_expr = resolve_expression(raw_expr, functions)
        terms.append(f"{weight:.15e}*({resolved_expr})")

    # Entropy contribution: R * T * Sconf_random
    terms.append(f"{R:.15e}*T*{sconf_random:.15e}")

    expr_str = "+".join(terms)
    return _eval_tdb_expr(expr_str, T_values)


def _eval_tdb_expr(expr_str: str, T: np.ndarray) -> np.ndarray:
    """Evaluate a resolved TDB expression string at temperatures T.

    Drops dependencies on sympy (~25 MB). Since resolve_expression() already
    replaces all SER references, the result is a pure arithmetic string in T
    with LN and EXP calls.  We replace LN/EXP with numpy log/exp and use a
    restricted eval — the namespace contains only numpy arrays and math
    functions, with no access to builtins.
    """
    expr = expr_str.replace("LN(", "log(").replace("EXP(", "exp(")
    ns = {"T": T, "log": np.log, "exp": np.exp, "__builtins__": {}}
    return eval(expr, ns)
