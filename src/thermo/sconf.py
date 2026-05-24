"""Configurational entropy — phase-agnostic, N-sublattice formula.

Sconf = -R * sum_k [ ratio_k * sum_i (y_i^(k) * ln(y_i^(k))) ]

Phase site ratios are resolved via ir.get_site_ratios() or PHASE_SITE_RATIOS.
"""

import numpy as np

R = 8.314  # J/(mol*K)


def calc_sconf(Y_subl: list[dict[str, float]], site_ratios: list[float]) -> float:
    """Configurational entropy for a phase with arbitrary sublattice count.

    Args:
        Y_subl: site fractions per sublattice, Y_subl[0] = subl1, Y_subl[1] = subl2, ...
        site_ratios: stoichiometric ratio per sublattice (must sum to 1).

    Returns:
        Sconf in J/(mol*K)
    """
    total = 0.0
    for sofs, ratio in zip(Y_subl, site_ratios):
        sub = 0.0
        for y in sofs.values():
            if y > 1e-30:
                sub += y * np.log(y)
        total += ratio * sub
    return -R * total


def calc_sconf_random(composition: dict[str, float], site_ratios: list[float]) -> float:
    """Configurational entropy for random mixing.

    Nominal composition is used as the site fraction on every sublattice:
      Sconf_random = -R * (sum_k ratio_k) * sum_i (x_i * ln(x_i))
    Since ratios sum to 1:
      Sconf_random = -R * sum_i (x_i * ln(x_i))
    """
    # Same composition on all sublattices → sum of ratios = 1 → single-lattice equivalent
    sub = sum(x * np.log(x) for x in composition.values() if x > 1e-30)
    total_ratio = sum(site_ratios)
    return -R * total_ratio * sub
