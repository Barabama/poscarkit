"""Plot output — 2×2 subplot (ΔG_f, ΔH_f, ΔS_f, S_conf)"""

from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from src.io.ir import SOFData


def draw_figure(
    ir: SOFData,
    G_random: np.ndarray,
    Sconf_real: np.ndarray,
    Sconf_random: float,
    H_random: np.ndarray,
    S_random: np.ndarray,
    figsize: tuple[float, float] = (14, 10),
    dpi: int = 100,
) -> Figure:
    """Draw the 2×2 figure: ΔG_f, ΔH_f, ΔS_f, S_conf."""
    has_G = ir.G_real is not None
    has_H = ir.H_real is not None
    has_S = ir.S_real is not None

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize, dpi=dpi)

    T = ir.T

    # === ΔG_f (top-left) ===
    if has_G:
        ax1.plot(T, ir.G_real, "b-", lw=1.5, label="ΔG_f (real)")
    ax1.plot(T, G_random, "r--", lw=1.5, label="ΔG_f (random)")
    ax1.set_xlabel("Temperature (K)")
    ax1.set_ylabel("ΔG_f (J/mol-atom)")
    ax1.set_title("ΔG_f")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # === ΔH_f (top-right) ===
    if has_H:
        ax2.plot(T, ir.H_real, "b-", lw=1.5, label="ΔH_f (real)")
    ax2.plot(T, H_random, "r--", lw=1.5, label="ΔH_f (random)")
    ax2.set_xlabel("Temperature (K)")
    ax2.set_ylabel("ΔH_f (J/mol-atom)")
    ax2.set_title("ΔH_f")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    # === ΔS_f (bottom-left) ===
    if has_S:
        ax3.plot(T, ir.S_real, "b-", lw=1.5, label="ΔS_f (real)")
    ax3.plot(T, S_random, "r--", lw=1.5, label="ΔS_f (random)")
    ax3.set_xlabel("Temperature (K)")
    ax3.set_ylabel("ΔS_f (J/mol-atom-K)")
    ax3.set_title("ΔS_f")
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)

    # === S_conf (bottom-right) ===
    ax4.plot(T, Sconf_real, "b-", lw=1.5, label="S_conf (real SOFs)")
    ax4.axhline(y=Sconf_random, color="r", ls="--", lw=1.5,
                label=f"S_conf (random) = {Sconf_random:.2f}")
    ax4.set_xlabel("Temperature (K)")
    ax4.set_ylabel("S_conf (J/mol-atom-K)")
    ax4.set_title("S_conf")
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)

    comp_parts = [f"{e}{v:.3g}" for e, v in sorted(ir.composition.items())]
    comp_str = "-".join(comp_parts)
    fig.suptitle(f"Thermodynamic Properties — {comp_str} ({ir.phase})", fontweight="bold")
    fig.tight_layout()
    return fig


def save_plot(
    ir: SOFData,
    G_random: np.ndarray,
    Sconf_real: np.ndarray,
    Sconf_random: float,
    H_random: np.ndarray,
    S_random: np.ndarray,
    outdir: Path,
) -> Path:
    """Save 2×2 chart as a PNG file."""
    fig = draw_figure(ir, G_random, Sconf_real, Sconf_random, H_random, S_random, dpi=150)
    fname_parts = [f"{e}{v:.3g}".replace(".", "_") for e, v in sorted(ir.composition.items())]
    fname = f"Thermo_{'-'.join(fname_parts)}.png"
    fpath = outdir / fname
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fpath
