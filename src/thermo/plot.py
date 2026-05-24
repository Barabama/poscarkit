"""Plot output — dual-panel PNG (Sconf-T + DeltaG-T)"""

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
    S_real: np.ndarray,
    S_random: float,
    figsize: tuple[float, float] = (14, 5),
    dpi: int = 100,
) -> Figure:
    """Draw the dual-panel Sconf-T + DeltaG-T figure.

    Can be embedded in a tkinter canvas or exported to PNG.
    Returns the Figure object.
    """
    has_G_real = ir.G_real is not None

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, dpi=dpi)

    # Left: Sconf-T
    ax1.plot(ir.T, S_real, "b-", lw=1.5, label="Sconf (real SOFs)")
    ax1.axhline(
        y=S_random, color="r", ls="--", lw=1.5, label=f"Sconf (random) = {S_random:.2f}"
    )
    ax1.set_xlabel("Temperature (K)")
    ax1.set_ylabel("Sconf (J/mol-K)")
    ax1.set_title("Configurational Entropy")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Right: DeltaG-T
    if has_G_real:
        ax2.plot(ir.T, ir.G_real, "b-", lw=1.5, label="ΔG (real SOFs)")
    ax2.plot(ir.T, G_random, "r--", lw=1.5, label="ΔG (random)")
    ax2.set_xlabel("Temperature (K)")
    ax2.set_ylabel("ΔG (J/mol)")
    ax2.set_title("Gibbs Free Energy")
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    comp_parts = [f"{e}{v:.3g}" for e, v in sorted(ir.composition.items())]
    comp_str = "-".join(comp_parts)
    fig.suptitle(f"Sconf & ΔG — {comp_str} ({ir.phase})", fontweight="bold")
    fig.tight_layout()
    return fig


def save_plot(
    ir: SOFData,
    G_random: np.ndarray,
    S_real: np.ndarray,
    S_random: float,
    outdir: Path,
) -> Path:
    """Save dual-panel chart as a PNG file alongside the source file.

    Returns:
        Output PNG file path
    """
    fig = draw_figure(ir, G_random, S_real, S_random, dpi=150)
    # out_dir = Path(ir.source_path).parent
    fname_parts = [f"{e}{v:.3g}".replace(".", "_") for e, v in sorted(ir.composition.items())]
    fname = f"Sconf_DeltaG_{'-'.join(fname_parts)}.png"
    fpath = outdir / fname
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return fpath
