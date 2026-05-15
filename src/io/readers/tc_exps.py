"""Reader for Thermo-Calc tc_exps format (single-sheet and multi-sheet XLSX).

Column convention (single-sheet or after merge):
  - T: temperature
  - Y(PHASE,EL): site fraction on sublattice 1
  - Y(PHASE,EL#N): site fraction on sublattice N
  - X(PHASE,EL): overall mole fraction
  - G(PHASE): phase Gibbs energy
"""

import re
import numpy as np
import pandas as pd

from src.io.ir import SOFData, PHASE_SITE_RATIOS


def read_tc_exps(path: str, phase_hint: str | None = None) -> SOFData:
    """Read Thermo-Calc exported data (XLSX or CSV) into SOFData.

    Handles both:
      - Single-sheet XLSX/CSV with all columns in one table.
      - Multi-sheet XLSX (merge-on-T handled by detect_format → re-read path).
    """
    if path.endswith(".csv"):
        df = pd.read_csv(path)
    else:
        df = pd.read_excel(path)
    df.columns = [str(c).strip() for c in df.columns]

    # Detect phase from column headers (e.g. "G(FCC)" or "Y(FCC,CO)")
    phase = phase_hint or _detect_phase_from_columns(df)
    phase = phase.upper()

    site_ratios = PHASE_SITE_RATIOS.get(phase, [0.5, 0.5])

    T = df["T"].values.astype(float)

    composition = _extract_composition(df, phase)
    elements = sorted(composition.keys())

    Y_subl = _extract_Y_columns(df, phase, elements)

    G_real = _extract_G_real(df, phase)

    return SOFData(
        source_path=path,
        phase=phase,
        site_ratios=site_ratios,
        T=T,
        Y_subl=Y_subl,
        composition=composition,
        elements=elements,
        G_real=G_real,
    )


def _detect_phase_from_columns(df: pd.DataFrame) -> str:
    """Extract phase name from G(PHASE) or H(PHASE) column headers."""
    for col in df.columns:
        m = re.match(r"[GH]\((\w+)\)", str(col), re.IGNORECASE)
        if m:
            return m.group(1).upper()
    return "FCC"


def _extract_composition(df: pd.DataFrame, phase: str) -> dict[str, float]:
    """Extract nominal composition from X(PHASE,EL) columns."""
    comp: dict[str, float] = {}
    pattern = re.compile(rf"X\({phase},\s*(\w+)\)", re.IGNORECASE)
    for col in df.columns:
        m = pattern.match(str(col))
        if m:
            elem = m.group(1).upper()
            val = float(df[col].iloc[0])
            if val > 1e-12:
                comp[elem] = val
    if not comp:
        raise ValueError(
            f"No X({phase},ELEM) columns found in TC data"
        )
    return comp


def _extract_Y_columns(
    df: pd.DataFrame, phase: str, elements: list[str]
) -> list[dict[str, np.ndarray]]:
    """Extract Y site fractions from Y(PHASE,EL#N) columns."""
    subl_cols: dict[int, dict[str, str]] = {}

    for col in df.columns:
        m = re.match(
            rf"Y\({phase},\s*(\w+)(?:#(\d+))?\)", str(col), re.IGNORECASE
        )
        if not m:
            continue
        elem = m.group(1).upper()
        if elem not in elements:
            continue
        if m.group(2):
            subl_idx = int(m.group(2)) - 1  # #2 → 1, #3 → 2
        else:
            subl_idx = 0
        subl_cols.setdefault(subl_idx, {})[elem] = col

    n_subl = max(subl_cols.keys()) + 1 if subl_cols else 0
    result: list[dict[str, np.ndarray]] = []
    for i in range(n_subl):
        d: dict[str, np.ndarray] = {}
        for elem in elements:
            if i in subl_cols and elem in subl_cols[i]:
                d[elem] = df[subl_cols[i][elem]].values.astype(float)
        result.append(d)

    return result


def _extract_G_real(df: pd.DataFrame, phase: str) -> np.ndarray | None:
    """Extract Gibbs energy from G(PHASE) column."""
    col_name = f"G({phase})"
    if col_name in df.columns:
        return df[col_name].values.astype(float)
    for col in df.columns:
        m = re.match(rf"G\({phase}\)", str(col), re.IGNORECASE)
        if m:
            return df[col].values.astype(float)
    return None
