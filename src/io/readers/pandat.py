"""Reader for Pandat flat-table format (CSV and XLSX).

Column conventions:
  - T: temperature
  - phase_name: phase identifier
  - x(ELEM): overall mole fraction
  - y(ELEM#N@PHASE): site fraction of ELEM on sublattice N of PHASE
  - G(@PHASE): phase Gibbs energy
"""

import re

import numpy as np
import pandas as pd

from src.io.ir import SOFData, get_site_ratios


def read_pandat(path: str, phase_hint: str | None = None) -> SOFData:
    """Read Pandat flat table (CSV or XLSX) into SOFData."""
    if path.endswith(".csv"):
        df = pd.read_csv(path)
    else:
        df = pd.read_excel(path)
    df.columns = [str(c).strip() for c in df.columns]

    phase = _detect_phase(df)
    if phase_hint:
        phase = phase_hint
    phase = phase.upper()

    site_ratios = get_site_ratios(phase)

    T = df["T"].values.astype(float)

    composition = _extract_composition(df)
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


def _detect_phase(df: pd.DataFrame) -> str:
    """Detect phase from phase_name column or y(...@PHASE) patterns."""
    if "phase_name" in df.columns:
        val = df["phase_name"].dropna().iloc[0]
        return str(val).strip().upper()
    for col in df.columns:
        m = re.search(r"@(\w+)\)", str(col), re.IGNORECASE)
        if m:
            return m.group(1).upper()
    return "FCC"


def _extract_composition(df: pd.DataFrame) -> dict[str, float]:
    """Extract nominal composition from x(ELEM) columns.

    Handles both mole fractions and percentages automatically.
    """
    comp: dict[str, float] = {}
    for col in df.columns:
        m = re.match(r"x\((\w+)\)", str(col), re.IGNORECASE)
        if m:
            elem = m.group(1).upper()
            val = float(df[col].iloc[0])
            if val > 1e-12:
                comp[elem] = val
    if not comp:
        raise ValueError("No x(ELEM) columns found in Pandat data")

    total = sum(comp.values())
    if abs(total - 100.0) < 50.0:
        comp = {e: v / 100.0 for e, v in comp.items()}
    elif total > 1.5:
        comp = {e: v / total for e, v in comp.items()}
    return comp


def _extract_Y_columns(
    df: pd.DataFrame, phase: str, elements: list[str]
) -> list[dict[str, np.ndarray]]:
    """Extract site fractions from y(ELEM#N@PHASE) columns."""
    subl_cols: dict[int, dict[str, str]] = {}

    for col in df.columns:
        m = re.match(
            rf"y\((\w+)#(\d+)@{phase}\)", str(col), re.IGNORECASE
        )
        if not m:
            continue
        elem = m.group(1).upper()
        if elem not in elements:
            continue
        subl_idx = int(m.group(2)) - 1
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
    """Extract Gibbs energy from G(@PHASE) column."""
    col_name = f"G(@{phase})"
    if col_name in df.columns:
        return df[col_name].values.astype(float)
    for col in df.columns:
        lower = str(col).lower()
        if lower == "g":
            return df[col].values.astype(float)
    return None
