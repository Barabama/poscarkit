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

from src.io.ir import SOFData, get_site_ratios

_RE_COL_HAS_PHASE = re.compile(r"[YGHX]\((\w+)", re.IGNORECASE)


def read_tc_exps(path: str, phase_hint: str | None = None) -> SOFData:
    """Read Thermo-Calc exported data (XLSX or CSV) into SOFData.

    Handles:
      - Multi-sheet XLSX: finds Y, G, X sheets by column patterns, aligns T.
      - Single-sheet XLSX/CSV: all columns in one table.
    """
    if path.endswith(".csv"):
        df = pd.read_csv(path)
        df.columns = [str(c).strip() for c in df.columns]
        df_y = df_g = df_x = df
    else:
        xl = pd.ExcelFile(path)
        if len(xl.sheet_names) > 1:
            df_y = _find_and_read(xl, path, "Y(")
            df_g = _find_and_read(xl, path, "G(")
            df_x = _find_and_read(xl, path, "X(")
        else:
            df = pd.read_excel(path)
            df.columns = [str(c).strip() for c in df.columns]
            df_y = df_g = df_x = df

    phase = phase_hint or _detect_phase_from_columns(df_g)
    phase = phase.upper()

    site_ratios = get_site_ratios(phase)

    T = _align_and_extract_T(df_y, df_g)

    composition = _extract_composition_from_x(df_x, phase)
    elements = sorted(composition.keys())

    Y_subl = _extract_Y_columns(df_y, phase, elements, T)

    G_real = _extract_G_real(df_g, phase, T)

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


def _find_and_read(xl: pd.ExcelFile, path: str, pattern: str) -> pd.DataFrame:
    """Find a sheet whose columns contain the given pattern, and read it."""
    for s in xl.sheet_names:
        df = pd.read_excel(xl, sheet_name=s, nrows=0)
        cols = [str(c).strip() for c in df.columns]
        if any(pattern in c for c in cols):
            result = pd.read_excel(path, sheet_name=s)
            result.columns = [str(c).strip() for c in result.columns]
            return result
    # Fallback: return first sheet
    df = pd.read_excel(path)
    df.columns = [str(c).strip() for c in df.columns]
    return df


def _find_T_col(df: pd.DataFrame) -> str:
    for col in df.columns:
        if str(col).strip().upper() == "T":
            return col
    return df.columns[0]


def _T_indices(df: pd.DataFrame, T_common: np.ndarray) -> np.ndarray:
    t_col = _find_T_col(df)
    return np.where(np.isin(df[t_col].values, T_common))[0]


def _align_and_extract_T(
    df_y: pd.DataFrame, df_g: pd.DataFrame
) -> np.ndarray:
    """Align T values between Y and G dataframes."""
    t_y = df_y[_find_T_col(df_y)].values
    t_g = df_g[_find_T_col(df_g)].values
    common = np.intersect1d(np.round(t_y, 6), np.round(t_g, 6))
    return np.sort(common)


def _detect_phase_from_columns(df: pd.DataFrame) -> str:
    """Extract phase name from G(PHASE) or H(PHASE) column headers."""
    for col in df.columns:
        m = re.match(r"[GH]\((\w+)\)", str(col), re.IGNORECASE)
        if m:
            return m.group(1).upper()
    return "FCC"


def _extract_composition_from_x(
    df_x: pd.DataFrame, phase: str
) -> dict[str, float]:
    """Extract nominal composition from X(PHASE,EL) columns."""
    comp: dict[str, float] = {}
    pattern = re.compile(rf"X\({phase},\s*(\w+)\)", re.IGNORECASE)
    for col in df_x.columns:
        m = pattern.match(str(col))
        if m:
            elem = m.group(1).upper()
            val = float(df_x[col].iloc[0])
            if val > 1e-12:
                comp[elem] = val
    if not comp:
        raise ValueError(f"No X({phase},ELEM) columns found in TC data")
    return comp


def _extract_Y_columns(
    df_y: pd.DataFrame,
    phase: str,
    elements: list[str],
    T: np.ndarray,
) -> list[dict[str, np.ndarray]]:
    """Extract Y site fractions from Y(PHASE,EL#N) columns."""
    idx_y = _T_indices(df_y, T)
    subl_cols: dict[int, dict[str, str]] = {}

    for col in df_y.columns:
        m = re.match(
            rf"Y\({phase},\s*(\w+)(?:#(\d+))?\)", str(col), re.IGNORECASE
        )
        if not m:
            continue
        elem = m.group(1).upper()
        if elem not in elements:
            continue
        if m.group(2):
            subl_idx = int(m.group(2)) - 1
        else:
            subl_idx = 0
        subl_cols.setdefault(subl_idx, {})[elem] = col

    n_subl = max(subl_cols.keys()) + 1 if subl_cols else 0
    result: list[dict[str, np.ndarray]] = []
    for i in range(n_subl):
        d: dict[str, np.ndarray] = {}
        for elem in elements:
            if i in subl_cols and elem in subl_cols[i]:
                d[elem] = df_y[subl_cols[i][elem]].values[idx_y]
        result.append(d)

    return result


def _extract_G_real(
    df_g: pd.DataFrame, phase: str, T: np.ndarray
) -> np.ndarray | None:
    """Extract Gibbs energy from G(PHASE) column, aligned to T."""
    col_name = f"G({phase})"
    target_col = None
    if col_name in df_g.columns:
        target_col = col_name
    else:
        for col in df_g.columns:
            if re.match(rf"G\({phase}\)", str(col), re.IGNORECASE):
                target_col = col
                break
    if target_col is None:
        return None
    idx_g = _T_indices(df_g, T)
    return df_g[target_col].values[idx_g]
