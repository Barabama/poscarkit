"""Reader for teacher-provided flat-table format (CSV and XLSX).

Column convention:
  - T: temperature
  - CO, CR, FE, NI, ...: site fractions for sublattice 1 (bare element names)
  - CO#2, CR#2, FE#2, ...: site fractions for sublattice 2
  - G, H, S, NP: thermodynamic properties (skipped)
"""

import re

import numpy as np
import pandas as pd

from src.io.ir import SOFData, get_site_ratios

_NON_SOF_COLS = frozenset(
    c.upper() for c in ("T", "G", "H", "S", "NP", "P", "n_kg", "n_mole")
)
_RE_ELEM = re.compile(r"^([A-Z][A-Z]?)(?:#(\d+))?$")


def read_teacher(path: str, phase_hint: str | None = None) -> SOFData:
    """Read teacher flat table (CSV or XLSX) into SOFData."""
    if path.endswith(".csv"):
        df = pd.read_csv(path)
    else:
        df = pd.read_excel(path)
    df.columns = [str(c).strip() for c in df.columns]

    phase = (phase_hint or "FCC").upper()
    site_ratios = get_site_ratios(phase)

    T = df["T"].values.astype(float)

    Y_subl, elements = _parse_Y_columns(df)

    composition = _derive_composition(Y_subl, site_ratios, elements)

    G_real = _find_G_column(df)
    H_real = _find_H_column(df)
    S_real = _find_S_column(df)

    return SOFData(
        source_path=path,
        phase=phase,
        site_ratios=site_ratios,
        T=T,
        Y_subl=Y_subl,
        composition=composition,
        elements=elements,
        G_real=G_real,
        H_real=H_real,
        S_real=S_real,
    )


def _parse_Y_columns(
    df: pd.DataFrame,
) -> tuple[list[dict[str, np.ndarray]], list[str]]:
    """Parse bare element columns into Y_subl."""
    subl_cols: dict[int, dict[str, str]] = {}
    elements: set[str] = set()

    for col in df.columns:
        upper = str(col).strip().upper()
        if upper in _NON_SOF_COLS:
            continue
        m = _RE_ELEM.match(str(col).strip())
        if not m:
            continue
        elem = m.group(1).upper()
        subl_idx = int(m.group(2)) - 1 if m.group(2) else 0
        elements.add(elem)
        subl_cols.setdefault(subl_idx, {})[elem] = col

    if not subl_cols:
        return [], []

    n_subl = max(subl_cols.keys()) + 1
    elem_list = sorted(elements)
    result: list[dict[str, np.ndarray]] = []
    for i in range(n_subl):
        d: dict[str, np.ndarray] = {}
        for elem in elem_list:
            if i in subl_cols and elem in subl_cols[i]:
                d[elem] = df[subl_cols[i][elem]].values.astype(float)
        result.append(d)

    return result, elem_list


def _derive_composition(
    Y_subl: list[dict[str, np.ndarray]],
    site_ratios: list[float],
    elements: list[str],
) -> dict[str, float]:
    """Derive overall composition from Y averages across T."""
    comp: dict[str, float] = {}
    for k, y_dict in enumerate(Y_subl):
        r = site_ratios[k] if k < len(site_ratios) else 0.0
        for elem in elements:
            if elem in y_dict:
                avg = float(np.mean(y_dict[elem]))
                comp[elem] = comp.get(elem, 0.0) + r * avg
    total = sum(comp.values())
    if total > 0:
        comp = {e: v / total for e, v in comp.items()}
    return comp


def _find_G_column(df: pd.DataFrame) -> np.ndarray | None:
    """Locate Gibbs free energy column: G, G(), G(@), gibbs, deltaG."""
    for col in df.columns:
        col_str = str(col).strip()
        if col_str.upper() == "G":
            return df[col].values.astype(float)
        if col_str.upper().startswith("G("):
            return df[col].values.astype(float)
        if col_str.upper().startswith("G(@"):
            return df[col].values.astype(float)
        if "gibbs" in col_str.lower() or "deltaG" in col_str.lower():
            return df[col].values.astype(float)
    return None


def _find_H_column(df: pd.DataFrame) -> np.ndarray | None:
    """Locate enthalpy column: H, H(), H(@), enthalpy."""
    for col in df.columns:
        col_str = str(col).strip()
        if col_str.upper() in ("H", "H(T)"):
            return df[col].values.astype(float)
        if col_str.upper().startswith("H("):
            return df[col].values.astype(float)
        if col_str.upper().startswith("H(@"):
            return df[col].values.astype(float)
        if "enthalpy" in col_str.lower() or "deltaH" in col_str.lower():
            return df[col].values.astype(float)
    return None


def _find_S_column(df: pd.DataFrame) -> np.ndarray | None:
    """Locate entropy column: S, S(), S(@), entropy."""
    for col in df.columns:
        col_str = str(col).strip()
        if col_str.upper() in ("S", "S(T)"):
            return df[col].values.astype(float)
        if col_str.upper().startswith("S("):
            return df[col].values.astype(float)
        if col_str.upper().startswith("S(@"):
            return df[col].values.astype(float)
        if "entropy" in col_str.lower() or "deltaS" in col_str.lower():
            return df[col].values.astype(float)
    return None
