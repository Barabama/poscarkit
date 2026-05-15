"""Format auto-detection and reader dispatch."""

import os
import re
import pandas as pd

from src.io.ir import SOFData

_RE_PANDAT = re.compile(r"y\([^)]*@[^)]*\)", re.IGNORECASE)
_RE_TC_EXPS = re.compile(r"[XY]\([^)]*,[^)]*\)")


def _sheet_has_T(xl: pd.ExcelFile, sheet: str) -> bool:
    """Check if a sheet has a 'T' column in its header."""
    df = pd.read_excel(xl, sheet_name=sheet, nrows=0)
    cols = [str(c).strip().upper() for c in df.columns]
    return "T" in cols


def _merge_sheets_on_T(xl: pd.ExcelFile) -> pd.DataFrame:
    """Merge all sheets on the T column (inner join).

    Deduplicates columns that appear in multiple sheets.
    """
    dfs = []
    for sheet in xl.sheet_names:
        df = pd.read_excel(xl, sheet_name=sheet)
        df.columns = [str(c).strip() for c in df.columns]
        dfs.append(df)

    merged = dfs[0]
    for df in dfs[1:]:
        on_cols = [c for c in df.columns if c in merged.columns]
        if not on_cols:
            merged = pd.concat([merged, df], axis=1)
        else:
            merged = merged.merge(df, on=on_cols, how="inner", suffixes=("", "_dup"))
        dup_cols = [c for c in merged.columns if c.endswith("_dup")]
        merged = merged.drop(columns=dup_cols)
    return merged


def _detect_from_columns(cols: list[str]) -> str:
    """Detect format from column name patterns."""
    for c in cols:
        if _RE_PANDAT.search(str(c)):
            return "pandat"
    for c in cols:
        if _RE_TC_EXPS.search(str(c)):
            return "tc_exps"
    return "teacher"


def detect_format(path: str) -> str:
    """Detect upstream format: 'tc_exps' | 'pandat' | 'teacher'.

    Logic:
      1. Multi-sheet xlsx where every sheet has a T column → merge on T → re-detect.
      2. Column pattern: y(...@...) → pandat, Y(...,...) → tc_exps, else teacher.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"File not found: {path}")

    if path.endswith((".xlsx", ".xls")):
        xl = pd.ExcelFile(path)
        if len(xl.sheet_names) > 1:
            all_have_T = all(_sheet_has_T(xl, s) for s in xl.sheet_names)
            if all_have_T:
                df = _merge_sheets_on_T(xl)
                return _detect_from_columns([str(c) for c in df.columns])
        df = pd.read_excel(path, nrows=0)
    else:
        df = pd.read_csv(path, nrows=0, dtype=str)

    return _detect_from_columns([str(c) for c in df.columns])


def read(path: str, phase_hint: str | None = None) -> SOFData:
    """Auto-detect format and read into SOFData.

    Args:
        path: Path to xlsx/xls/csv file.
        phase_hint: Optional phase override (e.g. 'FCC' for teacher format).

    Returns:
        SOFData dataclass.
    """
    fmt = detect_format(path)

    if fmt == "tc_exps":
        from src.io.readers.tc_exps import read_tc_exps
        return read_tc_exps(path, phase_hint=phase_hint)

    elif fmt == "pandat":
        from src.io.readers.pandat import read_pandat
        return read_pandat(path, phase_hint=phase_hint)

    elif fmt == "teacher":
        from src.io.readers.teacher import read_teacher
        return read_teacher(path, phase_hint=phase_hint)

    raise ValueError(f"Unknown format: {fmt}")
