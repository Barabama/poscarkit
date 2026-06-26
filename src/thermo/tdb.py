"""TDB file parser — self-contained, no dependency on GRS or em-tdb.

TDB (Thermodynamic Database) format:
  - Records are delimited by !
  - Comments run from $ to end of line
  - FUNCTION defines temperature-dependent functions (e.g. SER pure-element references)
  - PARAMETER defines end-member Gibbs free energies
  - PHASE / CONSTITUENT define sublattice structure per phase

This module extracts FUNCTION, PARAMETER (G type, order=0), and PHASE records
needed to build Gibbs free energy expressions.
"""

import re
import warnings


def parse_tdb(tdb_path: str, phase: str | None = None) -> dict:
    """Parse a TDB file, returning functions, phase parameters, and phase structure.

    When phase is None, extracts end-member parameters for ALL phases.
    When phase is given, filters to that phase only.

    Returns:
        dict with keys:
            functions: dict[str, str]       — {func_name: expression}
            params: dict[str, str]          — {"PHASE,ELEM1:ELEM2": expression}
            phase_ratios: dict[str, list[float]] — {phase_name: [ratio, ...]}
    """
    with open(tdb_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Strip comments and join (multi-line records become single-line)
    text = "".join(_strip_comment(line) for line in lines)

    functions: dict[str, str] = {}
    params: dict[str, str] = {}
    phase_info: dict[str, dict] = {}  # {phase_name: {"count": int, "ratios": [float]}}

    for record in text.split("!"):
        record = record.strip()
        if not record:
            continue

        upper = record.upper()

        if upper.startswith("FUNCTION"):
            m = _parse_function(record)
            if not m:
                continue
            name, expr = m
            if name in functions:
                warnings.warn(
                    f"TDB: FUNCTION '{name}' has multiple temperature ranges — "
                    f"keeping first expression only. DeltaG may be inaccurate "
                    f"outside the first range's validity interval."
                )
            else:
                functions[name] = expr

        elif upper.startswith("PARAMETER"):
            if m := _parse_parameter(record, phase):
                key, expr = m
                params[key] = expr

        elif upper.startswith("PHASE "):
            if m := _parse_phase(record):
                p_name, ratios = m
                phase_info[p_name.upper()] = ratios

    phase_ratios = {p: r["ratios"] for p, r in phase_info.items()}
    return {
        "functions": functions,
        "params": params,
        "phase_ratios": phase_ratios,
    }


def resolve_expression(expr: str, functions: dict[str, str]) -> str:
    """Replace all SERXX# / GHSERXX# / ETOT_SER_XX# references with function bodies.

    "0.25*SERCO#" → "0.25*((-6.815E+05+...))"
    "GHSERFE#"    → "((-6.815E+05+...))"
    "ETOT_SER_CO#" → "((-6.638E+05+...))"
    Function bodies are wrapped in parentheses to preserve operator precedence.
    """

    def _replace_ser(m: re.Match) -> str:
        ser_name = m.group(0)[:-1]  # strip trailing #
        if ser_name in functions:
            return f"({functions[ser_name]})"
        # Some TDBs name functions GHSERXX but params reference SERXX#
        ghser_name = "GH" + ser_name
        if ghser_name in functions:
            return f"({functions[ghser_name]})"
        raise KeyError(f"Function not found in TDB: {ser_name} (also tried {ghser_name})")

    # Match SERXX#, GHSERXX#, and ETOT_SER_XX# (case-insensitive element suffix)
    expr = re.sub(r"(?:(?:GH)?SER[A-Za-z]+|ETOT_SER_[A-Za-z]+)#", _replace_ser, expr)
    return expr


def _strip_comment(line: str) -> str:
    """Strip $ comment from a line"""
    return line.split("$", 1)[0].strip()


def _parse_function(record: str) -> tuple[str, str] | None:
    """Parse a FUNCTION record.

    FUNCTION SERCO   1.00 -6.815E+05+7.178E+01*T-...; 6000.00 N !
    """
    m = re.match(
        r"FUNCTION\s+(\S+)\s+\S+\s+(.+);\s*\S+\s*\S*",
        record,
        re.IGNORECASE,
    )
    if not m:
        return None
    name = m.group(1)
    expr = m.group(2).strip()
    expr = expr.replace("LN(", "ln(")  # sympy uses lowercase ln
    return name, expr


def _parse_parameter(record: str, target_phase: str | None) -> tuple[str, str] | None:
    """Parse a PARAMETER record, keeping only G parameters (order=0).

    When target_phase is None, extracts all phases.
    Otherwise filters to the target phase only.

    PARAMETER G(FCC,CO:CR;0) 1.00 ...expr...; 6000.00 N REF1 !
    """
    # Normalize inconsistent whitespace: some TDBs have spaces around
    # ":" and ";" in parameter signatures (e.g. G(FCC,V :AL;0)).
    # Removing them avoids fragile regex patterns.
    record_norm = re.sub(r"\s*([:;])", r"\1", record)

    m = re.match(
        r"PARAMETER\s+G\((\S+),(\S+);(\d+)\)\s+\S+\s+(.+);\s*\S+",
        record_norm,
        re.IGNORECASE,
    )
    if not m:
        return None
    p_phase = m.group(1)
    components = m.group(2)
    order_num = int(m.group(3))
    expr = m.group(4).strip()

    if target_phase is not None and p_phase.upper() != target_phase.upper():
        return None
    if order_num != 0:
        return None  # end-member only, skip interaction parameters

    expr = expr.replace("LN(", "ln(")
    key = f"{p_phase.upper()},{components}"
    return key, expr


def _parse_phase(record: str) -> tuple[str, dict] | None:
    """Parse a PHASE record to extract sublattice stoichiometry.

    PHASE FCC % 2 0.25 0.75 !
    PHASE BCC % 2 0.5 0.5 !
    Returns: (phase_name, {"ratios": [float, ...]})
    """
    m = re.match(
        r"PHASE\s+(\S+)\s+%\s+(\d+)\s+([\d.\s]+)",
        record,
        re.IGNORECASE,
    )
    if not m:
        return None
    phase_name = m.group(1)
    ratios = [float(x) for x in m.group(3).split()]
    if abs(sum(ratios) - 1.0) > 1e-4:
        return None  # stoichiometry line, not ratio line
    return phase_name.upper(), {"ratios": ratios}
