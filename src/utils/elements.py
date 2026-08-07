# src/utils/elements.py

"""Chemical element symbol canonicalization helpers."""


def canonicalize_element(symbol: str) -> str:
    """Return the canonical case of a chemical element symbol.

    Thermo-Calc / Pandat export files write element symbols in all caps
    (``CO``, ``CU``, ``NI`` for cobalt, copper, nickel).  ASE, however,
    requires the canonical chemical symbols (``Co``, ``Cu``, ``Ni``) and
    raises ``KeyError`` for anything else.  Every chemical element symbol
    is one or two letters with only the first letter capitalized, so
    title-casing recovers the canonical form from all-caps or mixed-case
    input.

    Args:
        symbol: An element symbol in any case (e.g. ``"CO"``, ``"co"``,
            ``"Co"``).

    Returns:
        The canonical symbol (e.g. ``"Co"``).  Empty strings pass through
        unchanged.
    """
    if not symbol:
        return symbol
    return symbol[0].upper() + symbol[1:].lower()
