

import re

def _parse_expression(line: str):
    digital = r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
    pattern_map = [
        ("A", fr"(?<!\*){digital}(?=\+|\-)", lambda m: f"{float(m.group(1)):+E}"),
        ("B", fr"{digital}\*T(?!\*)", lambda m: f"{float(m.group(1)):+E}*T"),
        ("C", fr"{digital}\*T\*LN\(T\)", lambda m: f"{float(m.group(1)):+E}*T*LN(T)"),
        ("D", fr"{digital}\*T\*\*2", lambda m: f"{float(m.group(1)):+E}*T**2"),
        ("E", fr"{digital}\*T\*\*3", lambda m: f"{float(m.group(1)):+E}*T**3"),
        ("F", fr"{digital}\*T\*\*\(-1\)", lambda m: f"{float(m.group(1)):+E}*T**(-1)"),
        ("X", fr"{digital}\*([^#\+\-][A-Za-z]+#)", lambda m: f"{float(m.group(1))}*{m.group(2)}"),
    ]
    new_line = ""
    for key, pattern, formatter in pattern_map:
        for match in re.finditer(pattern, line):
            new_line += formatter(match)
    
    print(f"old_line={line}, new_line={new_line}")
    return new_line



old_line = '-681500-14.54*T*LN(T)-0.009912*T**2+0.1033E-07*T**3-3.081E+04*T**(-1)+71.78*T-0.5*SERCO#-5.5*T**3'

print(_parse_expression(old_line))
