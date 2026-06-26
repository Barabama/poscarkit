# TDBManager.py

import re
from dataclasses import dataclass, asdict
from typing import TypedDict


class Elem(TypedDict):
    elem: str
    ref_state: str
    atomic_mass: float
    h298_h0: float
    s298: float


class Func(TypedDict):
    func: str
    elem: str
    temp_start: float
    temp_end: float
    expression: str
    is_continued: str


class Tdb(TypedDict):
    tdb: str
    description: str
    version: str


class Phase(TypedDict):
    phase: str
    sub_lattices: int
    stoichiometry: str
    components: str
    tdb: str


class Param(TypedDict):
    param: str
    ptype: str
    phase: str
    components: str
    order_num: int
    temp_start: float
    temp_end: float
    expression: str
    is_continued: str
    tdb: str


@dataclass
class ParsedData:
    elems: list[Elem]
    funcs: list[Func]
    phases: list[Phase]
    params: list[Param]
    tdb: str

    def to_dict(self):
        return asdict(self)

    def append(self, data: "ParsedData"):
        if data.tdb != self.tdb:
            raise ValueError("TDB mismatch")
        self.elems.extend(data.elems)
        self.funcs.extend(data.funcs)
        self.phases.extend(data.phases)
        self.params.extend(data.params)


class TDBParser:
    def _parse_expression(self, line: str):
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

    def _parse_elem(self, line: str):
        match = re.match(r"ELEMENT\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*\!", line)
        if not match:
            return None
        elem, ref_state, atomic_mass, h298_h0, s298 = match.groups()
        return Elem(elem=elem, ref_state=ref_state, atomic_mass=float(atomic_mass),
                    h298_h0=float(h298_h0), s298=float(s298))

    def _parse_func(self, line: str):
        match = re.match(r"FUNCTION\s+(\S+)\s+(\S+)\s+([^\;].+)\s*\;\s+(\S+)\s+([YN]?)\s*\!", line)
        if not match:
            return None
        func, temp_start, expression, temp_end, is_continued = match.groups()
        elem = func.split("SER")[1]
        return Func(func=func, elem=elem, temp_start=float(temp_start), temp_end=float(temp_end),
                    expression=self._parse_expression(expression), is_continued=is_continued)

    def _parse_param(self, line: str, tdb: str):
        pattern = r"PARAMETER\s+([GL])\((\S+)\,(\S+)\;(\d)\)\s+(\S+)\s+([^\;].+)\s*\;\s+(\S+)\s+([YN]?)\s*\!"
        match = re.match(pattern, line)
        if not match:
            return None
        ptype, phase, components, order_num, temp_start, expression, \
            temp_end, is_continues = match.groups()
        return Param(ptype=ptype, phase=phase, components=components, order_num=int(order_num),
                     temp_start=float(temp_start), temp_end=float(temp_end),
                     expression=self._parse_expression(expression),
                     is_continued=is_continues, tdb=tdb, param="")

    def _parse_phase(self, line: str):
        # parts = line.split("PHASE")[1].strip().split()
        match = re.match(r"PHASE\s+(\S+)\s+\%\s+(\d+)\s+([^\!].+)\!", line)
        if not match:
            return None
        phase, sub_lattices, stoichiometry = match.groups()
        sub_lattices = int(sub_lattices)
        metrics = stoichiometry.split()
        sub_lattices = len(metrics) if len(metrics) != sub_lattices else sub_lattices
        return {"phase": phase, "sub_lattices": sub_lattices, "stoichiometry": " ".join(metrics)}

    def _parse_const(self, line: str):
        match = re.match(r"CONSTITUENT\s+(\S+)\s+([^\!].+)\!", line)
        if not match:
            return None
        phase, components = match.groups()
        components = components.strip(": ")
        return {"phase": phase, "components": components}

    def _export_elements(self, elements: list[Elem]):
        lines = ["$ ELEMENT NAME REF_STATE ATOMIC_MASS H298-H0 S298 !"]
        for e in elements:
            s1 = f"ELEMENT {e['elem']:2s} {e['ref_state']:21s}"
            s2 = f"{e['atomic_mass']:E} {e['h298_h0']:E} {e['s298']:E}"
            lines.append(f"{s1} {s2} !")

        lines.extend(["$ END ELEMENT !", ""])
        return lines

    def _export_functions(self, functions: list[Func]):
        lines = ["$ FUNCTION FUNC TEMP_START EXPRESSION TEMP_END IS_CONTINUED !"]
        for f in functions:
            ex1, ex2 = f["expression"].split("*T*LN(T)")
            s1 = f"FUNCTION {f['func']:7s} {f['temp_start']:.2f} {ex1}*T*LN(T)"
            s2 = f"    {ex2}; {f['temp_end']:.2f} {f['is_continued']:1s}"
            lines.extend([s1, f"{s2} !"])

        lines.extend(["$ END FUNCTION !", ""])
        return lines

    def _export_phase_and_params(self, phase: Phase, params: list[Param]):
        lines = ["$ PHASE SUB_LATTICES STOICHIOMETRY !", "$ CONSTITUENT PHASE COMPONENTS !"]
        s1 = f"PHASE {phase['phase']:7s} % {phase['sub_lattices']} {phase['stoichiometry']} !"
        lines.append(s1)
        comps = phase["components"].split(":")
        cur = f"CONSTITUENT {phase['phase']:7s} :"

        # Split to new lines
        for c in comps:
            if len(c) + len(cur) < 70:
                cur += f" {c}:"
            else:
                lines.append(cur)
                cur = f"    {c}:"
        lines.append(f"{cur} !")
        # s2 = f"CONSTITUENT {phase['phase']:7s} :{phase['components']}:"
        # lines.extend([s1, s2])
        for p in params:
            ex1, ex2 = p["expression"].split("*T", 1)
            if "*T**3" in ex2:
                ex2, ex3 = ex2.split("*T**3", 1)
                ex2 += "*T**3"
            else:
                ex3 = ""
            s1 = f"PARAMETER {p['param']} {p['temp_start']:.2f} {ex1}*T"
            s2 = f"    {ex2}"
            s3 = f"    {ex3}; {p['temp_end']:.2f} {p['is_continued']:1s}"
            lines.extend([s1, s2, f"{s3} !"])

        lines.extend(["$ PHASE AND PARAMETER DATA END !", ""])
        return lines
