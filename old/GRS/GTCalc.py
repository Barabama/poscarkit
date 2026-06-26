# GRSFitter.py

from itertools import groupby, product
import os
import re
import sys

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sympy import lambdify, symbols, simplify

from TDBManager import Elem, Func, Param, Phase, ParsedData, TDBParser


def average_site_fracts(site_fracts: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
    averaged = {}
    for site, elem_fracts in site_fracts.items():
        elems = list(elem_fracts.keys())
        avg_value = 1.0 / len(elems)
        averaged[site] = {elem: avg_value for elem in elems}
    return averaged


class GTCalculator(TDBParser):
    parsed: ParsedData

    def __init__(self, tdb_path: str):
        self.parse_tdb(tdb_path=tdb_path)

    def parse_tdb(self, tdb_path: str, tdb: str = ""):
        """Parse TheromDynamics TDB file.

        Args:
            tdb_path (str): Tdb file path.
            tdb (str, optional): Tdb name. Defaults to "".

        Returns:
            Dict[str, Elem | Func | Phase | Param]: Parsed data.
        """
        tdb = tdb or os.path.splitext(os.path.basename(tdb_path))[0]
        with open(tdb_path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        # Clean notes
        text = "".join(s.split("$", 1)[0].strip() for s in lines)
        elems: list[Elem] = []
        funcs: list[Func] = []
        phases: list[Phase] = []
        params: list[Param] = []
        merged = []
        for line in [s.strip() for s in text.split("!")]:
            line += " !"
            if line.startswith("ELEMENT"):
                if e := self._parse_elem(line):
                    elems.append(e)
            elif line.startswith("FUNCTION"):
                if f := self._parse_func(line):
                    funcs.append(f)
            elif line.startswith("PARAMETER"):
                if p := self._parse_param(line, tdb):
                    params.append(p)
            elif line.startswith("PHASE"):
                if p := self._parse_phase(line):
                    merged.append(p)
            elif line.startswith("CONSTITUENT"):
                if c := self._parse_const(line):
                    merged.append(c)
            else:
                pass

        # Merge phases and constituents
        merged.sort(key=lambda x: x["phase"])
        for p, group in groupby(merged, lambda x: x["phase"]):
            data = {}
            for g in group:
                data.update({"tdb": tdb, **g})
            phases.append(Phase(**data))

        self.parsed = ParsedData(elems=elems, funcs=funcs, phases=phases, params=params, tdb=tdb)

    def _get_func_expr(self, func: str):
        funcs = {f["func"]: f for f in self.parsed.funcs}
        return funcs[func]["expression"]

    def _get_param_expr(self, phase: str, comps: list[str]):
        name = f"{phase},{':'.join(comps)}"
        params = {f"{p['phase']},{p['components']}": p for p in self.parsed.params}
        return params[name]["expression"]

    def merge_expr(self, site_fracts: dict[str, dict[str, float]], phase: str):
        # site_fracts = {'1a': {'A': 0.5, 'B': 0.5}, '1b': {'A': 0.5, 'B': 0.5}}
        if self.parsed is None:
            raise ValueError("No parsed data")
        metrics = [int(site[0]) for site, data in site_fracts.items()]
        m_sum = sum(m for m in metrics)
        metrics = [m / m_sum for m in metrics]
        groups = [[(elem, sof) for elem, sof in data.items()]
                  for site, data in site_fracts.items()]

        # Gibbs free energy of End-Members
        G_ems: list[str] = []
        for group in list(product(*groups)):
            comps = [elem.upper() for elem, sof in group]
            expr_param = self._get_param_expr(phase, comps)
            if matches := re.findall(r"(SER[A-Za-z]+\#)", expr_param):
                for match in matches:
                    expr_func = self._get_func_expr(match.split("#")[0])
                    expr_param = expr_param.replace(match, f"({expr_func})")

            sofs = [str(s) for e, s in group]
            G_ems.append(f"{'*'.join(sofs)}*({expr_param})")
        # print(G_ems)

        # Configurational entropy
        S_confs: list[str] = []
        for metry, group in zip(metrics, groups):
            s = "+".join(f"{s}*LN({s})" for e, s in group)
            S_confs.append(f"{metry}*({s})")
        # print(S_confs)

        expr_all = f"{'+'.join(G_ems)}+8.314*T*({'+'.join(S_confs)})"
        expr_all = simplify(expr_all.replace("LN", "ln"))
        # print(expr_all)

        return expr_all

    def check_sofs(self, site_fracts: dict[str, dict[str, float]]):
        for site, elem_fracts in site_fracts.items():
            if abs(sum(elem_fracts.values())) - 1 > 1e-6:
                raise ValueError(f"Site fractions for {site} do not sum to 1")
        return site_fracts

    def calc_formula(self, expr: str, t_values=range(10, 2010, 10)):
        T = symbols("T")
        f = lambdify(T, expr, "numpy")
        g_values = [f(t) for t in t_values]
        df = pd.DataFrame({"T": t_values, "G": g_values})
        return df

    def handle(self, site_fracts: dict[str, dict[str, float]], phase: str):
        sofs_des = self.check_sofs(site_fracts)
        sofs_avg = average_site_fracts(site_fracts)

        expr_avg = self.merge_expr(sofs_avg, phase)
        expr_des = self.merge_expr(sofs_des, phase)

        df_avg = self.calc_formula(expr_avg)
        df_des = self.calc_formula(expr_des)

        elems = sorted(list(set(e.upper() for d in site_fracts.values() for e in d)))
        elem_str = "".join(elems)
        output = os.path.join(os.path.dirname(sys.argv[0]), "output")
        if not os.path.exists(output):
            os.mkdir(output)
        
        # Save to excel
        with pd.ExcelWriter(os.path.join(output, f"G(T)-{elem_str}-{phase}.xlsx")) as writer:
            df_avg.to_excel(writer, sheet_name="Average", index=False)
            df_des.to_excel(writer, sheet_name="Designated", index=False)

        # Plot
        plt.figure(figsize=(8, 6))
        plt.plot(df_avg["T"], df_avg["G"], label="Average", color="red", marker="x")
        plt.plot(df_des["T"], df_des["G"], label="Designated", color="blue", marker="o")
        plt.title(f"G-T of {elem_str}_{phase}")
        plt.xlabel("T")
        plt.ylabel("G")
        plt.legend()
        plt.grid(True)

        plt.savefig(os.path.join(output, f"G(T)-{elem_str}-{phase}.png"))
        plt.show()
        plt.close()
