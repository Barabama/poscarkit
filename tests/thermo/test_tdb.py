"""Tests for src/thermo/tdb.py — TDB file parsing."""

import os
import tempfile
import unittest

from src.thermo.tdb import parse_tdb, resolve_expression


_TDB_SAMPLE = """
$ Sample TDB for testing
FUNCTION SERCO   1.00 -6.815E+5+7.178E+1*T-...; 6000.00 N !
FUNCTION SERCR   1.00 -6.0E+5+6.5E+1*T-...; 6000.00 N !
FUNCTION GHSERFE 1.00 -7.0E+5+7.0E+1*T-...; 6000.00 N !

PARAMETER G(FCC,CO:CR;0) 1.00 +5000-2*T; 6000.00 N !
PARAMETER G(FCC,CO:CO;0) 1.00 +10000-5*T; 6000.00 N !
PARAMETER G(FCC,CR:CR;0) 1.00 -10000+3*T; 6000.00 N !

PHASE FCC % 2 0.25 0.75 !
PHASE BCC % 2 0.5 0.5 !
"""


class TestTDB(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".TDB", delete=False, encoding="utf-8"
        )
        self.tmp.write(_TDB_SAMPLE)
        self.tmp.close()

    def tearDown(self):
        os.unlink(self.tmp.name)

    def test_parse_functions(self):
        tdb = parse_tdb(self.tmp.name)
        self.assertIn("SERCO", tdb["functions"])
        self.assertIn("SERCR", tdb["functions"])
        self.assertIn("GHSERFE", tdb["functions"])
        self.assertIn("-6.815E+5", tdb["functions"]["SERCO"])

    def test_parse_params_fcc(self):
        tdb = parse_tdb(self.tmp.name, phase="FCC")
        self.assertIn("FCC,CO:CR", tdb["params"])
        self.assertIn("FCC,CO:CO", tdb["params"])
        self.assertEqual(len(tdb["params"]), 3)

    def test_parse_params_filtered(self):
        tdb = parse_tdb(self.tmp.name, phase="BCC")
        self.assertEqual(len(tdb["params"]), 0)

    def test_parse_params_all_phases(self):
        tdb = parse_tdb(self.tmp.name)
        self.assertGreaterEqual(len(tdb["params"]), 3)

    def test_parse_phase_ratios(self):
        tdb = parse_tdb(self.tmp.name)
        self.assertIn("FCC", tdb["phase_ratios"])
        self.assertAlmostEqual(sum(tdb["phase_ratios"]["FCC"]), 1.0, places=6)
        self.assertEqual(tdb["phase_ratios"]["FCC"], [0.25, 0.75])
        self.assertEqual(tdb["phase_ratios"]["BCC"], [0.5, 0.5])

    def test_resolve_expression(self):
        tdb = parse_tdb(self.tmp.name)
        expr = resolve_expression("SERCO# + SERCR#", tdb["functions"])
        self.assertIn("-6.815E+5", expr)
        self.assertIn("-6.0E+5", expr)
        self.assertNotIn("SERCO#", expr)
        self.assertNotIn("SERCR#", expr)

    def test_resolve_expression_ghser_fallback(self):
        tdb = parse_tdb(self.tmp.name)
        expr = resolve_expression("GHSERFE#", tdb["functions"])
        self.assertIn("-7.0E+5", expr)
        self.assertNotIn("GHSERFE#", expr)


if __name__ == "__main__":
    unittest.main()
