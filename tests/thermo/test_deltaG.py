"""Tests for src/thermo/deltaG.py — case-insensitive TDB end-member lookup."""

import os
import tempfile
import unittest

import numpy as np

from src.thermo.tdb import parse_tdb
from src.thermo.deltaG import calc_deltaG_random


_TDB_SAMPLE = """
$ Sample TDB
FUNCTION SERCO 1.00 -6.815E+5+7.178E+1*T-14.54*T*LN(T); 6000.00 N !
FUNCTION SERCR 1.00 -6.0E+5+6.5E+1*T; 6000.00 N !
PARAMETER G(FCC,CO:CR;0) 1.00 +5000-2*T-.25*SERCO#-.75*SERCR#; 6000.00 N REF0 !
PARAMETER G(FCC,CR:CO;0) 1.00 +6000-2*T-.25*SERCR#-.75*SERCO#; 6000.00 N REF0 !
PARAMETER G(FCC,CO:CO;0) 1.00 +10000-5*T; 6000.00 N REF0 !
PARAMETER G(FCC,CR:CR;0) 1.00 -10000+3*T; 6000.00 N REF0 !
PHASE FCC % 2 0.25 0.75 !
"""


class TestDeltaGCaseInsensitive(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.NamedTemporaryFile(
            mode="w", suffix=".TDB", delete=False, encoding="utf-8"
        )
        self.tmp.write(_TDB_SAMPLE)
        self.tmp.close()
        self.addCleanup(os.unlink, self.tmp.name)
        self.tdb = parse_tdb(self.tmp.name)

    def test_canonical_composition_matches_uppercase_tdb_keys(self):
        """Canonical element symbols ('Co') must match all-caps TDB keys.

        Regression: SOFData canonicalizes reader elements to 'Co'/'Cr',
        while TDB parameter signatures keep the file's all-caps 'CO'/'CR'.
        calc_deltaG_random must be case-insensitive on both sides.
        """
        T = np.array([1000.0, 1100.0])
        comp = {"Co": 0.4, "Cr": 0.6}
        result = calc_deltaG_random(
            comp,
            self.tdb["functions"],
            self.tdb["params"],
            T,
            sconf_random=5.0,
            phase="FCC",
            site_ratios=[0.25, 0.75],
        )
        self.assertEqual(result.shape, T.shape)
        self.assertTrue(np.all(np.isfinite(result)))

    def test_uppercase_composition_still_works(self):
        """All-caps composition (legacy) also resolves."""
        T = np.array([1000.0])
        comp = {"CO": 0.4, "CR": 0.6}
        result = calc_deltaG_random(
            comp,
            self.tdb["functions"],
            self.tdb["params"],
            T,
            sconf_random=5.0,
            phase="FCC",
            site_ratios=[0.25, 0.75],
        )
        self.assertTrue(np.all(np.isfinite(result)))

    def test_parse_parameter_uppercases_components(self):
        """TDB parser stores element components all-caps."""
        self.assertIn("FCC,CO:CR", self.tdb["params"])
        self.assertIn("FCC,CR:CO", self.tdb["params"])


if __name__ == "__main__":
    unittest.main()
