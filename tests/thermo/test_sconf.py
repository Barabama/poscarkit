"""Tests for src/thermo/sconf.py — generalized N-sublattice entropy."""

import unittest

import numpy as np

from src.thermo.sconf import calc_sconf, calc_sconf_random
from src.io.ir import PHASE_SITE_RATIOS

R = 8.314


class TestSconf(unittest.TestCase):

    def test_sconf_random_quaternary(self):
        """Co0.3Cr0.2Fe0.1Ni0.4: -R * sum(x_i*ln(x_i))"""
        comp = {"CO": 0.3, "CR": 0.2, "FE": 0.1, "NI": 0.4}
        expected = -R * (0.3 * np.log(0.3) + 0.2 * np.log(0.2) +
                         0.1 * np.log(0.1) + 0.4 * np.log(0.4))
        result = calc_sconf_random(comp, PHASE_SITE_RATIOS["FCC"])
        self.assertAlmostEqual(result, expected, places=10)

    def test_sconf_random_equimolar_ternary(self):
        comp = {"FE": 1 / 3, "MN": 1 / 3, "NI": 1 / 3}
        expected = -R * np.log(1 / 3)
        result = calc_sconf_random(comp, PHASE_SITE_RATIOS["FCC"])
        self.assertAlmostEqual(result, expected, places=10)

    def test_sconf_random_pure(self):
        self.assertAlmostEqual(
            calc_sconf_random({"NI": 1.0}, PHASE_SITE_RATIOS["FCC"]), 0.0, places=10
        )

    def test_sconf_random_bcc_equals_fcc(self):
        """BCC same as FCC for random since ratios sum to 1."""
        comp = {"A": 0.5, "B": 0.5}
        s_fcc = calc_sconf_random(comp, PHASE_SITE_RATIOS["FCC"])
        s_bcc = calc_sconf_random(comp, PHASE_SITE_RATIOS["BCC"])
        self.assertAlmostEqual(s_fcc, s_bcc, places=10)

    def test_calc_sconf_pure(self):
        sofs = [{"NI": 1.0, "VA": 0.0}, {"NI": 1.0, "VA": 0.0}]
        self.assertAlmostEqual(
            calc_sconf(sofs, PHASE_SITE_RATIOS["FCC"]), 0.0, places=10
        )

    def test_calc_sconf_identical_sublattices(self):
        comp = {"A": 0.5, "B": 0.5}
        s = calc_sconf([comp, comp], PHASE_SITE_RATIOS["FCC"])
        expected = -R * (0.5 * np.log(0.5) + 0.5 * np.log(0.5))
        self.assertAlmostEqual(s, expected, places=10)

    def test_calc_sconf_near_zero(self):
        sofs_1 = {"A": 1.0 - 1e-40, "B": 1e-40}
        sofs_2 = {"A": 0.5, "B": 0.5}
        s = calc_sconf([sofs_1, sofs_2], PHASE_SITE_RATIOS["FCC"])
        expected = -R * 0.75 * (0.5 * np.log(0.5) + 0.5 * np.log(0.5))
        self.assertAlmostEqual(s, expected, places=4)


if __name__ == "__main__":
    unittest.main()
