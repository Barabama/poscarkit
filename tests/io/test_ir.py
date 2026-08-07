import unittest

import numpy as np

from src.io.ir import SOFData


def _make_sofdata(**overrides):
    """Build a SOFData with all-caps element symbols (as exported by
    Thermo-Calc / Pandat files)."""
    defaults = dict(
        source_path="dummy.csv",
        phase="FCC",
        site_ratios=[0.25, 0.75],
        T=np.array([1000.0, 1100.0]),
        Y_subl=[
            {"CO": np.array([1.0, 1.0])},
            {"CO": np.array([0.3, 0.3]), "CU": np.array([0.4, 0.4]),
             "NI": np.array([0.3, 0.3])},
        ],
        composition={"CO": 0.525, "CU": 0.3, "NI": 0.175},
        elements=["CO", "CU", "NI"],
    )
    defaults.update(overrides)
    return SOFData(**defaults)


class TestSOFDataCanonicalization(unittest.TestCase):
    """SOFData canonicalizes all-caps element symbols on construction."""

    def test_elements_canonicalized(self):
        ir = _make_sofdata()
        self.assertEqual(ir.elements, ["Co", "Cu", "Ni"])

    def test_composition_keys_canonicalized(self):
        ir = _make_sofdata()
        self.assertEqual(sorted(ir.composition.keys()), ["Co", "Cu", "Ni"])

    def test_ysubl_keys_canonicalized(self):
        ir = _make_sofdata()
        for y in ir.Y_subl:
            for elem in y:
                self.assertTrue(
                    elem.istitle() and len(elem) <= 2,
                    f"non-canonical symbol {elem!r} in Y_subl",
                )
        self.assertIn("Co", ir.Y_subl[0])
        self.assertIn("Ni", ir.Y_subl[1])

    def test_validations_still_run(self):
        """Canonicalization does not bypass structural validation."""
        with self.assertRaises(ValueError):
            # length mismatch: 3 site ratios vs 2 Y_subl entries
            _make_sofdata(site_ratios=[0.4, 0.3, 0.3])
        with self.assertRaises(ValueError):
            # does not sum to 1.0: 0.2 + 0.7 = 0.9
            _make_sofdata(site_ratios=[0.2, 0.7])


if __name__ == "__main__":
    unittest.main()
