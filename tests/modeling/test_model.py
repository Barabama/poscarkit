"""Unit tests for src/modeling/model.py — _integer_fractions, ModelStruct."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.modeling.model import _integer_fractions, _ask_normalize_fractions, ModelStruct
from src.modeling.base import Struct, Atom, SimplePoscar
from src.modeling.supercell import make_supercell


# ------------------------------------------------------------------ #
#  _integer_fractions                                                #
# ------------------------------------------------------------------ #

class TestIntegerFractions(unittest.TestCase):

    def test_single_element_exact(self):
        """Single element at 100% should produce exact integer count."""
        result, log = _integer_fractions(
            fracts={"Al": 1.0}, factors=(3, 3, 3), multipl=1
        )
        self.assertEqual(result, {"Al": 27})  # 1 * 1 * 27 = 27
        self.assertEqual(len(log), 1)
        self.assertEqual(log[0]["Element"], "Al")
        self.assertEqual(log[0]["Adjusted"], 27)
        self.assertEqual(log[0]["Delta"], 0)

    def test_multi_element_exact(self):
        """Multiple elements where integer conversion is exact."""
        result, log = _integer_fractions(
            fracts={"A": 0.5, "B": 0.5}, factors=(2, 2, 2), multipl=2
        )
        # 0.5 * 2 * 8 = 8 each, total 16 = 2*8
        self.assertEqual(result, {"A": 8, "B": 8})
        self.assertEqual(sum(result.values()), 16)
        for row in log:
            self.assertEqual(row["Delta"], 0)

    def test_rounding_adjustment(self):
        """Fractions that don't divide evenly should be adjusted."""
        result, log = _integer_fractions(
            fracts={"A": 0.3333, "B": 0.6667}, factors=(3, 3, 3), multipl=1
        )
        # 0.3333 * 27 = 8.9991 → rounds to 9
        # 0.6667 * 27 = 18.0009 → rounds to 18
        # Total = 27 = 1 * 27 ✓
        self.assertEqual(sum(result.values()), 27)
        for row in log:
            self.assertIn("Adjusted", row)
            self.assertIn("Delta", row)
            self.assertIn("Element", row)
            self.assertIn("SOF", row)

    def test_rounding_adjustment_needed(self):
        """When rounding total != target, adjustment corrects it."""
        result, log = _integer_fractions(
            fracts={"A": 1 / 3, "B": 1 / 3, "C": 1 / 3},
            factors=(2, 2, 2), multipl=3,
        )
        # Each: 1/3 * 3 * 8 = 8.0 → rounds to 8
        # Total = 24 = 3 * 8 ✓
        self.assertEqual(sum(result.values()), 24)
        self.assertEqual(len(result), 3)

    def test_empty_fracts_raises(self):
        """Empty fracts dict should raise ValueError."""
        with self.assertRaises(ValueError):
            _integer_fractions(fracts={}, factors=(2, 2, 2), multipl=1)

    def test_large_supercell(self):
        """Large supercell factors should still produce correct totals."""
        result, _ = _integer_fractions(
            fracts={"X": 0.4, "Y": 0.6}, factors=(10, 10, 10), multipl=2
        )
        # 0.4 * 2 * 1000 = 800; 0.6 * 2 * 1000 = 1200
        self.assertEqual(result["X"], 800)
        self.assertEqual(result["Y"], 1200)
        self.assertEqual(sum(result.values()), 2000)

    def test_log_rows_have_all_fields(self):
        """Verify log rows contain all expected fields."""
        _, log = _integer_fractions(
            fracts={"Fe": 0.7, "Cr": 0.3}, factors=(3, 3, 3), multipl=1
        )
        for row in log:
            for field in ["Element", "SOF", "Multiplicity", "Factor",
                          "Raw", "Rounded", "Adjusted", "Delta"]:
                self.assertIn(field, row, f"Missing field {field}")


# ------------------------------------------------------------------ #
#  _ask_normalize_fractions                                          #
# ------------------------------------------------------------------ #

class TestAskNormalizeFractions(unittest.TestCase):

    def test_already_normalized(self):
        """Fractions summing to 1.0 should be returned unchanged."""
        result = _ask_normalize_fractions("1a", {"A": 0.5, "B": 0.5})
        self.assertAlmostEqual(sum(result.values()), 1.0)
        self.assertAlmostEqual(result["A"], 0.5)

    def test_normalizes_when_sum_not_one(self):
        """Fractions not summing to 1.0 should be auto-normalized."""
        result = _ask_normalize_fractions("1a", {"A": 2.0, "B": 2.0})
        self.assertAlmostEqual(sum(result.values()), 1.0)
        self.assertAlmostEqual(result["A"], 0.5)
        self.assertAlmostEqual(result["B"], 0.5)

    def test_zero_sum_raises(self):
        """Zero-sum fractions should raise ValueError."""
        with self.assertRaises(ValueError):
            _ask_normalize_fractions("1a", {"A": 0.0, "B": 0.0})


# ------------------------------------------------------------------ #
#  ModelStruct integration                                           #
# ------------------------------------------------------------------ #

class TestModelStruct(unittest.TestCase):

    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())
        cell = np.array([
            [3.774, 0.0, 0.0],
            [0.0, 3.774, 0.0],
            [0.0, 0.0, 3.774],
        ])
        atom_list = [
            Atom(index=0, symbol="Au", coord=np.array([0.0, 0.0, 0.0]),
                 constr=["T", "T", "T"], note="1a-Au"),
            Atom(index=1, symbol="Cu", coord=np.array([0.0, 0.5, 0.5]),
                 constr=["T", "T", "T"], note="3c-Cu"),
            Atom(index=2, symbol="Cu", coord=np.array([0.5, 0.0, 0.5]),
                 constr=["T", "T", "T"], note="3c-Cu"),
            Atom(index=3, symbol="Cu", coord=np.array([0.5, 0.5, 0.0]),
                 constr=["T", "T", "T"], note="3c-Cu"),
        ]
        self.unitcell = Struct(cell=cell, is_direct=True, atom_list=atom_list)
        self.poscar_path = self.temp_dir / "POSCAR.vasp"
        SimplePoscar.write_poscar(self.poscar_path, self.unitcell, "test")

        self.structure_info = {
            "cell": [[3.774, 0.0, 0.0], [0.0, 3.774, 0.0], [0.0, 0.0, 3.774]],
            "1a": {
                "atoms": ("Au", [[0.0, 0.0, 0.0]]),
                "sofs": {"V": 1.0},
            },
            "3c": {
                "atoms": ("Cu", [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
                "sofs": {"Co": 0.4444, "Ni": 0.4444, "V": 0.1112},
            },
        }

    def tearDown(self):
        import shutil
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_gen_site_integers(self):
        """_gen_site_integers should produce correct integer mappings."""
        modeler = ModelStruct(
            name="test", poscar=self.poscar_path,
            factors=(3, 3, 3), struct_info=self.structure_info,
            batch_size=1,
        )
        si = modeler.site_integers
        self.assertIn(("1a", "Au"), si)
        self.assertIn(("3c", "Cu"), si)
        # 1a-Au: 1 atom * 27 = 27, 100% V → V:27
        self.assertEqual(si[("1a", "Au")], {"V": 27})
        # 3c-Cu: 3 atoms * 27 = 81
        self.assertEqual(sum(si[("3c", "Cu")].values()), 81)

    def test_gen_site_integers_empty_sofs(self):
        """Empty SOFs should keep the default element at 100%."""
        info = {
            "cell": [[3.774, 0.0, 0.0], [0.0, 3.774, 0.0], [0.0, 0.0, 3.774]],
            "1a": {
                "atoms": ("Au", [[0.0, 0.0, 0.0]]),
                "sofs": {},
            },
        }
        modeler = ModelStruct(
            name="test", poscar=self.poscar_path,
            factors=(2, 2, 2), struct_info=info,
            batch_size=1,
        )
        si = modeler.site_integers
        self.assertIn(("1a", "Au"), si)
        self.assertEqual(si[("1a", "Au")], {"Au": 8})  # 1 * 8 = 8

    def test_gen_site_integers_no_struct_info(self):
        """Missing structure_info should produce empty site_integers."""
        modeler = ModelStruct(
            name="test", poscar=self.poscar_path,
            factors=(2, 2, 2), struct_info={},
            batch_size=1,
        )
        self.assertEqual(dict(modeler.site_integers), {})

    def test_allocate_atoms_shuffle(self):
        """_allocate_atoms should assign symbols and meta correctly."""
        modeler = ModelStruct(
            name="test", poscar=self.poscar_path,
            factors=(2, 2, 2), struct_info=self.structure_info,
            batch_size=1,
        )
        supercell = make_supercell(modeler.unitcell, (2, 2, 2))
        site_substruct = modeler._allocate_atoms(
            supercell=supercell,
            site_integers=modeler.site_integers,
            seed=42,
        )
        self.assertIn(("1a", "Au"), site_substruct)
        self.assertIn(("3c", "Cu"), site_substruct)

        # Verify symbol assignment
        struct_1a = site_substruct[("1a", "Au")]
        # 1a-Au: 1 atom * 2*2*2 = 8 atoms, all should be "V" (SOF: V=1.0)
        self.assertEqual(len(struct_1a), 8)
        for atom in struct_1a:
            self.assertEqual(atom.symbol, "V")

        # Verify meta fields are set
        for atom in struct_1a:
            self.assertIsNotNone(atom.meta)
            self.assertTrue(atom.meta.startswith("s"))

    def test_allocate_atoms_reproducible(self):
        """Same seed produces same allocation."""
        modeler = ModelStruct(
            name="test", poscar=self.poscar_path,
            factors=(3, 3, 3), struct_info=self.structure_info,
            batch_size=1,
        )
        supercell = make_supercell(modeler.unitcell, (3, 3, 3))

        result1 = modeler._allocate_atoms(
            supercell=supercell.copy(),
            site_integers=modeler.site_integers,
            seed=42,
        )
        result2 = modeler._allocate_atoms(
            supercell=supercell.copy(),
            site_integers=modeler.site_integers,
            seed=42,
        )
        # Same seed should produce identical symbol sequences
        symbols1 = [a.symbol for a in result1[("3c", "Cu")]]
        symbols2 = [a.symbol for a in result2[("3c", "Cu")]]
        self.assertEqual(symbols1, symbols2)

    def test_sof_integer_log(self):
        """ModelStruct should produce SOF integer conversion log."""
        modeler = ModelStruct(
            name="test", poscar=self.poscar_path,
            factors=(3, 3, 3), struct_info=self.structure_info,
            batch_size=1,
        )
        log = modeler._sof_integer_log
        self.assertGreater(len(log), 0)
        for row in log:
            self.assertIn("Site", row)
            self.assertIn("Element", row)

    def test_gen_site_integers_missing_atoms_info(self):
        """Missing atoms_info should raise ValueError."""
        bad_info = {
            "cell": [[3.774, 0.0, 0.0], [0.0, 3.774, 0.0], [0.0, 0.0, 3.774]],
            "1a": {"sofs": {"V": 1.0}},
            # missing "atoms" key
        }
        with self.assertRaises(ValueError):
            ModelStruct(
                name="test", poscar=self.poscar_path,
                factors=(2, 2, 2), struct_info=bad_info,
                batch_size=1,
            )

    def test_gen_site_integers_missing_sofs_info(self):
        """Missing sofs_info should raise ValueError."""
        bad_info = {
            "cell": [[3.774, 0.0, 0.0], [0.0, 3.774, 0.0], [0.0, 0.0, 3.774]],
            "1a": {"atoms": ("Au", [[0.0, 0.0, 0.0]])},
            # missing "sofs" key
        }
        with self.assertRaises(ValueError):
            ModelStruct(
                name="test", poscar=self.poscar_path,
                factors=(2, 2, 2), struct_info=bad_info,
                batch_size=1,
            )


if __name__ == "__main__":
    unittest.main()
