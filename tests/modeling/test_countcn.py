import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.modeling.countcn import CNCounter
from src.modeling.base import Struct, Atom, SimplePoscar


class TestCNCounter(unittest.TestCase):
    """Test cases for CNCounter class."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a simple cubic structure for testing
        cell = np.array([[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]])

        atom_list = [
            Atom(index=0, symbol="Al", coord=np.array([0.0, 0.0, 0.0])),
            Atom(index=1, symbol="Al", coord=np.array([0.5, 0.5, 0.0])),
            Atom(index=2, symbol="Al", coord=np.array([0.5, 0.0, 0.5])),
            Atom(index=3, symbol="Al", coord=np.array([0.0, 0.5, 0.5])),
        ]

        self.test_struct = Struct(cell=cell, is_direct=True, atom_list=atom_list)

        # Create temporary directory and POSCAR file
        self.temp_dir = Path(tempfile.mkdtemp())
        self.poscar_path = self.temp_dir / "POSCAR"
        SimplePoscar.write_poscar(self.poscar_path, self.test_struct, "Test structure")

        self.counter = CNCounter(name="test", poscar=self.poscar_path)

    def test_countCN2files(self):
        """Test countCN2files method."""
        # Run the countCN2files method with by_ase=False (default)
        output_dir = self.counter.countCN2files(
            outdir=self.temp_dir, cutoff_mult=1.1, parallel=2, by_ase=False
        )

        # Check that output directory was created
        self.assertTrue(output_dir.exists())
        self.assertTrue(output_dir.is_dir())

        # Check that some output files were created
        csv_file = output_dir / "test-d1nn-cn-count.csv"
        self.assertTrue(csv_file.exists())

        png_files = list(output_dir.glob("*.png"))
        self.assertGreater(len(png_files), 0)

        # Check cumulative POSCAR files exist (all Al atoms have CN=3)
        cum_files = list(output_dir.glob("*-d1nn-*.vasp"))
        # At minimum, cumulative file test-d1nn-Al-Al-3+.vasp should exist
        cum_file = output_dir / "test-d1nn-Al-Al-3+.vasp"
        self.assertTrue(cum_file.exists(),
                       f"Cumulative file {cum_file} not found among {cum_files}")

        # Check CSV has Cumulative column
        import pandas as pd
        df = pd.read_csv(csv_file)
        self.assertIn("Cumulative", df.columns)

    def test_cumulative_cn_bi_element(self):
        """Test cumulative CN with a bi-element structure with varying CN.

        Structure: Fe-O in a 20x20x20 cubic cell.
        Fe1 at (0,0,0) has 6 O neighbours (CN=6).
        Fe2 at (10,10,10) has 4 O neighbours (CN=4).
        O-O distances ~4.24 > cutoff=3.3, so no O-O pairs.
        Cluster separation ~17.3 >> cutoff.
        """
        cell = np.eye(3) * 20.0
        atoms = [
            # Cluster 1: Fe1 + 6 O
            Atom(index=0, symbol="Fe", coord=np.array([0.0, 0.0, 0.0])),
            Atom(index=1, symbol="O", coord=np.array([3.0, 0.0, 0.0])),
            Atom(index=2, symbol="O", coord=np.array([0.0, 3.0, 0.0])),
            Atom(index=3, symbol="O", coord=np.array([0.0, 0.0, 3.0])),
            Atom(index=4, symbol="O", coord=np.array([-3.0, 0.0, 0.0])),
            Atom(index=5, symbol="O", coord=np.array([0.0, -3.0, 0.0])),
            Atom(index=6, symbol="O", coord=np.array([0.0, 0.0, -3.0])),
            # Cluster 2: Fe2 + 4 O
            Atom(index=7, symbol="Fe", coord=np.array([10.0, 10.0, 10.0])),
            Atom(index=8, symbol="O", coord=np.array([13.0, 10.0, 10.0])),
            Atom(index=9, symbol="O", coord=np.array([10.0, 13.0, 10.0])),
            Atom(index=10, symbol="O", coord=np.array([10.0, 10.0, 13.0])),
            Atom(index=11, symbol="O", coord=np.array([7.0, 10.0, 10.0])),
        ]
        struct = Struct(cell=cell, is_direct=False, atom_list=atoms)

        tmp = Path(tempfile.mkdtemp())
        poscar = tmp / "POSCAR"
        SimplePoscar.write_poscar(poscar, struct, "Fe-O test")
        counter = CNCounter(name="fe-o-cum", poscar=poscar)

        try:
            out = counter.countCN2files(outdir=tmp, cutoff_mult=1.1)
            self.assertTrue(out.is_dir())

            # 1. Cumulative files exist with correct names
            fe_o_6p = out / "fe-o-cum-d1nn-Fe-O-6+.vasp"
            fe_o_4p = out / "fe-o-cum-d1nn-Fe-O-4+.vasp"
            self.assertTrue(fe_o_6p.exists(), f"Expected {fe_o_6p.name}")
            self.assertTrue(fe_o_4p.exists(), f"Expected {fe_o_4p.name}")
            # No CN=5 level should exist
            self.assertFalse((out / "fe-o-cum-d1nn-Fe-O-5+.vasp").exists())

            # 2. Atom counts in cumulative files
            struct_6p = SimplePoscar.read_poscar(fe_o_6p)
            self.assertEqual(len(struct_6p.atom_list), 7,
                           "6+ should have Fe1 + 6 O = 7 atoms")
            struct_4p = SimplePoscar.read_poscar(fe_o_4p)
            self.assertEqual(len(struct_4p.atom_list), 12,
                           "4+ should have Fe1 + Fe2 + 10 O = 12 atoms")

            # 3. Cumulative >= exact at same CN level
            exact_6 = out / "fe-o-cum-d1nn-Fe-O-6.vasp"
            exact_4 = out / "fe-o-cum-d1nn-Fe-O-4.vasp"
            struct_exact6 = SimplePoscar.read_poscar(exact_6)
            self.assertGreaterEqual(len(struct_6p.atom_list),
                                   len(struct_exact6.atom_list))
            struct_exact4 = SimplePoscar.read_poscar(exact_4)
            self.assertGreaterEqual(len(struct_4p.atom_list),
                                   len(struct_exact4.atom_list))

            # 4. CSV columns
            import pandas as pd
            csv_file = out / "fe-o-cum-d1nn-cn-count.csv"
            self.assertTrue(csv_file.exists())
            df = pd.read_csv(csv_file)
            self.assertEqual(list(df.columns), ["CN", "Count", "Cumulative"])

            # 5. Fe-O rows with correct values
            fe_o_rows = df[df["CN"].str.contains("Fe.*O")]
            self.assertEqual(len(fe_o_rows), 2)

            row6 = fe_o_rows[fe_o_rows["CN"] == "Fe*-6O"]
            row4 = fe_o_rows[fe_o_rows["CN"] == "Fe*-4O"]
            self.assertFalse(row6.empty)
            self.assertFalse(row4.empty)
            row6 = row6.iloc[0]
            row4 = row4.iloc[0]

            self.assertEqual(row6["Count"], 1)
            self.assertEqual(row6["Cumulative"], 1)   # only Fe1 has CN >= 6
            self.assertEqual(row4["Count"], 1)
            self.assertEqual(row4["Cumulative"], 2)   # both Fe atoms have CN >= 4

            # 6. Cumulative monotonically non-decreasing as CN decreases
            self.assertLessEqual(row6["Cumulative"], row4["Cumulative"])

            # 7. At max CN, Cumulative == Count
            self.assertEqual(row6["Cumulative"], row6["Count"])

            # 8. At min CN for Fe-O, Cumulative == total Fe atoms
            self.assertEqual(row4["Cumulative"], 2)
        finally:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)

    def test_cumulative_empty(self):
        """Cumulative structs: empty cndata_list returns empty dict."""
        # Single atom — no pairs possible, cndata_list stays empty.
        # We set cndata_list directly to avoid detect_cutoff ValueError
        # on a single-atom structure.
        cell = np.eye(3) * 10.0
        atoms = [Atom(index=0, symbol="He", coord=np.array([0.0, 0.0, 0.0]))]
        struct = Struct(cell=cell, is_direct=False, atom_list=atoms)

        tmp = Path(tempfile.mkdtemp())
        poscar = tmp / "POSCAR"
        SimplePoscar.write_poscar(poscar, struct, "Empty test")
        counter = CNCounter(name="empty-cum", poscar=poscar)

        try:
            counter.struct = struct
            counter.cndata_list = []
            cum_structs = counter.generate_cn_cumulative_structs()
            self.assertEqual(len(cum_structs), 0,
                           "Empty cndata_list should produce no cumulative structs")
        finally:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)

    def test_countCN2files_ase(self):
        """Test countCN2files method with ASE."""
        # Run the countCN2files method with by_ase=True
        output_dir = self.counter.countCN2files(
            outdir=self.temp_dir, cutoff_mult=1.1, parallel=2, by_ase=True
        )

        # Check that output directory was created
        self.assertTrue(output_dir.exists())
        self.assertTrue(output_dir.is_dir())

        # Check that some output files were created
        csv_file = output_dir / "test-d1nn-cn-count.csv"
        self.assertTrue(csv_file.exists())

    def tearDown(self):
        """Clean up test fixtures."""
        # Clean up temporary directory
        import shutil

        shutil.rmtree(self.temp_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
