import sys
import io
import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch

from src.cli.poscarkit import main, cmd_help, cmd_supercell, cmd_modeling, cmd_compare, cmd_merge, cmd_separate
import argparse

poscar_content = """SYSTEM=Au1Cu3-FCC
 1.0000000000000000
    3.7740000000000000    0.0000000000000000    0.0000000000000000
    0.0000000000000000    3.7740000000000000    0.0000000000000000
    0.0000000000000000    0.0000000000000000    3.7740000000000000
  Au  Cu
   1   3
Selective dynamics
Direct
 0.0000000000000000 0.0000000000000000 0.0000000000000000 T T T # 1a-Au-#1
 0.0000000000000000 0.5000000000000000 0.5000000000000000 T T T # 3c-Cu-#1
 0.5000000000000000 0.0000000000000000 0.5000000000000000 T T T # 3c-Cu-#2
 0.5000000000000000 0.5000000000000000 0.0000000000000000 T T T # 3c-Cu-#3
"""

poscar_content2 = """SYSTEM=Au1Cu3-FCC
 1.0000000000000000
    3.7740000000000000    0.0000000000000000    0.0000000000000000
    0.0000000000000000    3.7740000000000000    0.0000000000000000
    0.0000000000000000    0.0000000000000000    3.7740000000000000
  Au  Cu
   1   3
Selective dynamics
Direct
 0.0000000000000000 0.0000000000000000 0.0000000000000000 T T T # 1a-Au-#1
 0.0000000000000000 0.5000000000000000 0.5000000000000000 T T T # 3c-Cu-#1
 0.5000000000000000 0.0000000000000000 0.5000000000000000 T T T # 3c-Cu-#2
 0.5000000000000000 0.5000000000000000 0.2500000000000000 T T T # 3c-Cu-#3
"""

poscar_content3 = """SYSTEM=Au1Cu3-FCC
 1.0000000000000000
    3.7740000000000000    0.0000000000000000    0.0000000000000000
    0.0000000000000000    3.7740000000000000    0.0000000000000000
    0.0000000000000000    0.0000000000000000    3.7740000000000000
  Au  Cu
   2   2
Selective dynamics
Direct
 0.0000000000000000 0.0000000000000000 0.0000000000000000 T T T # 1a-Au-#1
 0.0000000000000000 0.5000000000000000 0.5000000000000000 T T T # 3c-Cu-#1
 0.5000000000000000 0.0000000000000000 0.5000000000000000 T T T # 3c-Cu-#2
 0.5000000000000000 0.5000000000000000 0.0000000000000000 T T T # 1a-Au-#2
"""


class TestCLI(unittest.TestCase):
    """Test CLI commands."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(tempfile.mkdtemp())
        self.poscar_file = self.test_dir / "POSCAR.vasp"
        with open(self.poscar_file, "w") as f:
            f.write(poscar_content)
        self.poscar_file2 = self.test_dir / "POSCAR2.vasp"
        with open(self.poscar_file2, "w") as f:
            f.write(poscar_content2)
        self.poscar_file3 = self.test_dir / "POSCAR3.vasp"
        with open(self.poscar_file3, "w") as f:
            f.write(poscar_content3)

    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_cmd_help(self):
        """Test help command."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            cmd_help(argparse.Namespace())
            output = fake_out.getvalue()
            self.assertIn("POSCARKIT", output)
            self.assertIn("help", output)
            self.assertIn("modeling", output)
            self.assertIn("supercell", output)
            self.assertIn("EXAMPLES", output)

    def test_cmd_supercell(self):
        """Test supercell command."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            poscar=str(self.poscar_file),
            outdir=str(outdir),
            factors=[2, 2, 2],
            by_ase=False
        )

        cmd_supercell(args)

        expected_file = outdir / "Au1Cu3-supercell-2x2x2.vasp"
        self.assertTrue(expected_file.exists(), f"Supercell file {expected_file} should be created")

    def test_cmd_supercell_with_ase(self):
        """Test supercell command with ASE."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            poscar=str(self.poscar_file),
            outdir=str(outdir),
            factors=[2, 2, 2],
            by_ase=True
        )

        cmd_supercell(args)

        expected_file = outdir / "Au1Cu3-supercell-2x2x2.vasp"
        self.assertTrue(expected_file.exists(), f"Supercell file {expected_file} should be created with ASE")

    def test_cmd_supercell_invalid_poscar(self):
        """Test supercell command with invalid POSCAR file."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            poscar=str(self.test_dir / "nonexistent.vasp"),
            outdir=str(outdir),
            factors=[2, 2, 2],
            by_ase=False
        )

        exit_code = cmd_supercell(args)
        self.assertEqual(exit_code, 1)

    def test_cmd_modeling(self):
        """Test modeling command."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            name="test_modeling",
            poscar=str(self.poscar_file),
            outdir=str(outdir),
            factors=[2, 2, 2],
            phase=None,
            config=None,
            seeds=[42],
            batch_size=1,
            enable_sqs=False,
            iterations=10000000
        )

        cmd_modeling(args)

        expected_dir = outdir / "test_modeling"
        self.assertTrue(expected_dir.exists(), f"Modeling directory {expected_dir} should be created")

        vasp_files = list(expected_dir.glob("*.vasp"))
        self.assertGreater(len(vasp_files), 0, "At least one VASP file should be generated")

    def test_cmd_modeling_with_seeds(self):
        """Test modeling command with multiple seeds."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            name="test_modeling_seeds",
            poscar=str(self.poscar_file),
            outdir=str(outdir),
            factors=[2, 2, 2],
            phase=None,
            config=None,
            seeds=[1, 2, 3],
            batch_size=1,
            enable_sqs=False,
            iterations=10000000
        )

        cmd_modeling(args)

        expected_dir = outdir / "test_modeling_seeds"
        self.assertTrue(expected_dir.exists(), f"Modeling directory {expected_dir} should be created")

        vasp_files = list(expected_dir.glob("*.vasp"))
        self.assertGreater(len(vasp_files), 0, "At least one VASP file should be generated")

    def test_cmd_modeling_missing_poscar_and_config(self):
        """Test modeling command with missing POSCAR and config."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            name="test_modeling",
            poscar=None,
            outdir=str(outdir),
            factors=[2, 2, 2],
            phase=None,
            config=None,
            seeds=[42],
            batch_size=1,
            enable_sqs=False,
            iterations=10000000
        )

        exit_code = cmd_modeling(args)
        self.assertEqual(exit_code, 1)

    def test_cmd_modeling_invalid_poscar(self):
        """Test modeling command with invalid POSCAR file."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            name="test_modeling",
            poscar=str(self.test_dir / "nonexistent.vasp"),
            outdir=str(outdir),
            factors=[2, 2, 2],
            phase=None,
            config=None,
            seeds=[42],
            batch_size=1,
            enable_sqs=False,
            iterations=10000000
        )

        exit_code = cmd_modeling(args)
        self.assertEqual(exit_code, 1)

    @patch('sys.argv', ['poscarkit', 'help'])
    def test_main_help_command(self):
        """Test main function with help command."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            main()
            output = fake_out.getvalue()
            self.assertIn("POSCARKIT", output)

    @patch('sys.argv', ['poscarkit'])
    def test_main_no_command(self):
        """Test main function with no command."""
        with patch('sys.stdout', new=io.StringIO()) as fake_out:
            exit_code = main()
            self.assertEqual(exit_code, 0)
            output = fake_out.getvalue()
            self.assertIn("usage:", output)

    def test_cmd_compare_identical(self):
        """Test compare command with identical structures."""
        args = argparse.Namespace(
            poscar1=str(self.poscar_file),
            poscar2=str(self.poscar_file)
        )

        cmd_compare(args)

    def test_cmd_compare_different(self):
        """Test compare command with different structures."""
        args = argparse.Namespace(
            poscar1=str(self.poscar_file),
            poscar2=str(self.poscar_file2)
        )

        cmd_compare(args)

    def test_cmd_compare_invalid_poscar1(self):
        """Test compare command with invalid first POSCAR file."""
        args = argparse.Namespace(
            poscar1=str(self.test_dir / "nonexistent1.vasp"),
            poscar2=str(self.poscar_file)
        )

        exit_code = cmd_compare(args)
        self.assertEqual(exit_code, 1)

    def test_cmd_compare_invalid_poscar2(self):
        """Test compare command with invalid second POSCAR file."""
        args = argparse.Namespace(
            poscar1=str(self.poscar_file),
            poscar2=str(self.test_dir / "nonexistent2.vasp")
        )

        exit_code = cmd_compare(args)
        self.assertEqual(exit_code, 1)

    def test_cmd_merge(self):
        """Test merge command."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            poscar1=str(self.poscar_file),
            poscar2=str(self.poscar_file3),
            outdir=str(outdir)
        )

        cmd_merge(args)

        expected_file = outdir / "POSCAR-merged-POSCAR-POSCAR3.vasp"
        self.assertTrue(expected_file.exists(), f"Merged file {expected_file} should be created")

    def test_cmd_merge_invalid_poscar1(self):
        """Test merge command with invalid first POSCAR file."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            poscar1=str(self.test_dir / "nonexistent1.vasp"),
            poscar2=str(self.poscar_file),
            outdir=str(outdir)
        )

        exit_code = cmd_merge(args)
        self.assertEqual(exit_code, 1)

    def test_cmd_merge_invalid_poscar2(self):
        """Test merge command with invalid second POSCAR file."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            poscar1=str(self.poscar_file),
            poscar2=str(self.test_dir / "nonexistent2.vasp"),
            outdir=str(outdir)
        )

        exit_code = cmd_merge(args)
        self.assertEqual(exit_code, 1)

    def test_cmd_separate_by_note(self):
        """Test separate command by note."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            poscar=str(self.poscar_file),
            key="note",
            outdir=str(outdir)
        )

        cmd_separate(args)

        group1_file = outdir / "POSCAR-group-1a-Au.vasp"
        group2_file = outdir / "POSCAR-group-3c-Cu.vasp"
        self.assertTrue(group1_file.exists(), f"Group file {group1_file} should be created")
        self.assertTrue(group2_file.exists(), f"Group file {group2_file} should be created")

    def test_cmd_separate_by_symbol(self):
        """Test separate command by symbol."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            poscar=str(self.poscar_file),
            key="symbol",
            outdir=str(outdir)
        )

        cmd_separate(args)

        au_file = outdir / "POSCAR-group-Au.vasp"
        cu_file = outdir / "POSCAR-group-Cu.vasp"
        self.assertTrue(au_file.exists(), f"Group file {au_file} should be created")
        self.assertTrue(cu_file.exists(), f"Group file {cu_file} should be created")

    def test_cmd_separate_invalid_poscar(self):
        """Test separate command with invalid POSCAR file."""
        outdir = self.test_dir / "output"
        args = argparse.Namespace(
            poscar=str(self.test_dir / "nonexistent.vasp"),
            key="note",
            outdir=str(outdir)
        )

        exit_code = cmd_separate(args)
        self.assertEqual(exit_code, 1)


if __name__ == "__main__":
    unittest.main()
