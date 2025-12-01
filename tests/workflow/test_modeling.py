import os
import tempfile
import unittest
from pathlib import Path

from src.workflow.modeling import run_modeling

poscar = """SYSTEM=Au1Cu3-FCC
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

structure_info = {
    "cell": [[3.774, 0.0, 0.0], [0.0, 3.774, 0.0], [0.0, 0.0, 3.774]],
    "1a": {"atoms": ("Au", [[0.0, 0.0, 0.0]]), "sofs": {"V": 1.0}},
    "3c": {
        "atoms": ("Cu", [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
        "sofs": {"Co": 4.444e-1, "Ni": 0.4444, "V": 0.1112},
    },
}


class TestModelingWorkflow(unittest.TestCase):
    """Test modeling workflow."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test outputs
        # self.test_dir = Path(tempfile.mkdtemp())  # Using actual temp directory now
        self.test_dir = Path("output")
        self.output_dir = self.test_dir.joinpath("output")
        self.output_dir.mkdir(exist_ok=True, parents=True)  # Ensure parents are created

        # Write POSCAR content to a temporary file for testing
        self.poscar_file = self.test_dir / "POSCAR.vasp"
        with open(self.poscar_file, "w") as f:
            f.write(poscar)

    def test_run_modeling_with_structure_info(self):
        """Test run_modeling function with structure_info parameter."""
        # Define a simple cubic structure info for testing

        # Run modeling with simple supercell factors
        output_files = run_modeling(
            name="test-shuffle",
            poscar=Path(""),  # Empty filepath to use structure_info
            outdir=self.output_dir,
            supercell_factors=(10, 10, 10),  # Smaller supercell for faster testing
            structure_info=structure_info,
            shuffle_seeds=[42],  # Fixed seed for reproducible test
            batch_size=1,
        )

        # Check that files were generated
        self.assertTrue(len(output_files) > 0, "Should generate at least one file")

        # Check that output files exist
        for filepath in output_files:
            self.assertTrue(filepath.exists(), f"File {filepath} should exist")

        # Check that files are saved in the correct directory
        for filepath in output_files:
            self.assertTrue(
                str(self.output_dir) in str(filepath.parent),
                f"File {filepath} should be in the correct output directory",
            )

        print(
            f"Generated {len(output_files)} structure files from POSCAR in {self.output_dir}"
        )

    def test_run_modeling_with_poscar_file(self):
        """Test run_modeling function with POSCAR file parameter."""
        output_files = run_modeling(
            name="test-poscar",
            poscar=self.poscar_file,
            outdir=self.output_dir,
            supercell_factors=(2, 2, 2),  # Smaller supercell for faster testing
            structure_info=structure_info,
            shuffle_seeds=[42],  # Fixed seed for reproducible test
            batch_size=1,
        )

        # Check that files were generated
        self.assertTrue(len(output_files) > 0, "Should generate at least one file")

        # Check that output files exist
        for filepath in output_files:
            self.assertTrue(filepath.exists(), f"File {filepath} should exist")

        print(f"Generated {len(output_files)} structure files in {self.output_dir}")

    def test_run_modeling_with_sqsgen(self):
        """Test run_modeling function with sqsgen parameter."""
        try:
            output_files = run_modeling(
                name="test-sqsgen",
                poscar=self.poscar_file,
                outdir=self.output_dir,
                supercell_factors=(2, 2, 2),  # Smaller supercell for faster testing
                structure_info=structure_info,
                shuffle_seeds=[42],  # Fixed seed for reproducible test
                batch_size=1,
                enable_sqs=True,
            )

            # Check that files were generated
            self.assertTrue(len(output_files) > 0, "Should generate at least one file")

            # Check that output files exist
            for filepath in output_files:
                self.assertTrue(filepath.exists(), f"File {filepath} should exist")

            print(
                f"Generated {len(output_files)} structure files with SQS in {self.output_dir}"
            )
        except ImportError:
            # sqsgenerator may not be installed, skip test
            self.skipTest("sqsgenerator not installed")

    def tearDown(self):
        """Clean up test fixtures."""
        # Clean up the temporary directory and all its contents
        import shutil

        shutil.rmtree(self.test_dir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
