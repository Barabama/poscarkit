import tempfile
import unittest
import shutil
from pathlib import Path

from src.modeling.base import SimplePoscar
from src.workflow.import_to_model import run_import_to_model

# Sample Pandat export uses all-caps element symbols in headers:
#   x(CO), y(CO#1@FCC), ...  representing Co, Cu, Ni.
_SAMPLE = (
    Path(__file__).resolve().parents[2]
    / "ref_sof_data"
    / "pandat-data.csv"
)


class TestImportToModelSymbolCase(unittest.TestCase):
    """Regression: Excel-imported POSCARs must use canonical element symbols.

    https://github.com/Barabama/poscarkit issue: POSCARs generated from
    Pandat/Thermo-Calc Excel data carried all-caps symbols ('CO', 'CU',
    'NI'), which made ASE backends (countcn/supercell/surface/slice) raise
    KeyError.  The generated POSCAR symbol line must be canonical.
    """

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.test_dir, ignore_errors=True)

    def test_generated_poscar_has_canonical_symbols(self):
        output = run_import_to_model(
            csv_path=str(_SAMPLE),
            config_path="config.toml",
            phase="FCC",
            temperatures=[1000.0],
            factors=(2, 2, 2),
            seeds=[42],
            batch_size=1,
            outdir=str(self.test_dir),
            name="symcase",
        )
        self.assertTrue(output, "expected generated POSCAR files")

        # The merged shuffle POSCAR (filename without '-site').
        mains = [
            p for p in output
            if "shuffle1.vasp" in p.name and "-site" not in p.name
        ]
        self.assertEqual(len(mains), 1, f"expected one merged POSCAR in {output}")
        poscar = mains[0]

        # Symbol line is the 6th line of the file.
        sym_tokens = set(poscar.read_text().splitlines()[5].split())
        self.assertEqual(sym_tokens, {"Co", "Cu", "Ni"})

        # And it must round-trip through the ASE bridge without error.
        struct = SimplePoscar.read_poscar(poscar)
        atoms = SimplePoscar.struct2atoms(struct)
        self.assertEqual(set(atoms.get_chemical_symbols()), {"Co", "Cu", "Ni"})


if __name__ == "__main__":
    unittest.main()
