"""Unit tests for src/gui/app.py and src/gui/forms.py."""

import argparse
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch, MagicMock


# ================================================================== #
#  Pure logic tests — no Tk dependency                               #
# ================================================================== #


class TestTomlValue(unittest.TestCase):
    def setUp(self):
        from src.gui.app import _toml_value
        self.f = _toml_value

    def test_bool_true(self):
        self.assertEqual(self.f(True), "true")

    def test_bool_false(self):
        self.assertEqual(self.f(False), "false")

    def test_int(self):
        self.assertEqual(self.f(42), "42")

    def test_float(self):
        self.assertEqual(self.f(3.14), "3.14")

    def test_string(self):
        self.assertEqual(self.f("hello"), '"hello"')

    def test_empty_string(self):
        self.assertEqual(self.f(""), '""')

    def test_list_of_ints(self):
        self.assertEqual(self.f([1, 2, 3]), "[1, 2, 3]")

    def test_list_of_strings(self):
        self.assertEqual(self.f(["a", "b"]), '["a", "b"]')

    def test_list_of_bools(self):
        self.assertEqual(self.f([True, False]), "[true, false]")

    def test_list_mixed(self):
        result = self.f([3, "x", True])
        self.assertEqual(result, '[3, "x", true]')

    def test_empty_list(self):
        self.assertEqual(self.f([]), "[]")


class TestColorForName(unittest.TestCase):
    def setUp(self):
        from src.gui.app import _color_for_name
        self.f = _color_for_name

    def test_known_name(self):
        self.assertEqual(self.f("Modeling"), "#4A90D9")
        self.assertEqual(self.f("Config"), "#7F8C8D")

    def test_unknown_name(self):
        self.assertEqual(self.f("NonExistent"), "#34495E")


class TestParseTemps(unittest.TestCase):
    def setUp(self):
        from src.gui.forms import _parse_temps
        self.f = _parse_temps

    def test_empty(self):
        self.assertEqual(self.f(""), [])
        self.assertEqual(self.f("   "), [])

    def test_space_separated(self):
        self.assertEqual(self.f("300 600 900"), [300.0, 600.0, 900.0])

    def test_comma_separated(self):
        self.assertEqual(self.f("300,600,900"), [300.0, 600.0, 900.0])

    def test_chinese_comma(self):
        self.assertEqual(self.f("300，600，900"), [300.0, 600.0, 900.0])

    def test_semicolon(self):
        self.assertEqual(self.f("300;600;900"), [300.0, 600.0, 900.0])

    def test_k_suffix(self):
        self.assertEqual(self.f("300K 600k 900"), [300.0, 600.0, 900.0])

    def test_mixed_separators(self):
        self.assertEqual(self.f("300K,600；900、1200"), [300.0, 600.0, 900.0, 1200.0])

    def test_float_values(self):
        self.assertEqual(self.f("298.15 773.15"), [298.15, 773.15])


class TestParseSeeds(unittest.TestCase):
    def setUp(self):
        from src.gui.forms import _parse_seeds
        self.f = _parse_seeds

    def test_empty(self):
        self.assertIsNone(self.f(""))
        self.assertIsNone(self.f("   "))

    def test_space_separated(self):
        self.assertEqual(self.f("42 123 7"), [42, 123, 7])

    def test_comma_separated(self):
        self.assertEqual(self.f("42,123,7"), [42, 123, 7])

    def test_chinese_comma(self):
        self.assertEqual(self.f("42，123，7"), [42, 123, 7])

    def test_semicolon(self):
        self.assertEqual(self.f("42;123"), [42, 123])

    def test_single_value(self):
        self.assertEqual(self.f("42"), [42])

    def test_mixed_separators(self):
        self.assertEqual(self.f("42,123 7；99"), [42, 123, 7, 99])


class TestWriteSofsToConfig(unittest.TestCase):
    def setUp(self):
        from src.gui.forms import _write_sofs_to_config
        self.f = _write_sofs_to_config
        self.test_dir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_replace_existing_section(self):
        cfg = self.test_dir / "config.toml"
        cfg.write_text(
            "[FCC.1a.sofs]\nAu = 1.0\n\n[FCC.3c.sofs]\nCo = 0.5\nNi = 0.5\n",
            encoding="utf-8",
        )
        self.f(str(cfg), "FCC", {"1a": {"Ag": 0.6, "Au": 0.4}})
        content = cfg.read_text(encoding="utf-8")
        self.assertIn("Ag = 0.6", content)
        self.assertIn("Au = 0.4", content)
        self.assertIn("[FCC.1a.sofs]", content)

    def test_append_new_section(self):
        cfg = self.test_dir / "config.toml"
        cfg.write_text("[FCC.1a.sofs]\nAu = 1.0\n", encoding="utf-8")
        self.f(str(cfg), "FCC", {"1a": {"Au": 1.0}, "3c": {"Cu": 0.5, "Ni": 0.5}})
        content = cfg.read_text(encoding="utf-8")
        self.assertIn("[FCC.3c.sofs]", content)
        self.assertIn("Cu = 0.5", content)

    def test_preserves_other_content(self):
        cfg = self.test_dir / "config.toml"
        cfg.write_text(
            "name = \"test\"\n\n[FCC.1a.sofs]\nAu = 1.0\n\n[BCC]\ncell = [2.7]\n",
            encoding="utf-8",
        )
        self.f(str(cfg), "FCC", {"1a": {"Ag": 1.0}})
        content = cfg.read_text(encoding="utf-8")
        self.assertIn('name = "test"', content)
        self.assertIn("[BCC]", content)


class TestSaveToConfig(unittest.TestCase):
    """Test the file-writing logic of _save_to_config."""

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())
        self.orig_cwd = os.getcwd()
        os.chdir(self.test_dir)

    def tearDown(self):
        os.chdir(self.orig_cwd)
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def _make_gui_stub(self):
        """Create a minimal stub with _save_to_config and _load_config."""
        from src.gui.app import PoscaKitGUI
        stub = object.__new__(PoscaKitGUI)
        stub._cfg = {}
        stub.root = MagicMock()
        return stub

    def test_creates_new_file(self):
        stub = self._make_gui_stub()
        args = argparse.Namespace(
            name="test", poscar="a.vasp", outdir="out",
            phase="FCC", factors=[3, 3, 3], seeds=None,
            batch_size=2, enable_sqs=False, iterations=1000000,
            cutoff_mult=1.1, parallel=4, by_ase=False, pbc=True,
            slice_direction=[0, 0, 1],
        )
        stub._save_to_config(args)
        cfg_path = Path("config.toml")
        self.assertTrue(cfg_path.is_file())
        content = cfg_path.read_text(encoding="utf-8")
        self.assertIn('name = "test"', content)
        self.assertIn("batch_size = 2", content)
        self.assertIn("pbc = true", content)

    def test_replaces_existing_key(self):
        Path("config.toml").write_text(
            'name = "old"\nposcar = "old.vasp"\n\n[FCC]\ncell = [3.7]\n',
            encoding="utf-8",
        )
        stub = self._make_gui_stub()
        args = argparse.Namespace(
            name="new", poscar="new.vasp", outdir=None,
            phase=None, factors=None, seeds=None,
            batch_size=None, enable_sqs=None, iterations=None,
            cutoff_mult=None, parallel=None, by_ase=None, pbc=None,
            slice_direction=None,
        )
        stub._save_to_config(args)
        content = Path("config.toml").read_text(encoding="utf-8")
        self.assertIn('name = "new"', content)
        self.assertIn('poscar = "new.vasp"', content)
        self.assertNotIn('"old"', content)

    def test_uncomments_key(self):
        Path("config.toml").write_text(
            '# name = "MyWork"\n# poscar = "./a.vasp"\n\n[FCC]\ncell = [3.7]\n',
            encoding="utf-8",
        )
        stub = self._make_gui_stub()
        args = argparse.Namespace(
            name="updated", poscar=None, outdir=None,
            phase=None, factors=None, seeds=None,
            batch_size=None, enable_sqs=None, iterations=None,
            cutoff_mult=None, parallel=None, by_ase=None, pbc=None,
            slice_direction=None,
        )
        stub._save_to_config(args)
        content = Path("config.toml").read_text(encoding="utf-8")
        self.assertIn('name = "updated"', content)
        self.assertNotIn("# name", content)

    def test_inserts_before_section(self):
        Path("config.toml").write_text(
            "# comment\n\n[FCC]\ncell = [3.7]\n",
            encoding="utf-8",
        )
        stub = self._make_gui_stub()
        args = argparse.Namespace(
            name="work", poscar=None, outdir=None,
            phase=None, factors=None, seeds=None,
            batch_size=None, enable_sqs=None, iterations=None,
            cutoff_mult=None, parallel=None, by_ase=None, pbc=None,
            slice_direction=None,
        )
        stub._save_to_config(args)
        content = Path("config.toml").read_text(encoding="utf-8")
        # name should appear before [FCC]
        name_pos = content.index('name = "work"')
        fcc_pos = content.index("[FCC]")
        self.assertLess(name_pos, fcc_pos)

    def test_roundtrip_save_load(self):
        """Save then load — values should be recoverable."""
        import tomllib
        from src.config import normalize_config_keys

        stub = self._make_gui_stub()
        args = argparse.Namespace(
            name="roundtrip", poscar="test.vasp", outdir="output",
            phase="BCC", factors=[2, 2, 2], seeds=[42, 7],
            batch_size=3, enable_sqs=True, iterations=5000000,
            cutoff_mult=1.2, parallel=8, by_ase=True, pbc=False,
            slice_direction=[1, 1, 0],
        )
        stub._save_to_config(args)

        with open("config.toml", "rb") as f:
            cfg = normalize_config_keys(tomllib.load(f))

        self.assertEqual(cfg["name"], "roundtrip")
        self.assertEqual(cfg["poscar"], "test.vasp")
        self.assertEqual(cfg["supercell_factors"], [2, 2, 2])
        self.assertEqual(cfg["shuffle_seeds"], [42, 7])
        self.assertEqual(cfg["batch_size"], 3)
        self.assertEqual(cfg["enable_sqs"], True)
        self.assertEqual(cfg["by_ase"], True)
        self.assertEqual(cfg["pbc"], False)
        self.assertEqual(cfg["slice_direction"], [1, 1, 0])


# ================================================================== #
#  Widget integration tests — require Tk display                     #
# ================================================================== #

try:
    import tkinter as tk
    _test_root = tk.Tk()
    _test_root.withdraw()
    _test_root.destroy()
    _HAS_TK = True
except Exception:
    _HAS_TK = False


@unittest.skipUnless(_HAS_TK, "No display available for Tk tests")
class TestFormWidgets(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = tk.Tk()
        cls.root.withdraw()

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()

    def setUp(self):
        self.frame = tk.Frame(self.root)

    def tearDown(self):
        self.frame.destroy()

    def test_entry_default(self):
        from src.gui.forms import _entry
        var = _entry(self.frame, "Test", {}, "key", "default_val")
        self.assertEqual(var.get(), "default_val")

    def test_entry_from_cfg(self):
        from src.gui.forms import _entry
        var = _entry(self.frame, "Test", {"key": "from_cfg"}, "key", "default")
        self.assertEqual(var.get(), "from_cfg")

    def test_entry_list_value_joined(self):
        from src.gui.forms import _entry
        var = _entry(self.frame, "Test", {"seeds": [42, 7, 123]}, "seeds", "")
        self.assertEqual(var.get(), "42 7 123")

    def test_checkbox_default_false(self):
        from src.gui.forms import _checkbox
        var = _checkbox(self.frame, "Test", False)
        self.assertFalse(var.get())

    def test_checkbox_default_true(self):
        from src.gui.forms import _checkbox
        var = _checkbox(self.frame, "Test", True)
        self.assertTrue(var.get())

    def test_combo_default(self):
        from src.gui.forms import _combo
        _, var = _combo(self.frame, "Phase", ["FCC", "BCC"], {}, "phase", "FCC")
        self.assertEqual(var.get(), "FCC")

    def test_combo_from_cfg(self):
        from src.gui.forms import _combo
        _, var = _combo(self.frame, "Phase", ["FCC", "BCC"], {"phase": "BCC"}, "phase", "FCC")
        self.assertEqual(var.get(), "BCC")

    def test_int_entries_defaults(self):
        from src.gui.forms import _int_entries
        vars_ = _int_entries(self.frame, "Factors", 3, {}, "factors", (3, 3, 3))
        self.assertEqual(len(vars_), 3)
        self.assertEqual([v.get() for v in vars_], ["3", "3", "3"])

    def test_int_entries_from_cfg(self):
        from src.gui.forms import _int_entries
        vars_ = _int_entries(self.frame, "Factors", 3, {"factors": [2, 4, 6]}, "factors", (3, 3, 3))
        self.assertEqual([v.get() for v in vars_], ["2", "4", "6"])


@unittest.skipUnless(_HAS_TK, "No display available for Tk tests")
class TestFormBuilders(unittest.TestCase):
    """Test that form builders produce valid get_args() callables."""

    @classmethod
    def setUpClass(cls):
        cls.root = tk.Tk()
        cls.root.withdraw()
        cls.test_dir = Path(tempfile.mkdtemp())
        cls.orig_cwd = os.getcwd()
        os.chdir(cls.test_dir)
        # Write a minimal config.toml for forms that need it
        from src.config import DEFAULT_CONFIG
        (cls.test_dir / "config.toml").write_text(DEFAULT_CONFIG, encoding="utf-8")

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.orig_cwd)
        shutil.rmtree(cls.test_dir, ignore_errors=True)
        cls.root.destroy()

    def setUp(self):
        self.frame = tk.Frame(self.root)

    def tearDown(self):
        self.frame.destroy()

    def test_countcn_form(self):
        from src.gui.forms import countcn_form
        get_args = countcn_form(self.frame, {"name": "test_cn", "poscar": "a.vasp"})
        args = get_args()
        self.assertEqual(args.name, "test_cn")
        self.assertEqual(args.poscar, "a.vasp")
        self.assertIsInstance(args.cutoff_mult, float)
        self.assertIsInstance(args.parallel, int)
        self.assertIsInstance(args.by_ase, bool)
        self.assertIsInstance(args.pbc, bool)

    def test_slice_form(self):
        from src.gui.forms import slice_form
        get_args = slice_form(self.frame, {"name": "test_slice"})
        args = get_args()
        self.assertEqual(args.name, "test_slice")
        self.assertEqual(len(args.miller_index), 3)

    def test_supercell_form(self):
        from src.gui.forms import supercell_form
        get_args = supercell_form(self.frame, {"poscar": "cell.vasp"})
        args = get_args()
        self.assertEqual(args.poscar, "cell.vasp")
        self.assertEqual(len(args.factors), 3)
        self.assertIsInstance(args.by_ase, bool)

    def test_compare_form(self):
        from src.gui.forms import compare_form
        get_args = compare_form(self.frame, {})
        args = get_args()
        self.assertHasAttr(args, "poscar1")
        self.assertHasAttr(args, "poscar2")

    def test_separate_form(self):
        from src.gui.forms import separate_form
        get_args = separate_form(self.frame, {"outdir": "sep_out"})
        args = get_args()
        self.assertEqual(args.outdir, "sep_out")
        self.assertIn(args.key, ["note", "symbol", "x", "y", "z"])

    def test_merge_form(self):
        from src.gui.forms import merge_form
        get_args = merge_form(self.frame, {})
        args = get_args()
        self.assertIsInstance(args.poscars, list)
        self.assertIsInstance(args.outdir, str)

    def assertHasAttr(self, obj, attr):
        self.assertTrue(hasattr(obj, attr), f"Missing attribute: {attr}")


@unittest.skipUnless(_HAS_TK, "No display available for Tk tests")
class TestLoadConfig(unittest.TestCase):
    """Test PoscaKitGUI._load_config behavior."""

    @classmethod
    def setUpClass(cls):
        cls.root = tk.Tk()
        cls.root.withdraw()

    @classmethod
    def tearDownClass(cls):
        cls.root.destroy()

    def setUp(self):
        self.test_dir = Path(tempfile.mkdtemp())
        self.orig_cwd = os.getcwd()
        os.chdir(self.test_dir)

    def tearDown(self):
        os.chdir(self.orig_cwd)
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def _make_stub(self):
        from src.gui.app import PoscaKitGUI
        stub = object.__new__(PoscaKitGUI)
        stub._cfg = {}
        stub.root = MagicMock()
        return stub

    def test_generates_default_when_missing(self):
        stub = self._make_stub()
        stub._load_config()
        self.assertTrue(Path("config.toml").is_file())
        self.assertIn("FCC", stub._cfg)

    def test_loads_existing_file(self):
        Path("config.toml").write_text(
            'name = "loaded"\n\n[FCC]\ncell = [3.7, 3.7, 3.7]\n',
            encoding="utf-8",
        )
        stub = self._make_stub()
        stub._load_config()
        self.assertEqual(stub._cfg.get("name"), "loaded")

    def test_handles_invalid_toml(self):
        Path("config.toml").write_text("invalid = [unclosed", encoding="utf-8")
        stub = self._make_stub()
        stub._load_config()
        self.assertEqual(stub._cfg, {})


if __name__ == "__main__":
    unittest.main()
