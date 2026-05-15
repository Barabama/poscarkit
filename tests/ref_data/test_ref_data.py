"""Tests for src.io — format detection, reading, and SOF extraction."""

import os
import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from src.io.readers import detect_format, read
from src.io.ir import (
    SOFData,
    PHASE_SITE_RATIOS,
    PHASE_SUBLATTICE_MAP,
    parse_sublattice_map,
    get_sofs_at,
    build_structure_info,
)


def _write_csv(dirpath, name, header, rows):
    import csv
    path = dirpath / name
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for row in rows:
            w.writerow(row)
    return str(path)


class TestDetectFormat(unittest.TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_detect_pandat_csv(self):
        path = _write_csv(self.tmpdir,
            "test.csv",
            ["T", "P", "phase_name", "x(CO)", "y(CO#1@FCC)", "y(CO#2@FCC)"],
            [["273", "1", "FCC", "0.33", "0.5", "0.5"]],
        )
        self.assertEqual(detect_format(path), "pandat")

    def test_detect_tc_exps_csv(self):
        path = _write_csv(self.tmpdir,
            "test.csv",
            ["T", "G(FCC)", "X(FCC,CO)", "Y(FCC,CO)", "Y(FCC,CO#2)"],
            [["1073", "-10000", "0.33", "0.5", "0.5"]],
        )
        self.assertEqual(detect_format(path), "tc_exps")

    def test_detect_teacher_csv(self):
        path = _write_csv(self.tmpdir,
            "test.csv",
            ["T", "G", "H", "S", "NP", "CO", "CO#2", "CR", "CR#2"],
            [["1073", "-10000", "5000", "10", "1", "0.5", "0.5", "0.5", "0.5"]],
        )
        self.assertEqual(detect_format(path), "teacher")


class TestParseSublatticeMap(unittest.TestCase):

    def test_single(self):
        self.assertEqual(parse_sublattice_map("1:1a"), {1: "1a"})

    def test_multiple(self):
        self.assertEqual(
            parse_sublattice_map("1:1a,2:3c"), {1: "1a", 2: "3c"}
        )

    def test_invalid(self):
        with self.assertRaises(ValueError):
            parse_sublattice_map("abc")


class TestGetSofsAt(unittest.TestCase):

    def setUp(self):
        self.ir = SOFData(
            source_path="test",
            phase="FCC",
            site_ratios=[0.25, 0.75],
            T=np.array([1073.0, 1273.0]),
            Y_subl=[
                {"CO": np.array([0.5, 0.6]), "CU": np.array([0.5, 0.4])},
                {"CR": np.array([0.3, 0.4]), "CU": np.array([0.7, 0.6])},
            ],
            composition={"CO": 0.2, "CU": 0.4, "CR": 0.4},
            elements=["CO", "CR", "CU"],
        )

    def test_exact_match(self):
        result = get_sofs_at(self.ir, 1073.0)
        self.assertIn("1a", result)
        self.assertAlmostEqual(result["1a"]["CO"], 0.5, places=4)

    def test_nearest_match(self):
        result = get_sofs_at(self.ir, 1080.0)
        self.assertAlmostEqual(result["1a"]["CO"], 0.5, places=4)

    def test_custom_map(self):
        result = get_sofs_at(self.ir, 1073.0, sublattice_map={1: "3c", 2: "1a"})
        self.assertIn("3c", result)
        self.assertAlmostEqual(result["3c"]["CO"], 0.5, places=4)


class TestBuildStructureInfo(unittest.TestCase):

    def test_merge(self):
        cfg = {
            "cell": [3.774, 3.774, 3.774],
            "1a": {"atoms": ["Au", [[0, 0, 0]]]},
            "3c": {"atoms": ["Cu", [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]]},
        }
        sofs = {
            "1a": {"CO": 0.5, "CU": 0.5},
            "3c": {"CR": 0.3, "CU": 0.7},
        }
        result = build_structure_info(cfg, sofs)
        self.assertEqual(result["cell"], [3.774, 3.774, 3.774])
        self.assertEqual(result["1a"]["sofs"], {"CO": 0.5, "CU": 0.5})
        self.assertEqual(result["3c"]["sofs"], {"CR": 0.3, "CU": 0.7})


class TestReadTeacher(unittest.TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_read(self):
        path = _write_csv(self.tmpdir,
            "teacher_test.csv",
            ["T", "G", "H", "S", "NP", "CO", "CO#2", "CR", "CR#2", "CU", "CU#2"],
            [
                ["73", "5000", "6000", "12", "1", "0.16", "0.0", "0.84", "0.213", "0.0", "0.253"],
                ["83", "5100", "6100", "13", "1", "0.15", "0.0", "0.85", "0.22", "0.0", "0.25"],
            ],
        )
        ir = read(str(path), phase_hint="FCC")
        self.assertEqual(ir.phase, "FCC")
        self.assertEqual(len(ir.T), 2)
        self.assertIn("CO", ir.elements)
        # Sublattice 1 has CO and CU
        self.assertAlmostEqual(float(ir.Y_subl[0]["CO"][0]), 0.16, places=4)


class TestReadTcExps(unittest.TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_read_csv(self):
        path = _write_csv(self.tmpdir,
            "tc_test.csv",
            ["T", "G(FCC)", "X(FCC,CO)", "X(FCC,CU)",
             "Y(FCC,CO)", "Y(FCC,CO#2)", "Y(FCC,CU)", "Y(FCC,CU#2)"],
            [["1073", "-10000", "0.33", "0.67", "0.5", "0.2", "0.5", "0.8"]],
        )
        ir = read(path, phase_hint="FCC")
        self.assertEqual(ir.phase, "FCC")
        self.assertAlmostEqual(float(ir.Y_subl[0]["CO"][0]), 0.5, places=4)
        self.assertAlmostEqual(float(ir.Y_subl[1]["CO"][0]), 0.2, places=4)


class TestReadPandat(unittest.TestCase):

    def setUp(self):
        self.tmpdir = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_read_csv(self):
        path = _write_csv(self.tmpdir,
            "pandat_test.csv",
            ["T", "P", "phase_name", "x(CO)", "x(CU)", "x(NI)",
             "y(CO#1@FCC)", "y(CU#1@FCC)", "y(NI#1@FCC)",
             "y(CO#2@FCC)", "y(CU#2@FCC)", "y(NI#2@FCC)"],
            [["773", "1", "FCC", "0.333", "0.333", "0.334",
              "0.1", "0.2", "0.7", "0.4", "0.4", "0.2"]],
        )
        ir = read(str(path))
        self.assertEqual(ir.phase, "FCC")
        self.assertAlmostEqual(float(ir.Y_subl[0]["CO"][0]), 0.1, places=4)
        self.assertAlmostEqual(float(ir.Y_subl[1]["CO"][0]), 0.4, places=4)


if __name__ == "__main__":
    unittest.main()
