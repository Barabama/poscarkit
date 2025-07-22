
import logging
import os

from PoscarTools.AtomSupercell import supercell2file
from PoscarTools.AtomAllocate import allocate2files
from PoscarTools.AtomCountCN import countCN2files
from PoscarTools.AtomSlice import slice2file


logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

outdir = "examples/out"
os.makedirs(outdir, exist_ok=True)


def test_supercell():
    filepath = "examples/POSCAR-fcc.vasp"
    factors = (3, 3, 3)
    supercell2file(filepath=filepath, outdir=outdir, factors=factors)


def test_allocate():
    filepath = "examples/out/POSCAR-supercell-3x3x3.vasp"
    factors = (3, 3, 3)
    info = {"1a": {"atoms": ["Au", [[0, 0, 0]]], "sofs": {"V": 1}},
            "3c": {"atoms": ["Cu", [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]],
                   "sofs": {"Co": 0.4444, "Ni": 0.4444, "V": 0.1112}}}
    # info = None
    allocate2files(filepath=filepath, outdir=outdir, factors=factors, struct_info=info)


def test_countCN():
    filepath = "examples/out/POSCAR-allocate1-CoNiV.vasp"
    countCN2files(filepath=filepath, outdir=outdir)


def test_slice():
    filepath = "examples/out/POSCAR-allocate1-CoNiV.vasp"
    direction = (1, 1, 1)
    slice2file(filepath=filepath, outdir=outdir, miller_index=direction)


if __name__ == "__main__":
    test_supercell()
    test_allocate()
    test_countCN()
    test_slice()
