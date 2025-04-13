


import logging

from PoscarTools.Supercell import supercell2file
from PoscarTools.AtomSlice import slice2file
from PoscarTools.AtomShuffle import shuffle2files
from PoscarTools.AtomAllocate import allocate2file
from PoscarTools.AtomCountCN import countCN2files


logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def test_supercell():
    filepath = "examples/POSCAR-fcc.vasp"
    factors = (2, 2, 2)
    supercell2file(filepath, factors)
    # make_supercell(filepath, factors)
    
def test_slice():
    filepath = "examples/POSCAR-fcc.vasp"
    direction = (1, 1, 0)
    slice2file(filepath, direction)
    

def test_shuffle():
    filepath = "examples/POSCAR-fcc-222.vasp"
    structure = {"cell": [3.5,3.5,3.5],
                 "1a": {"atoms": ["Au", [[0,0,0]]]},
                 "3c": {"atoms": ["Cu", [[0,0.5,0.5],[0.5,0,0.5], [0.5,0.5,0]]]}}
    seeds = [7, 42]
    shuffle2files(filepath, structure, seeds)


def test_allocate():
    filepath = "examples/POSCAR-fcc-222.vasp"
    structure = {"cell": [3.5,3.5,3.5],
                 "1a": {"atoms": ["Au", [[0,0,0]]], "sofs": {"V": 1}},
                 "3c": {"atoms": ["Cu", [[0,0.5,0.5],[0.5,0,0.5], [0.5,0.5,0]]],
                        "sofs": {"Co": 0.4444, "Cr": 0.4444, "V": 0.1112}}}
    allocate2file(filepath, structure, (2,2,2))

if __name__ == "__main__":
    test_supercell()
    # test_slice()
    test_shuffle()
    test_allocate()
    