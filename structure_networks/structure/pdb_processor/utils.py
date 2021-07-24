import os
import warnings
import tempfile
from Bio.PDB import PDBIO
from Bio.PDB import PDBParser
from .pdb_structure import CaTrace
from .trajectory_reader import read_pdb
from Bio.PDB.Structure import Structure
from .pdb_structure import PDBStructure

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['convert_to_biopython',
           'convert_from_biopython',
           'mark_by_bfactor']


def convert_to_biopython(structure, **kwargs):
    id = kwargs.get('id', 'X')
    assert isinstance(structure, (PDBStructure, CaTrace))
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_file = os.path.join(tmp_dir, 'pdb_file.pdb')
        with open(tmp_file, 'w') as fh:
            structure.write(fh)
        with warnings.catch_warnings():
            pdb = PDBParser().get_structure(id, tmp_file)
    return pdb


def convert_from_biopython(structure):
    assert isinstance(structure, Structure)
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_file = os.path.join(tmp_dir, 'pdb_file.pdb')
        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(file=tmp_file)
        pdb = read_pdb(pdb_file=tmp_file)
    return pdb


def mark_by_bfactor(pdb_structure, residues, **kwargs):
    mark = kwargs.get('mark', 10)
    other = kwargs.get('other', 0)
    skip_other = kwargs.get('skip_other', False)
    assert isinstance(pdb_structure, (PDBStructure, CaTrace))
    if kwargs.get('inplace', False):
        structure = pdb_structure.copy()
    else:
        structure = pdb_structure
    residues = list(residues) if hasattr(residues, '__iter__') else [residues]
    all_residues = structure.residue_ids
    residues = [r for r in residues if r in all_residues]
    if not skip_other:
        for r in all_residues:
            if isinstance(structure, PDBStructure):
                for atom in structure.atom_list(residue_id=r):
                    structure.b_factor(residue_id=r, atom_name=atom, value=other)
            elif isinstance(structure, CaTrace):
                structure.b_factor(resid=r, value=other)

    for r in residues:
        if isinstance(structure, PDBStructure):
            for atom in structure.atom_list(residue_id=r):
                structure.b_factor(residue_id=r, atom_name=atom, value=mark)
        elif isinstance(structure, CaTrace):
            structure.b_factor(resid=r, value=mark)
    return structure

