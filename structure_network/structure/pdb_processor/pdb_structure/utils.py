import logging
from copy import deepcopy
from .catrace import CaTrace
from .pdb_structure import PDBStructure

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['do_mutation', 'pdb_to_catrace']


def do_mutation(ca_trace, res_id, to_residue):
    assert isinstance(ca_trace, CaTrace)
    residue_ids = ca_trace.residue_ids
    assert res_id in residue_ids
    structure_copy = ca_trace.copy()
    structure_copy.set_amino(res_id, to_residue)
    return structure_copy


def pdb_to_catrace(pdb_struct):
    assert isinstance(pdb_struct, PDBStructure)
    residue_ids = pdb_struct.residue_ids
    entries = list()
    logger = logging.getLogger(name='pdb_processor.pdb_to_catrace')
    for r in residue_ids:
        resname = pdb_struct.residue_name(residue_id=r)
        x, y, z = pdb_struct.xyz(residue_id=r, atom_name='CA')
        bfactor = pdb_struct.b_factor(residue_id=r, atom_name='CA')
        if (x is None) or (y is None) or (z is None):
            logger.error('Missing CA in residue (%d)' % r)
            raise Exception('Missing CA in the protein for residue (%d)' % r)
        entries.append({'resname': resname,
                        'resid': r,
                        'x': x,
                        'y': y,
                        'z': z,
                        'bfactor': bfactor})
    return CaTrace(pdb_struct.name, entries, pdb_struct.chain)

