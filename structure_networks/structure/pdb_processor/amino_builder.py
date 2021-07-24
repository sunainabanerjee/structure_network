import math
import numpy as np
from copy import deepcopy
from .pdb_structure import *
from structure_networks.geometry import *
from structure_networks.structure.amino_acids import *


__version__ = '1.0'
__all__ = ['reconstruction_units',
           'AminoBuilder',
           'ProteinReconstruction']

__backbone_model_desc__ = [[False, None, 1.46, 2.87, None],
                           [False, None, 1.53, 1.20, None],
                           [False, None, 1.23, 1.03, None]]

__backbone_dihedral__ = [['CA:-2', 'CA:-1', 'CA', 'N'],
                         ['CA:-1', 'N', 'CA', 'C'],
                         ['N', 'CA', 'C', 'O']]

__side_chain_model_desc__ = {
    'A': [[False, None, 1.53, 1.21, None]],
    'C': [[False, None, 1.53, 1.21, None], [False, None, 1.81, 1.15, None]],
    'D': [[False, None, 1.53, 1.21, None], [False, None, 1.53, 1.07, None],
          [False, None, 1.25, 1.07, None], [True, -1, 1.25, 1.07, math.pi]],
    'E': [[False, None, 1.53, 1.21, None], [False, None, 1.53, 1.21, None],
          [False, None, 1.53, 1.15, None], [False, None, 1.25, 1.07, None],
          [True, -1, 1.25, 1.07, math.pi]],
    'F': [[False, None, 1.53, 1.21, None], [False, None, 1.53, 1.15, None],
          [False, None, 1.39, 1.03, None], [True, -1, 1.39, 1.03, math.pi],
          [False, None, 1.39, 1.03, None], [True, -1, 1.39, 1.03, 0.],
          [False, None, 1.38, 1.05, 0.]],
    'G': [],
    'H': [[False, None, 1.53, 1.21, None], [False, None, 1.53, 1.15, None],
          [False, None, 1.39, 1.00, None], [True, -1, 1.35, 0.85, math.pi],
          [False, None, 1.32, 1.23, math.pi], [False, None, 1.37, 1.27, math.pi]],
    'I': [[False, None, 1.53, 1.21, None],  [False, None, 1.53, 1.21, None],
          [True, -1, 1.53, 1.21, 0.667 * math.pi], [False, None, 1.53, 1.21, None]],
    'K': [[False, None, 1.53, 1.21, None], [False, None, 1.53, 1.15, None],
          [False, None, 1.53, 1.19, None], [False, None, 1.53, 1.19, None],
          [False, None, 1.49, 1.19, None]],
    'L': [[False, None, 1.53, 1.21, None], [False, None, 1.53, 1.10, None],
          [False, None, 1.53, 1.21, None], [True, -1, 1.53, 1.21, 0.667*math.pi]],
    'M': [[False, None, 1.53, 1.21, None], [False, None, 1.53, 1.15, None],
          [False, None, 1.81, 1.18, None], [False, None, 1.79, 1.38, None]],
    'N': [[False, None, 1.53, 1.21, None], [False, None, 1.53, 1.17, None],
          [False, None, 1.23, 1.03, None], [True, -1, 1.32, 1.10, math.pi]],
    'P': [[False, None, 1.53, 1.34, None], [False, None, 1.50, 1.32, None],
          [False, None, 1.50, 1.30, None]],
    'Q': [[False, None, 1.53, 1.21, None], [False, None, 1.53, 1.15, None],
          [False, None, 1.52, 1.17, None], [False, None, 1.32, 1.10, None],
          [True, -1, 1.23, 1.03, math.pi]],
    'R': [[False, None, 1.53, 1.21, None], [False, None, 1.53, 1.15, None],
          [False, None, 1.53, 1.21, None], [False, None, 1.46, 1.19, None],
          [False, None, 1.33, 0.96, None], [False, None, 1.32, 1.04, 0.],
          [False, None, 1.32, 1.04, math.pi]],
    'S': [[False, None, 1.53, 1.21, None], [False, None, 1.41, 1.21, None]],
    'T': [[False, None, 1.54, 1.19, None], [False, None, 1.43, 1.23, None],
          [True, -1, 1.53, 1.20, 0.667*math.pi]],
    'V': [[False, None, 1.54, 1.19, None], [False, None, 1.53, 1.21, None],
          [True, -1, 1.53, 1.21, 0.667 * math.pi]],
    'W': [[False, None, 1.53, 1.21, None], [False, None, 1.50, 1.15, None],
          [False, None, 1.38, 0.92, None], [True, -1, 1.43, 0.93, math.pi],
          [False, None, 1.38, 1.22, math.pi], [False, None, 1.41, 1.27, math.pi],
          [False, None, 1.41, 0.81, 0], [False, None, 1.40, 1.00, math.pi],
          [False, None, 1.39, 1.07, math.pi], [False, None, 1.37, 1.09, 0]],
    'Y': [[False, None, 1.53, 1.21, None], [False, None, 1.51, 1.15, None],
          [False, None, 1.39, 1.03, None], [True, -1, 1.39, 1.03, math.pi],
          [False, None, 1.39, 1.03, math.pi], [False, None, 1.39, 1.03, math.pi],
          [False, None, 1.39, 1.05, 0], [False, None, 1.38, 1.05, math.pi]]
}


__side_chain_dihedral__ = {
    'A': [['CA:-1', 'N', 'CA', 'CB']],
    'C': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'SG']],
    'D': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'OD1'], ['CA', 'CB', 'CG', 'OD2']],
    'E': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'CD'],  ['CB', 'CG', 'CD', 'OE1'],
          ['CB', 'CG', 'CD', 'OE2']],
    'F': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'CD1'], ['CA', 'CB', 'CG', 'CD2'],
          ['CB', 'CG', 'CD1', 'CE1'], ['CB', 'CG', 'CD2', 'CE2'],
          ['CG', 'CD1', 'CE1', 'CZ']],
    'G': [],
    'H': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'ND1'], ['CA', 'CB', 'CG', 'CD2'],
          ['CB', 'CG', 'ND1', 'CE1'], ['CB', 'CG', 'CD2', 'NE2']],
    'I': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG1'],
          ['N', 'CA', 'CB', 'CG2'], ['CA', 'CB', 'CG1', 'CD1']],
    'K': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'CD'], ['CB', 'CG', 'CD', 'CE'],
          ['CG', 'CD', 'CE', 'NZ']],
    'L': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'CD1'], ['CA', 'CB', 'CG', 'CD2']],
    'M': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'SD'], ['CB', 'CG', 'SD', 'CE']],
    'N': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'OD1'], ['CA', 'CB', 'CG', 'ND2']],
    'P': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'CD']],
    'Q': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'CD'], ['CB', 'CG', 'CD', 'NE2'],
          ['CB', 'CG', 'CD', 'OE1']],
    'R': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'CD'], ['CB', 'CG', 'CD', 'NE'],
          ['CG', 'CD', 'NE', 'CZ'], ['CD', 'NE', 'CZ', 'NH1'],
          ['CD', 'NE', 'CZ', 'NH2']],
    'S': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'OG']],
    'T': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'OG1'],
          ['N', 'CA', 'CB', 'CG2']],
    'V': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG1'],
          ['N', 'CA', 'CB', 'CG2']],
    'W': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'CD1'], ['CA', 'CB', 'CG', 'CD2'],
          ['CB', 'CG', 'CD1', 'NE1'], ['CB', 'CG', 'CD2', 'CE2'],
          ['CB', 'CG', 'CD2', 'CE3'], ['CG', 'CD2', 'CE2', 'CZ2'],
          ['CG', 'CD2', 'CE3', 'CZ3'], ['CD2', 'CE2', 'CZ2', 'CH2']],
    'Y': [['CA:-1', 'N', 'CA', 'CB'], ['N', 'CA', 'CB', 'CG'],
          ['CA', 'CB', 'CG', 'CD1'], ['CA', 'CB', 'CG', 'CD2'],
          ['CB', 'CG', 'CD1', 'CE1'], ['CB', 'CG', 'CD2', 'CE2'],
          ['CG', 'CD1', 'CE1', 'CZ'], ['CD1', 'CE1', 'CZ', 'OH']]
}


class AminoBuilder:
    @staticmethod
    def supported_amino_acids():
        return [get_amino(aa) for aa in __side_chain_dihedral__.keys()]

    @staticmethod
    def parse_position_key(key):
        res = key.split(':')
        assert len(res) <= 2
        if len(res) == 1:
            res = res + [0]
        return res[0], int(res[1])

    @staticmethod
    def sidechain_dihedrals(amino):
        if isinstance(amino, str):
            amino = get_amino(amino.upper())
        assert isinstance(amino, AminoAcid)
        diheds = __side_chain_dihedral__[amino.name(one_letter_code=True)]
        return [[AminoBuilder.parse_position_key(elem) for elem in dihed] for dihed in diheds]

    @staticmethod
    def all_dihedrals(amino, model_desc = False):
        if isinstance(amino, str):
            amino = get_amino(amino.upper())
        assert isinstance(amino, AminoAcid)
        diheds = __backbone_dihedral__ + __side_chain_dihedral__[amino.name(one_letter_code=True)]
        all_dihed = [[AminoBuilder.parse_position_key(elem) for elem in dihed] for dihed in diheds]
        if model_desc:
            model_desc = __backbone_model_desc__ + __side_chain_model_desc__[amino.name(one_letter_code=True)]
            assert len(all_dihed) == len(model_desc)
            return all_dihed, model_desc
        return all_dihed

    @staticmethod
    def build_dependency_order(amino):
        if isinstance(amino, str):
            amino = get_amino(amino.upper())
        assert isinstance(amino, AminoAcid)
        return ["CA"] + [diheds[-1][0] for diheds in AminoBuilder.all_dihedrals(amino)]

    @staticmethod
    def get_dihedral(amino, atm, model_desc=False):
        if isinstance(amino, str):
            amino = get_amino(amino.upper())
        assert isinstance(amino, AminoAcid) and (atm in amino.atom_names())
        if atm == "CA":
            return None
        if model_desc:
            diheds, desc = AminoBuilder.all_dihedrals(amino, model_desc=True)
            for i, dihed in enumerate(diheds):
                if dihed[-1][0] == atm:
                    return dihed, desc[i]
        else:
            diheds = AminoBuilder.all_dihedrals(amino, model_desc=False)
            for dihed in diheds:
                if dihed[-1][0] == atm:
                    return dihed

    @staticmethod
    def model_dihedrals(amino):
        if isinstance(amino, str):
            amino = get_amino(amino.upper())
        assert isinstance(amino, AminoAcid)
        diheds, desc = AminoBuilder.all_dihedrals(amino, model_desc=True)
        return [diheds[i] for i, d in enumerate(desc) if d[-1] is None]

    @staticmethod
    def fixed_dihedrals(amino):
        if isinstance(amino, str):
            amino = get_amino(amino.upper())
        assert isinstance(amino, AminoAcid)
        dihed, desc = AminoBuilder.all_dihedrals(amino, model_desc=True)
        return [dihed[i] for i, d in enumerate(desc) if (d[0] is False) and (d[-1] is not None)]

    @staticmethod
    def coupled_dihedrals(amino):
        if isinstance(amino, str):
            amino = get_amino(amino.upper())
        assert isinstance(amino, AminoAcid)
        dihed, desc = AminoBuilder.all_dihedrals(amino, model_desc=True)
        return [(dihed[i + d[1]], dihed[i]) for i, d in enumerate(desc) if (d[0] is True) and (d[-1] is not None)]


class ProteinReconstruction:
    def __init__(self, ca_trace):
        assert isinstance(ca_trace, CaTrace)
        self.__ca_trace = deepcopy(ca_trace)
        self.__protein = {}
        self.__dihedral_full_list = None
        self.__dihedral_desc_full_list = None
        self.__dihedral_fix_list = None
        self.__index = -1
        self.__build_recons_unit()

    def __build_recons_unit(self):
        all_residues = self.__ca_trace.residue_ids
        for i, r in enumerate(all_residues):
            aa = self.__ca_trace.get_amino(r)
            if i < 2:
                self.__protein[r] = {'amino': aa,
                                     'atoms': {
                                         "CA":
                                             Coordinate3d(*self.__ca_trace.xyz(r))
                                     }
                                     }
                continue
            self.__protein[r] = {'amino': aa, 'atoms': {}}
            for atm in aa.atom_names():
                if atm != "CA":
                    self.__protein[r]['atoms'][atm] = None
                else:
                    self.__protein[r]['atoms'][atm] = Coordinate3d(*self.__ca_trace.xyz(r))

    def dihedral_build_list(self, model_desc=False, fix_filter=True):
        if self.__dihedral_full_list is None:
            all_residues = self.__ca_trace.residue_ids
            self.__dihedral_full_list = []
            self.__dihedral_desc_full_list = []
            self.__dihedral_fix_list = []
            counter = 0
            for r in all_residues[2:-1]:
                diheds, desc = AminoBuilder.all_dihedrals(self.__ca_trace.get_amino(r), model_desc=True)
                for i, dihed in enumerate(diheds):
                    atm = dihed[-1][0]
                    if self.__protein[r]['atoms'][atm] is None:
                        idx = all_residues.index(r)
                        dlist = [(pair[0], all_residues[idx + pair[1]]) for pair in dihed]
                        fix = (desc[i][-1] is None)
                        self.__dihedral_full_list.append(dlist)
                        self.__dihedral_desc_full_list.append(desc[i])
                        if fix:
                            self.__dihedral_fix_list.append(counter)
                        counter += 1
            self.__dihedral_fix_list = np.array(self.__dihedral_fix_list)
            self.__index = 0

        if self.__index == len(self.__dihedral_full_list):
            if model_desc:
                return [], []
            else:
                return []

        if not fix_filter:
            if model_desc:
                return self.__dihedral_full_list[self.__index:], self.__dihedral_desc_full_list[self.__index:]
            else:
                return self.__dihedral_full_list[self.__index:]
        else:
            res = np.where(self.__dihedral_fix_list >= self.__index)[0]
            if len(res) > 0:
                idx = res[0]
                if model_desc:
                    return [self.__dihedral_full_list[i] for i in self.__dihedral_fix_list[idx:]], \
                           [self.__dihedral_desc_full_list[i] for i in self.__dihedral_fix_list[idx:]]
                else:
                    return [self.__dihedral_full_list[i] for i in self.__dihedral_fix_list[idx:]]
            else:
                if model_desc:
                    return [], []
                else:
                    return []

    def is_complete(self):
        diheds = self.dihedral_build_list(model_desc=False)
        return len(diheds) == 0

    def curr_fix_dihedral(self):
        if not self.is_complete():
            return self.__dihedral_full_list[self.__index]
        return None

    def curr_fix_residue(self):
        if not self.is_complete():
            return self.__dihedral_full_list[self.__index][-1][1]
        return None

    def curr_fix_atom(self):
        if not self.is_complete():
            return self.__dihedral_full_list[self.__index][-1][0]
        return None

    def fix(self, d):
        diheds, desc = self.dihedral_build_list(model_desc=True, fix_filter=False)
        if len(diheds) > 0:
            dihed = diheds[0]
            dist, ang = desc[0][2], desc[0][3]
            coords = []
            for i in range(3):
                atom, residue = dihed[i][0], dihed[i][1]
                c = self.__protein[residue]["atoms"][atom]
                assert c is not None
                coords.append(c)
            atom, residue = dihed[-1][0], dihed[-1][1]
            c = reconstruct_coordinate(coords[0], coords[1], coords[2], dist, ang, d)
            self.__protein[residue]["atoms"][atom] = c
            self.__index += 1
            idx = 1
            while (idx < len(desc)) and (desc[idx][-1] is not None):
                d = d + desc[idx][-1] if desc[idx][0] is True else desc[idx][-1]
                dihed = diheds[idx]
                dist, ang = desc[idx][2], desc[idx][3]
                coords = []
                for i in range(3):
                    atom, residue = dihed[i][0], dihed[i][1]
                    c = self.__protein[residue]["atoms"][atom]
                    assert c is not None
                    coords.append(c)
                atom, residue = dihed[-1][0], dihed[-1][1]
                c = reconstruct_coordinate(coords[0], coords[1], coords[2], dist, ang, d)
                self.__protein[residue]["atoms"][atom] = c
                self.__index += 1
                idx += 1

    def get_pdb(self):
        entity, counter = [], 1
        for r in self.__ca_trace.residue_ids[2:]:
            aa = self.__ca_trace.get_amino(r)
            atoms = aa.atom_names()
            full = all([self.__protein[r]["atoms"][atm] is not None for atm in atoms])
            if full:
                for atm in atoms:
                    c = self.__protein[r]["atoms"][atm]
                    entity.append({'resid': r,
                                   'resname': aa.name(one_letter_code=False),
                                   'atomname': atm,
                                   'atomid': counter,
                                   'x': c.x,
                                   'y': c.y,
                                   'z': c.z })
                    counter = counter + 1
            else:
                break
        if len(entity) > 0:
            return PDBStructure(name=self.__ca_trace.name,
                                entry=entity,
                                chain_id=self.__ca_trace.chain)
        return None


def reconstruction_units(pdb):
    assert isinstance(pdb, PDBStructure)
    all_residues = pdb.residue_ids
    assert len(all_residues) > 2
    recons = {}
    for r in all_residues[2:]:
        amino = pdb.residue_name(r, one_letter_code=False)
        dihed_lists = AminoBuilder.all_dihedrals(amino)
        rIdx = all_residues.index(r)
        recons[r] = {}
        for i, dihed in enumerate(dihed_lists):
            assert len(dihed) == 4
            aa1, offset1 = dihed[0][0], dihed[0][1]
            aa2, offset2 = dihed[1][0], dihed[1][1]
            aa3, offset3 = dihed[2][0], dihed[2][1]
            aa4, offset4 = dihed[3][0], dihed[3][1]

            c1 = Coordinate3d(*pdb.xyz(all_residues[rIdx + offset1], aa1))
            c2 = Coordinate3d(*pdb.xyz(all_residues[rIdx + offset2], aa2))
            c3 = Coordinate3d(*pdb.xyz(all_residues[rIdx + offset3], aa3))
            c4 = Coordinate3d(*pdb.xyz(all_residues[rIdx + offset4], aa4))
            recons[r][aa4] = (amino,
                              distance(c3, c4),
                              angle(c2, c3, c4),
                              dihedral(c1, c2, c3, c4))
    return recons


