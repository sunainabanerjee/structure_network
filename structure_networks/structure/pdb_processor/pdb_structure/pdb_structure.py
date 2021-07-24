"""
This program file defines the basic class object
for storing pdb 3D structure. The basic class structure
only allows to read and access the coordinates and related
information in a systematic way.
"""

import logging
import numpy as np
from copy import deepcopy
from structure_networks.structure.amino_acids import get_amino
from structure_networks.structure.force_field import FFManager

__version__ = "1.0"
__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['PDBStructure']


class PDBStructure:
    def __init__(self, name, entry, chain_id='A'):
        self.__name = name
        self.__chain_id = chain_id
        self.__structure = dict()
        self.__residues = dict()
        self.__bfactor = dict()

        logger = logging.getLogger(name='pdb_processor.PDBStructure')
        assert isinstance(entry, list)
        ff = FFManager()
        keys = ['resid', 'resname', 'atomname', 'atomid', 'x', 'y', 'z']
        for item in entry:
            assert isinstance(item, dict)
            for k in keys:
                assert k in item
            residue_id, residue_name, atom_name, atom_id = item['resid'], item['resname'], item['atomname'], item[
                'atomid']
            if residue_id not in self.__structure:
                self.__structure[residue_id] = dict()
            if not ff.is_valid_atom(residue_name, atom_name):
                logger.warning('Invalid atom pair (%s, %s) ignoring!!' % (residue_name, atom_name))
                continue
            charge = ff.get_charge(residue_name, atom_name) if 'q' not in item else item['q']
            self.__structure[residue_id][atom_name] = {'x': item['x'],
                                                       'y': item['y'],
                                                       'z': item['z'],
                                                       'id': atom_id,
                                                       'q': charge}
            self.__residues[residue_id] = get_amino(residue_name)
            if residue_id not in self.__bfactor:
                self.__bfactor[residue_id] = {}
            self.__bfactor[residue_id][atom_name] = item['bfactor'] if 'bfactor' in item else 100

    @property
    def name(self):
        return self.__name

    @property
    def chain(self):
        return self.__chain_id

    @property
    def size(self):
        return len(self.__structure)

    def get_amino(self, residue_id):
        return self.__residues[residue_id]

    @property
    def is_broken_resid_sequence(self):
        return len(np.where(np.diff(self.residue_ids) != 1)[0]) > 0

    @property
    def is_improper_calphas(self):
        min_alpha, max_alpha = 2.9, 4.0
        ca_xyz = np.array([[*self.xyz(r, 'CA')] for r in self.residue_ids], dtype=np.float)
        ca_dist = np.sqrt(np.sum(np.square(ca_xyz[:-1] - ca_xyz[1:]), axis=1))
        return len(np.where((ca_dist < min_alpha) | (ca_dist > max_alpha))[0]) > 0

    @property
    def has_missing_atoms(self):
        complete = True
        for res_id in self.residue_ids:
            a_list = self.atom_names(residue_id=res_id)
            a_full = self.get_amino(res_id).atom_names
            complete = complete and all([atm in a_list for atm in a_full])
            if not complete:
                break
        return not complete

    def missing_atom_list(self):
        residue_list = {}
        for res_id in self.residue_ids:
            a_list = self.atom_names(residue_id=res_id)
            a_full = self.get_amino(res_id).atom_names
            missing = [atom for atom in a_full if atom not in a_list]
            if len(missing) > 0:
                residue_list[res_id] = missing
        return residue_list

    @property
    def sequence(self):
        s = list()
        for k in sorted([int(r) for r in self.__residues.keys()]):
            s.append(self.__residues[k])
        return s

    @property
    def residue_ids(self):
        return sorted([int(k) for k in self.__structure.keys()])

    def find_residue(self, aa_type):
        aa = get_amino(aa_type)
        return [rid for rid in self.__residues.keys() if self.__residues[rid] == aa]

    def atom_list(self, residue_id):
        if residue_id in self.__structure:
            return [self.__structure[residue_id][atm] for atm in self.atom_names(residue_id)]
        return []

    def atom_names(self, residue_id):
        if residue_id in self.__structure:
            return [atm for atm in
                    sorted(self.__structure[residue_id].keys(), key=lambda x: self.__structure[residue_id][x]['id'])]
        return []

    def charge(self, residue_id, atm):
        if (residue_id in self.__structure) and (atm in self.__structure[residue_id]):
            return self.__structure[residue_id][atm]['q']

    def key(self, residue_id):
        if residue_id in self.__residues:
            return '%s%d' % (self.__residues[residue_id], residue_id)
        return ''

    def write(self, fh):
        if not hasattr(fh, 'write'):
            raise Exception('Invalid object type [expects file object]')
        res_ids = self.residue_ids
        for i, r in enumerate(res_ids):
            for j, aa in enumerate(self.atom_names(r)):
                line = "ATOM  %5d  %-3s %3s %1s%4d    " \
                       "%8.3f%8.3f%8.3f%6.2f%6.2f          %-2s" % (self.__structure[r][aa]['id'],
                                                                    aa,
                                                                    self.__residues[r],
                                                                    self.__chain_id,
                                                                    int(r),
                                                                    float(self.__structure[r][aa]['x']),
                                                                    float(self.__structure[r][aa]['y']),
                                                                    float(self.__structure[r][aa]['z']),
                                                                    1.0,
                                                                    float(self.__bfactor[r][aa]),
                                                                    aa[0])
                fh.write("%s\n" % line)

    def __len__(self):
        return len(self.__structure)

    def __str__(self):
        s = ''
        for aa in self.sequence:
            s = s + aa.name(one_letter_code=True)
        return s

    def residue_name(self, residue_id, one_letter_code=False):
        if residue_id in self.__residues:
            return self.__residues[residue_id].name(one_letter_code=one_letter_code)
        return None

    def xyz(self, residue_id, atom_name):
        if (residue_id in self.__structure) and (atom_name in self.__structure[residue_id]):
            atom_details = self.__structure[residue_id][atom_name]
            return atom_details['x'], atom_details['y'], atom_details['z']
        return None, None, None

    def b_factor(self, residue_id, atom_name, value=None):
        if (residue_id in self.__bfactor) and (atom_name in self.__bfactor[residue_id]):
            if value is None:
                return self.__bfactor[residue_id][atom_name]
            else:
                self.__bfactor[residue_id][atom_name] = float(value)
                return self
        return None

    @property
    def coordinate_min(self):
        resId = self.residue_ids[0]
        atom = list(self.__structure[resId].keys())[0]
        x = self.__structure[resId][atom]['x']
        y = self.__structure[resId][atom]['y']
        z = self.__structure[resId][atom]['z']
        for r in self.__structure.keys():
            for a in self.__structure[r].keys():
                x = self.__structure[r][a]['x'] if x > self.__structure[r][a]['x'] else x
                y = self.__structure[r][a]['y'] if y > self.__structure[r][a]['y'] else y
                z = self.__structure[r][a]['z'] if z > self.__structure[r][a]['z'] else z
        return x, y, z

    @property
    def coordinate_max(self):
        resId = self.residue_ids[0]
        atom = list(self.__structure[resId].keys())[0]
        x = self.__structure[resId][atom]['x']
        y = self.__structure[resId][atom]['y']
        z = self.__structure[resId][atom]['z']
        for r in self.__structure.keys():
            for a in self.__structure[r].keys():
                x = self.__structure[r][a]['x'] if x < self.__structure[r][a]['x'] else x
                y = self.__structure[r][a]['y'] if y < self.__structure[r][a]['y'] else y
                z = self.__structure[r][a]['z'] if z < self.__structure[r][a]['z'] else z
        return x, y, z

    def copy(self):
        return deepcopy(self)



