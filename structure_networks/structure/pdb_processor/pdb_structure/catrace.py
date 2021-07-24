import numpy as np
from copy import deepcopy
from structure_networks.structure.amino_acids import get_amino

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['CaTrace']


class CaTrace:
    """
    This is a basic pdb structure class. This class reads only
    CA trace atom trace from a given pdb file. This class relies
    on the ATOM tags of the record. In case of modified residues
    with HETATM tag, this class may results in missing residues
    in the trace.
    """

    def __init__(self, name, entry, chainId='A'):
        self.__name = name
        self.__bfactor = dict()
        self.__chain_id = chainId
        self.__structure = dict()
        self.__residues = dict()

        assert isinstance(entry, list)
        keys = ['resid', 'resname', 'x', 'y', 'z']
        for item in entry:
            assert isinstance(item, dict)
            for k in keys:
                assert k in item
            self.__structure[item['resid']] = {'x': item['x'],
                                               'y': item['y'],
                                               'z': item['z']}
            self.__residues[item['resid']] = get_amino(item['resname'])
            self.__bfactor[item['resid']] = item['bfactor'] if 'bfactor' in item else 100

    @property
    def name(self):
        return self.__name

    @property
    def chain(self):
        return self.__chain_id

    @property
    def size(self):
        return len(self.__structure)

    @property
    def sequence(self):
        s = list()
        for k in sorted([int(r) for r in self.__residues.keys()]):
            s.append(self.__residues[k])
        return s

    @property
    def residue_ids(self):
        return sorted([int(k) for k in self.__structure.keys()])

    def key(self, residue_id):
        if residue_id in self.__residues:
            return '%s%d' % (self.__residues[residue_id], residue_id)
        return ''

    def xyz(self, resid):
        residue_id = int(resid)
        if residue_id not in self.__structure:
            raise KeyError('Invalid residue id: %d' % residue_id)
        return self.__structure[residue_id]['x'], self.__structure[residue_id]['y'], self.__structure[residue_id]['z']

    def b_factor(self, resid, value=None):
        residue_id = int(resid)
        if residue_id not in self.__bfactor:
            raise KeyError('Invalid residue id: %d' % residue_id)
        if value is None:
            return self.__bfactor[residue_id]
        else:
            self.__bfactor[residue_id] = value
            return value

    def get_amino(self, res_id):
        residue_id = int(res_id)
        if residue_id not in self.__residues:
            raise KeyError('Invalid residue id: {}'.format(residue_id))
        return self.__residues[residue_id]

    def set_amino(self, res_id, aa_type='A'):
        residue_id = int(res_id)
        if residue_id not in self.__residues:
            return self
        self.__residues[residue_id] = get_amino(aa_type)
        return self

    def write(self, fh):
        if not hasattr(fh, 'write'):
            raise Exception('Invalid object type [expects file object]')
        res_ids = self.residue_ids
        for i, r in enumerate(res_ids):
            line = "ATOM  %5d  %-3s %3s %1s%4d    " \
                   "%8.3f%8.3f%8.3f%6.2f%6.2f          %-2s" % (i,
                                                                'CA',
                                                                self.__residues[r],
                                                                self.__chain_id,
                                                                int(r),
                                                                float(self.__structure[r]['x']),
                                                                float(self.__structure[r]['y']),
                                                                float(self.__structure[r]['z']),
                                                                1.0,
                                                                float(self.__bfactor[r]),
                                                                'C')
            fh.write("%s\n" % line)

    def __len__(self):
        return len(self.__structure)

    def __str__(self):
        s = ''
        for aa in self.sequence:
            s = s + aa.name(one_letter_code=True)
        return s

    @property
    def coordinate_min(self):
        resId = self.residue_ids[0]
        x, y, z = self.__structure[resId]['x'], \
                  self.__structure[resId]['y'], \
                  self.__structure[resId]['z']
        for r in self.__structure.keys():
            x = self.__structure[r]['x'] if x > self.__structure[r]['x'] else x
            y = self.__structure[r]['y'] if y > self.__structure[r]['y'] else y
            z = self.__structure[r]['z'] if z > self.__structure[r]['z'] else z
        return x, y, z

    @property
    def coordinate_max(self):
        resId = self.residue_ids[0]
        x, y, z = self.__structure[resId]['x'], \
                  self.__structure[resId]['y'], \
                  self.__structure[resId]['z']
        for r in self.__structure.keys():
            x = self.__structure[r]['x'] if x < self.__structure[r]['x'] else x
            y = self.__structure[r]['y'] if y < self.__structure[r]['y'] else y
            z = self.__structure[r]['z'] if z < self.__structure[r]['z'] else z
        return x, y, z

    def coordinate_array(self, residue_list=None):
        if residue_list is None:
            residue_list = self.residue_ids
        assert isinstance(residue_list, list) and len(residue_list) > 0
        return np.array([[*self.xyz(r)] for r in residue_list], dtype=np.float)

    def copy(self):
        return deepcopy(self)
