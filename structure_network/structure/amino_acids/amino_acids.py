"""
This code base reports all supported amino acid residue types
and also supports their associated property index.
"""

from .properties import *

__version__ = 1.0
__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['AminoAcid',
           'get_amino',
           'valid_amino_acids',
           'get_vdw_radius']


def valid_amino_acids(one_letter=True):
    if one_letter is True:
        return sorted(amino_acids.keys())
    else:
        return [amino_acids[a] for a in sorted(amino_acids.keys())]


class AminoAcid:
    def __init__(self, aa_type="G"):
        assert aa_type in amino_acids
        self.aa = aa_type

    def __str__(self):
        return '%s' % amino_acids[self.aa]

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, amino_acids[self.aa])

    def __format__(self, **kwargs):
        return '%s' % amino_acids[self.aa]

    def __hash__(self):
        return hash(self.name())

    def __eq__(self, other):
        if isinstance(other, AminoAcid):
            return self.aa == other.aa
        elif isinstance(other, str):
            return self.aa == other
        else:
            return False

    def __lt__(self, other):
        return (type(other) == type(self)) and (self.aa < other.aa)

    def name(self, one_letter_code=False):
        if one_letter_code is True:
            return self.aa
        else:
            return amino_acids[self.aa]

    @property
    def molecular_weight(self):
        return aa_mw[self.aa]

    @property
    def sasa_free(self):
        return aa_sasa_free[self.aa]

    @property
    def sasa_folded(self):
        return aa_sasa_folded[self.aa]

    @property
    def volume(self):
        return aa_volume[self.aa]

    @property
    def flexibility_index(self):
        return aa_flexibility[self.aa]

    @property
    def hydropathy_index(self):
        return aa_hydropathy[self.aa]

    @property
    def buried_propensity(self):
        return aa_buried_propensity[self.aa]

    @property
    def pIsoelectric(self):
        return aa_pI[self.aa]

    @property
    def atom_names(self):
        return atoms[self.aa].copy()

    @property
    def charge_atoms(self):
        return [aa for aa in self.atom_names if aa[0] in {'N', 'O', 'S'}]

    @property
    def hb_donors(self):
        return hb_donor[self.aa].copy()

    @property
    def hb_acceptors(self):
        return hb_acceptor[self.aa].copy()

    @property
    def bonds(self):
        return [(v1.split(":")[0],
                 v2.split(":")[0]) for v1, v2 in bonds[self.aa]
                if v1.split(":")[1] == v2.split(":")[1]]

    @property
    def angles(self):
        bond_map = {}
        for v1, v2 in self.bonds:
            if v1 not in bond_map:
                bond_map[v1] = []
            if v2 not in bond_map:
                bond_map[v2] = []
            bond_map[v1].append(v2)
            bond_map[v2].append(v1)

        angle_list = []

        def find_bonded_triplet(links):
            if len(links) < 3:
                assert len(links) > 0
                aa = links[-1]
                if aa in bond_map:
                    for ab in bond_map[aa]:
                        if ab not in links:
                            find_bonded_triplet(links + [ab])
            elif len(links) == 3:
                angle_list.append(links)

        for atom in self.atom_names:
            find_bonded_triplet([atom])
        return angle_list


def get_amino(amino_name):
    if isinstance(amino_name, AminoAcid):
        return amino_name
    assert isinstance(amino_name, str)
    amino_name = amino_name.strip()
    if len(amino_name) == 3:
        assert amino_name in amino_acids.values()
        for key, value in amino_acids.items():
            if amino_name == value:
                return AminoAcid(aa_type=key)
    else:
        assert len(amino_name) == 1 and amino_name in amino_acids
        return AminoAcid(aa_type=amino_name)
    raise KeyError("Invalid amino name: {}".format(amino_name))


def get_vdw_radius(atom_name):
    assert isinstance(atom_name, str)
    atom_name = atom_name.strip()
    assert atom_name[0] in vdw_radius.keys()
    return vdw_radius[atom_name[0]]
