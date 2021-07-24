import os
import logging
import numpy as np
from .charmm import *

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__all__ = ["FFManager"]


class FFManager:
    class ProteinFF:
        def __init__(self):
            cwd = os.path.dirname(os.path.abspath(__file__))
            rtf_file = os.path.join(cwd, 'data', 'prot_top_charmm.rtf')
            prm_file = os.path.join(cwd, 'data', 'prot_nonbonded_charmm.prm')
            atoms, bonds, donors, acceptors = parse_amino_charmff(rtf_file)
            self.__logger = logging.getLogger(name="FFManager.ProteinFF")
            self.__atoms = atoms
            self.__bonds = bonds
            self.__donors = donors
            self.__acceptors = acceptors
            self.__ccelec = 331.843
            self.__name_map = {}
            for n in self.__atoms.keys():
                if n == 'HSE':
                    self.__name_map['HIS'] = n
                self.__name_map[n] = n
            self.__parameters = parse_nb_charmff(prm_file)

        def get_residue_names(self):
            return list(self.__name_map.keys())

        def is_valid_residue(self, residue_name):
            return residue_name in self.__name_map

        def __get_mapped_name(self, residue_name):
            if residue_name in self.__name_map:
                return self.__name_map[residue_name]
            return None

        def get_atom_names(self, residue_name, noh=False ):
            if self.is_valid_residue(residue_name):
                mapped_name = self.__get_mapped_name(residue_name)
                if noh:
                    return [aa for aa in self.__atoms[mapped_name].keys() if not self.__atoms[mapped_name][aa][0].startswith('H')]
                else:
                    return [aa for aa in self.__atoms[mapped_name].keys()]
            return []

        def is_valid_atom(self, residue_name, atom_name):
            return self.is_valid_residue(residue_name) and (atom_name in self.__atoms[self.__get_mapped_name(residue_name)])

        def get_charge(self, residue_name, atom_name):
            if self.is_valid_atom(residue_name, atom_name):
                return self.__atoms[self.__get_mapped_name(residue_name)][atom_name][1]
            return 0

        def get_atom_type(self, residue_name, atom_name):
            if self.is_valid_atom(residue_name, atom_name):
                return self.__atoms[self.__get_mapped_name(residue_name)][atom_name][0]
            return None

        def __get_epsilon(self, residue_name, atom_name, type=""):
            atype = self.get_atom_type(residue_name, atom_name)
            if (atype is not None) and (atype in self.__parameters):
                return self.__parameters[atype][2] if type == "14" else self.__parameters[atype][0]
            return 0

        def __get_radius(self, residue_name, atom_name, type=""):
            atype = self.get_atom_type(residue_name, atom_name)
            if (atype is not None) and (atype in self.__parameters):
                return self.__parameters[atype][3] if type == "14" else self.__parameters[atype][1]
            return 0

        def get_donor_atoms(self, residue_name):
            if self.is_valid_residue(residue_name):
                mapped_name = self.__get_mapped_name(residue_name)
                if mapped_name in self.__donors:
                    return self.__donors[mapped_name]
            return []

        def get_acceptor_atoms(self, residue_name):
            if self.is_valid_residue(residue_name):
                mapped_name = self.__get_mapped_name(residue_name)
                if mapped_name in self.__acceptors:
                    return self.__acceptors[mapped_name]
            return []

        def is_donor(self, residue_name, atom_name):
            if self.is_valid_atom(residue_name, atom_name):
                return atom_name in self.get_donor_atoms(residue_name)
            return False

        def is_acceptor(self, residue_name, atom_name):
            if self.is_valid_atom(residue_name, atom_name):
                return atom_name in self.get_acceptor_atoms(residue_name)
            return False

        def calculate_elec_energy(self,
                                  residue_name1,
                                  atom_name1,
                                  residue_name2,
                                  atom_name2,
                                  distance,
                                  epsilon=0.8):
            distance = np.clip(distance, a_min=1e-3, a_max=None)
            q1 = self.get_charge(residue_name1, atom_name1)
            q2 = self.get_charge(residue_name2, atom_name2)
            e = (q1 * q2 * self.__ccelec)/(distance * epsilon)
            self.__logger.debug("Electrostatic energy: (%s:%s <=>  %s:%s) is [%f]" %(residue_name1,
                                                                                     atom_name1,
                                                                                     residue_name2,
                                                                                     atom_name2,
                                                                                     e))
            return e

        def calculate_nb_energy(self,
                                residue_name1,
                                atom_name1,
                                residue_name2,
                                atom_name2,
                                distance):
            assert isinstance(distance, np.float) and distance >= 0.
            d = np.clip(distance, a_min=1e-3, a_max=None)

            if self.is_valid_atom(residue_name1, atom_name1) and self.is_valid_atom(residue_name2, atom_name2) and d < 12:
                rmin = self.__get_radius(residue_name1, atom_name1) + self.__get_radius(residue_name2, atom_name2)
                eps = np.sqrt(self.__get_epsilon(residue_name1, atom_name1) * self.__get_epsilon(residue_name2, atom_name2))
                e = eps * ((rmin/d)**12 - 2*(rmin/d)**6)
                self.__logger.debug("Nonbonded energy between (%s:%s <=> %s:%s) [%f]" % (residue_name1,
                                                                                         atom_name1,
                                                                                         residue_name2,
                                                                                         atom_name2,
                                                                                         e))
                return e
            return 0

        def calculate_energy(self,
                             residue_name1,
                             atom_name1,
                             residue_name2,
                             atom_name2,
                             distance,
                             epsilon=1.,
                             elec_only=False,
                             summed=True):
            elec, nb = 0, 0
            if elec_only or summed:
                elec = self.calculate_elec_energy(residue_name1, atom_name1, residue_name2, atom_name2, distance, epsilon)
            if not elec_only or summed:
                nb = self.calculate_nb_energy(residue_name1, atom_name1, residue_name2, atom_name2, distance)
            return elec + nb

    __instance = None

    def __init__(self):
        if FFManager.__instance is None:
            FFManager.__instance = FFManager.ProteinFF()

    def is_valid_residue(self, residue_name):
        return self.__instance.is_valid_residue(residue_name)

    def get_charge(self, residue_name, atom_name):
        return self.__instance.get_charge(residue_name, atom_name)

    def is_valid_atom(self, residue_name, atom_name):
        return self.__instance.is_valid_atom(residue_name, atom_name)

    @property
    def residue_names(self):
        return FFManager.__instance.get_residue_names()

    def atom_names(self, residue_name, noh=False):
        return FFManager.__instance.get_atom_names(residue_name, noh=noh)

    def get_acceptors(self, residue_name):
        return FFManager.__instance.get_acceptor_atoms(residue_name)

    def is_acceptor(self, residue_name, atom_name):
        return FFManager.__instance.is_acceptor(residue_name, atom_name)

    def get_donors(self, residue_name):
        return FFManager.__instance.get_donor_atoms(residue_name)

    def is_donor(self, residue_name, atom_name):
        return FFManager.__instance.is_donor(residue_name, atom_name)

    def nb_energy(self, residue1, atom1, residue2, atom2, distance):
        return FFManager.__instance.calculate_nb_energy(residue1, atom1, residue2, atom2, distance)

    def elec_energy(self, residue1, atom1, residue2, atom2, distance, epsilon=1.):
        return FFManager.__instance.calculate_elec_energy(residue1, atom1, residue2, atom2, distance, epsilon)

    def energy(self, residue1, atom1, residue2, atom2, distance, epsilon, elec_only=False, summed=True):
        return FFManager.__instance.calculate_energy(residue1,
                                                     atom1,
                                                     residue2,
                                                     atom2,
                                                     distance,
                                                     epsilon,
                                                     elec_only,
                                                     summed)

