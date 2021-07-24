import numpy as np
import networkx as nx
from ..graph_algorithms import diameter
from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist
from structure_network.structure import get_pair_potential
from structure_network.structure import CaTrace
from structure_network.structure import PDBStructure
from structure_network.structure import pdb_to_catrace
from structure_network.structure import PairPotential
from .structure_network import ProteinStructureNetwork
from .contact_categories import ContactTypes, ThresholdTypes

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__version__ = "1.0"
__all__ = ["CANetwork"]


class CANetwork(ProteinStructureNetwork):
    def __init__(self,
                 structure,
                 **kwargs):
        if isinstance(structure, PDBStructure):
            structure = pdb_to_catrace(structure)

        distance_cutoff = kwargs.get('distance_threshold', 6.5)
        energy_cutoff = kwargs.get('energy_threshold', 0.0)
        threshold_type = kwargs.get('threshold_type', ThresholdTypes.distance_energy)
        contact_type = kwargs.get('contact_type', ContactTypes.all_contacts)
        potential = kwargs.get('potential', PairPotential.charmm)

        assert isinstance(structure, CaTrace)
        assert distance_cutoff > 0
        self.__contact_type = contact_type
        self.__contact_potential = potential
        self.__threshold_type = threshold_type
        self.__distance_threshold = float(distance_cutoff)
        self.__energy_threshold = float(energy_cutoff)
        self.__structure = structure
        self.__g = None
        self.__build()
        self.__max_wt = np.max(list(nx.get_edge_attributes(self.__g, 'weight').values()))
        self.__min_wt = np.min(list(nx.get_edge_attributes(self.__g, 'weight').values()))

    def __len__(self):
        return self.order

    @property
    def order(self):
        return self.__g.order()

    @property
    def size(self):
        return self.__g.size()

    @property
    def max_edge_weight(self):
        return self.__max_wt

    @property
    def min_edge_weight(self):
        return self.__min_wt

    @property
    def threshold(self):
        return {'distance': self.__distance_threshold,
                'energy': self.__energy_threshold}

    @property
    def structure(self):
        return self.__structure

    def set_distance_threshold(self, d):
        assert d > 0
        self.__distance_threshold = float(d)
        self.__build()
        return self

    def set_energy_threshold(self, d):
        assert d >= 0
        self.__energy_threshold = float(d)
        self.__build()
        return self

    @property
    def contact_type(self):
        return self.__contact_type

    def __build(self):
        self.__g = nx.Graph()
        residue_ids = self.__structure.residue_ids
        residue_names = ["%s%d" % (self.__structure.get_amino(r).name(), r) for r in residue_ids]
        xyz = np.array([list(self.__structure.xyz(r)) for r in residue_ids])
        dist_mat = cdist(xyz, xyz)
        from_x, to_x = np.where(dist_mat <= self.__distance_threshold)
        for i, r in enumerate(residue_names):
            self.__g.add_node(r,
                              name=self.__structure.get_amino(residue_ids[i]).name(),
                              resid=residue_ids[i],
                              xyz=self.__structure.xyz(residue_ids[i]))
        for i, u in enumerate(from_x):
            v = to_x[i]
            if (self.__contact_type == ContactTypes.backbone_only) and abs(u - v) > 1:
                continue
            if (self.__contact_type == ContactTypes.non_bonded_only) and abs(u - v) == 1:
                continue

            u_name = residue_names[u]
            v_name = residue_names[v]
            wt = get_pair_potential(u_name[:3], v_name[:3],
                                    distance=dist_mat[u, v],
                                    pot_type=self.__contact_potential)

            add_edge = (u_name != v_name)
            if self.__threshold_type in (ThresholdTypes.distance_energy, ThresholdTypes.energy):
                add_edge = add_edge and (wt > self.__energy_threshold)

            if add_edge:
                self.__g.add_edge(u_name, v_name, weight=wt)

    @property
    def residue_ids(self):
        return self.__structure.residue_ids

    @property
    def g(self):
        return self.__g

    @property
    def residue_names(self):
        return ["%s%d" % (self.__structure.get_amino(r).name(), r)
                for r in self.__structure.residue_ids]

    def xyz(self, r=None):
        if r is None:
            return np.array([tuple(self.__g.nodes[r]['xyz']) for r in self.residue_names])
        residue_ids = self.__structure.residue_ids
        residue_names = self.residue_names
        if r in residue_ids:
            r = residue_names[residue_ids.index(r)]
        if r not in residue_names:
            raise RuntimeError("Error: unknown "
                               "residue reference {}".format(r))
        return self.__g.nodes[r]['xyz']

    @property
    def centroid(self):
        xyz = self.xyz()
        wts = np.array([self.__structure.get_amino(r).molecular_weight
                        for r in self.residue_ids])
        return tuple(np.sum((xyz.T * wts).T, axis=0)/np.sum(wts))

    @property
    def radius_of_gyration(self):
        xc, yc, zc = self.centroid
        cntr = np.array([xc, yc, zc])
        wts = np.array([self.__structure.get_amino(r).molecular_weight
                        for r in self.residue_ids])
        xyz = self.xyz()
        return np.sqrt(np.sum(np.sum(np.square((xyz - cntr)),
                                     axis=1) * wts)/np.sum(wts))

    def get_amino(self, r):
        r_names = self.residue_names
        if r in r_names:
            r = self.residue_ids[r_names.index(r)]
        return self.__structure.get_amino(r)

    @property
    def threshold_type(self):
        return self.__threshold_type

    @property
    def structure_diameter(self):
        return np.max(pdist(self.xyz()))

    @property
    def network_diameter(self):
        return diameter(self.g)

    @property
    def params(self):
        return {'distance_threshold': self.__distance_threshold,
                'energy_threshold': self.__energy_threshold,
                'threshold_type': self.__threshold_type,
                'contact_type': self.__contact_type,
                'potential': self.__contact_potential}

    def copy(self):
        return CANetwork(structure=self.__structure,
                         **self.params)


