import logging
import networkx as nx
from structure_networks.structure import PairPotential
from structure_networks.structure import get_pair_potential
from structure_networks.networks.network_stat.interfaces import PairLookup

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['DistanceDependentPotentialLookup']


class DistanceDependentPotentialLookup(PairLookup):
    def __init__(self, **kwargs):
        self._g = kwargs.get('g', None)
        self._potential = kwargs.get('potential', PairPotential.charmm)
        self._lookup_name = kwargs.get('node_attr', 'name')
        self._distance_attr = kwargs.get('edge_attr', 'distance')
        self._logger = logging.getLogger(self.__class__.__name__)
        self._check_ready = False
        self._check_status = False

    @property
    def is_ready(self):
        if not self._check_status:
            self._check_ready = isinstance(self._g, nx.Graph) and \
                                (len(nx.get_node_attributes(self._g, self._lookup_name)) == self._g.order()) and \
                                (len(nx.get_edge_attributes(self._g, self._distance_attr)) == self._g.size())
            self._check_status = True
        return self._check_ready

    def set_network_context(self, g):
        if isinstance(g, nx.Graph):
            self._g = g
            self._check_status = False
        return self

    @property
    def is_random(self):
        return False

    @property
    def is_stateful(self):
        return False

    def reset(self):
        return self

    def __call__(self, x, y, **kwargs):
        if not self.is_ready:
            raise RuntimeError("Error: graph context is not setup!")
        assert (x in self._g.nodes) and (y in self._g.nodes)
        d = kwargs.get('distance', self._g[x][y][self._distance_attr])
        return get_pair_potential(amino1=self._g.nodes[x][self._lookup_name],
                                  amino2=self._g.nodes[y][self._lookup_name],
                                  distance=d,
                                  pot_type=self._potential)

