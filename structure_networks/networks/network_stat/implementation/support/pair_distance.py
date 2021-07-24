import logging
import numpy as np
import networkx as nx
from structure_networks.networks.network_stat.interfaces import SupportFunction

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['DistanceSupport']


class DistanceSupport(SupportFunction):
    def __init__(self, g, **kwargs):
        assert isinstance(g, nx.Graph)
        self._g = g
        self._name = kwargs.get('name', 'distance')
        self._attr = kwargs.get('attribute', 'xyz')
        self._logger = logging.getLogger(self.__class__.__name__)
        assert len(nx.get_node_attributes(self._g, self._attr)) > 0, \
            'Error: no node attribute found {}'.format(self._attr)

    @property
    def name(self):
        return self._name

    @property
    def attribute(self):
        return self._attr

    def __call__(self, x, y):
        self._logger.debug("Queried for node pair {} and {}".format(x, y))
        assert self._g.has_node(x) and self._g.has_node(y)
        xc = np.array(self._g.nodes[x][self.attribute])
        yc = np.array(self._g.nodes[y][self.attribute])
        return np.sqrt(np.sum(np.square(xc - yc)))

