import logging
import networkx as nx
from structure_network.networks.network_stat.interfaces import SupportFunction

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['ShortestDistanceSupport']


class ShortestDistanceSupport(SupportFunction):
    def __init__(self, g, **kwargs):
        assert isinstance(g, nx.Graph)
        self._g = g
        self._logger = logging.getLogger(self.__class__.__name__)
        self._name = kwargs.get('name', 'st')
        self._attr = kwargs.get('attribute', 'weight')

    @property
    def name(self):
        return self._name

    @property
    def attribute(self):
        return self._attr

    def __call__(self, x, y):
        self._logger.debug("Queried for node pair {} and {}".format(x, y))
        assert self._g.has_node(x) and self._g.has_node(y)
        return nx.shortest_path_length(self._g, x, y, weight=self._attr)
