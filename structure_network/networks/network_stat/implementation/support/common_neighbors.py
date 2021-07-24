import logging
import networkx as nx
from structure_network.networks.network_stat.interfaces import SupportFunction

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ["CommonNeighborSupport"]


class CommonNeighborSupport(SupportFunction):
    def __init__(self, g, **kwargs):
        assert isinstance(g, nx.Graph)
        self._g = g
        self._name = kwargs.get('name', 'common_nbrs')
        self._logger = logging.getLogger(self.__class__.__name__)

    @property
    def name(self):
        return self._name

    @property
    def attribute(self):
        return None

    def __call__(self, x, y):
        self._logger.debug("Queried for node pair {} and {}".format(x, y))
        assert self._g.has_node(x) and self._g.has_node(y)
        return len(set(nx.neighbors(self._g, x)).intersection(nx.neighbors(self._g, y)))
