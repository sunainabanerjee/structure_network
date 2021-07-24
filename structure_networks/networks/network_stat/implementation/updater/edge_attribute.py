import networkx as nx
from structure_networks.networks.core import to_list
from structure_networks.networks.network_stat.interfaces import UpdaterTypes
from structure_networks.networks.network_stat.interfaces import NetworkUpdater
from structure_networks.networks.network_stat.interfaces import PairLookup
from structure_networks.networks.network_stat.interfaces import SupportFunction

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NetworkEdgeAttributeUpdater']


class NetworkEdgeAttributeUpdater(NetworkUpdater):
    def __init__(self,
                 weight_lookup,
                 support_func=None):
        self.__lookup_fn = weight_lookup
        self.__support_fn = to_list(support_func)
        assert isinstance(self.__lookup_fn, PairLookup)
        assert all([isinstance(fn, SupportFunction)
                    for fn in self.__support_fn])

    @property
    def type(self):
        return UpdaterTypes.edge_attribute

    @property
    def attribute(self):
        return

    def __call__(self, g, **kwargs):
        assert isinstance(g, nx.Graph)
        nodes = [x for x in to_list(kwargs.get('nodes', None))
                 if g.has_node(x)]
        if len(nodes) == 0:
            edges = g.edges
        else:
            edges = []
            for x in nodes:
                for n in g.neighbors(x):
                    edges.append((x, n))
        edge_attr = kwargs.get('edge_attr', 'weight')
        self.__lookup_fn.set_network_context(g)
        for u, v in edges:
            args = {}
            for fn in self.__support_fn:
                args[fn.name] = fn(u, v)
            g[u][v][edge_attr] = self.__lookup_fn(u, v, **args)
        return g


