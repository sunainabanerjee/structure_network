import threading
import numpy as np
import networkx as nx
from structure_networks.networks.network_stat.interfaces import NetworkGroupSelector

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NodeNeighborSelector']


class NodeNeighborSelector(NetworkGroupSelector):
    def __init__(self, g, **kwargs):
        assert isinstance(g, nx.Graph)
        self._g = g
        self._depth = 0
        self._repeat_after = kwargs.get('replace_after', self._g.order())
        self._visited_node = {}
        self._curr_node = kwargs.get('start_node', np.random.choice(self._g.nodes))
        if self._curr_node not in self._g.nodes:
            raise RuntimeError("Error: invalid node {} "
                               "specified as start node".format(self._curr_node))
        self._lock = threading.Lock()

    @property
    def n_pair(self):
        return 1

    def current_node(self, node=None):
        if node is None:
            return self._curr_node
        else:
            self._lock.acquire()
            try:
                assert node in self._g.nodes
                self._curr_node = node
            finally:
                self._lock.release()
        return self

    def reset(self):
        self._lock.acquire()
        try:
            self._depth = 0
            self._visited_node.clear()
        finally:
            self._lock.release()
        return self

    def next(self):
        self._lock.acquire()
        try:
            nbr_list = list(nx.neighbors(self._g, self._curr_node))
            nbr_list = [ n for n in nbr_list
                         if (n not in self._visited_node) or
                         ((self._depth - self._visited_node[n]) > self._repeat_after)]
            if len(nbr_list) == 0:
                raise StopIteration
            node_visit = self._curr_node
            self._visited_node[node_visit] = self._depth
            self._curr_node = np.random.choice(nbr_list)
            self._depth += 1
        finally:
            self._lock.release()
        return node_visit
