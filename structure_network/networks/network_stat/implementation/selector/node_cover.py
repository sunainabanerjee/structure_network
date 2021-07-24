import threading
import numpy as np
import networkx as nx
from structure_network.networks.network_stat.interfaces import NetworkGroupSelector

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NodeCoverSelector']


class NodeCoverSelector(NetworkGroupSelector):
    def __init__(self, g, node_group=None, k=None):
        assert isinstance(g, nx.Graph)
        self._g = g
        if node_group is None:
            node_group = list(self._g.nodes)
        self._node_list = [n for n in node_group if n in self._g.nodes]
        self._counter = 0
        self._visited = set()
        self._k = len(node_group) if k is None else min(len(node_group), k)
        self._lock = threading.Lock()
        self.reset()

    @property
    def n_pair(self):
        return 1

    def reset(self):
        self._lock.acquire()
        try:
            np.random.shuffle(self._node_list)
            self._counter = 0
            self._visited = set()
        finally:
            self._lock.release()
        return self

    def next(self):
        if (len(self._visited) == self._k) or (self._counter == len(self._node_list)):
            raise StopIteration
        self._lock.acquire()
        try:
            found, v = False, self._node_list[self._counter]
            while not found:
                v = self._node_list[self._counter]
                self._counter += 1
                if v not in self._visited:
                    found = all([not (self._g.has_edge(u, v) or self._g.has_edge(v, u)) for u in self._visited])
                if self._counter == len(self._node_list):
                    raise StopIteration
            self._visited.add(v)
        finally:
            self._lock.release()
        return v
