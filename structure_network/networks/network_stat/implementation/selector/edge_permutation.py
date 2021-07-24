import threading
import numpy as np
import networkx as nx
from structure_network.networks.network_stat.interfaces import NetworkGroupSelector

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['EdgePermutationSelector']


class EdgePermutationSelector(NetworkGroupSelector):
    def __init__(self, g, node_group, k=None):
        assert isinstance(g, nx.Graph)
        node_group = [n for n in node_group if n in g.nodes]
        assert len(node_group) > 1
        self._edge_list = []
        for u in node_group:
            for v in node_group:
                if g.has_edge(u, v):
                    self._edge_list.append((u, v))
        assert len(self._edge_list) > 0
        self._counter = 0
        self._k = len(self._edge_list) if k is None else min(k, len(self._edge_list))
        self._lock = threading.Lock()
        self.reset()

    @property
    def n_pair(self):
        return 2

    def reset(self):
        self._lock.acquire()
        try:
            np.random.shuffle(self._edge_list)
            self._counter = 0
        finally:
            self._lock.release()
        return self

    def next(self):
        self._lock.acquire()
        try:
            if self._counter == self._k:
                raise StopIteration
            u, v = self._edge_list[self._counter]
            self._counter += 1
        finally:
            self._lock.release()
        return u, v
