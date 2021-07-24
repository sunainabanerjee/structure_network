import threading
import numpy as np
import networkx as nx
from structure_network.networks.network_stat.interfaces import NetworkGroupSelector

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NodePermutationSelector']


class NodePermutationSelector(NetworkGroupSelector):
    def __init__(self, g, node_group, k=None):
        assert isinstance(g, nx.Graph)
        self._g = g
        self._nodes = [n for n in node_group if n in self._g.nodes]
        self._counter = 0
        self._k = len(self._nodes) if k is None else min(k, len(self._nodes))
        self._lock = threading.Lock()
        self.reset()

    @property
    def n_pair(self):
        return 1

    def reset(self):
        self._lock.acquire()
        try:
            np.random.shuffle(self._nodes)
            self._counter = 0
        finally:
            self._lock.release()
        return self

    def next(self):
        self._lock.acquire()
        try:
            if self._counter == self._k:
                raise StopIteration
            val = self._nodes[self._counter]
            self._counter += 1
        finally:
            self._lock.release()
        return val
