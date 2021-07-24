import threading
import numpy as np
import networkx as nx
from .edge_permutation import EdgePermutationSelector

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['EdgeCoverSelector']


class EdgeCoverSelector(EdgePermutationSelector):
    def __init__(self, g, **kwargs):
        assert isinstance(g, nx.Graph)
        node_group = kwargs.get('nodes', list(g.nodes))
        k = kwargs.get('k', g.size())
        super(EdgeCoverSelector, self).__init__(g, node_group, k)
        self._visited = set()
        self._counter = 0
        self._child_lock = threading.Lock()
        self.reset()

    def reset(self):
        self._child_lock.acquire()
        try:
            np.random.shuffle(self._edge_list)
            self._counter = 0
            self._visited = set()
        finally:
            self._child_lock.release()
        return self

    def next(self):
        self._child_lock.acquire()
        try:
            if (self._counter == len(self._edge_list)) or (len(self._visited)//2 >= self._k):
                raise StopIteration
            found = False
            u, v = None, None
            while not found:
                u, v = self._edge_list[self._counter]
                self._counter += 1
                found = (u not in self._visited) and (v not in self._visited)
                if self._counter == len(self._edge_list):
                    raise StopIteration
            self._visited.add(u)
            self._visited.add(v)
        finally:
            self._child_lock.release()
        return u, v


