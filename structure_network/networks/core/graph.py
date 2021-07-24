"""
Basic utility for graph representation
"""

import numpy as np
import pandas as pd
import networkx as nx
from structure_network.geometry import *

__version__ = "1.0"
__all__ = ['SimpleGraph',
           'WeightedGraph',
           'GeometricGraph3d']


class SimpleGraph(object):
    def __init__(self, directed=False):
        self.__g = nx.DiGraph() if directed else nx.Graph()

    def add_node(self, v):
        if v not in self.__g.nodes:
            self.__g.add_node(v)
        return self

    def add_edge(self, v1, v2):
        self.__g.add_edge(v1, v2)
        return self

    def has_node(self, v):
        return self.__g.has_node(v)

    def is_vertex(self, v):
        return self.__g.has_node(v)

    def has_edge(self, v1, v2):
        return self.__g.has_edge(v1, v2)

    def del_edge(self, v1, v2):
        if self.has_edge(v1, v2):
            self.__g.remove_edge(v1, v2)
        return self

    def del_vertex(self, v):
        if self.is_vertex(v):
            self.__g.remove_node(v)
        return self

    @property
    def order(self):
        return self.__g.order()

    @property
    def size(self):
        return self.__g.size()

    @property
    def edges(self):
        return list(self.__g.edges)

    @property
    def vertices(self):
        return list(sorted(self.__g.nodes))

    @property
    def nodes(self):
        return self.vertices

    @property
    def adjacency(self):
        vs = self.vertices
        adj = pd.DataFrame(0, index=vs, columns=vs, dtype=np.int)
        for v1 in vs:
            for v2 in self.neighbors(v1):
                adj.loc[v1, v2] = 1
        return adj

    def in_degree(self, v):
        if not self.is_directed:
            raise RuntimeError("In degree not defined for undirected graph!")
        return self.__g.in_degree(v) if self.has_node(v) else 0

    def out_degree(self, v):
        if not self.is_directed:
            raise RuntimeError("Out degree not defined for undirected graph!")
        return self.__g.out_degree(v) if self.has_node(v) else 0

    def neighbors(self, v):
        if self.has_node(v):
            return list(self.__g.neighbors(v))
        return []

    def successors(self, v):
        if not self.is_directed:
            raise RuntimeError("Successor not defined for "
                               "undirected graph!")
        if self.has_node(v):
            return list(self.__g.successors(v))
        return []

    def predecessors(self, v):
        if not self.is_directed:
            raise RuntimeError("Predecessor not defined for "
                               "undirected graph!")
        if self.has_node(v):
            return list(self.__g.predecessors(v))
        return []

    def in_neighbors(self, v):
        return self.predecessors(v) if self.is_directed \
            else self.neighbors(v)

    def out_neighbors(self, v):
        return self.successors(v) if self.is_directed \
            else self.neighbors(v)

    @property
    def is_directed(self):
        return self.__g.is_directed()

    @property
    def g(self):
        return self.__g.copy()


class WeightedGraph(SimpleGraph):
    def __init__(self, directed=False, def_weight=1.0):
        super(self.__class__, self).__init__(directed=directed)
        self.__def_wt = def_weight
        self.__node_attrib = dict()
        self.__edge_attrib = dict()

    def add_node(self, v, attribute=None):
        super(self.__class__, self).add_node(v)
        self.__node_attrib[v] = attribute
        return self

    def add_edge(self, u, v, weight=None):
        super(self.__class__, self).add_edge(u, v)
        if u not in self.__edge_attrib:
            self.__edge_attrib[u] = dict()
        self.__edge_attrib[u][v] = weight
        if not self.is_directed:
            if v not in self.__edge_attrib:
                self.__edge_attrib[v] = dict()
            self.__edge_attrib[v][u] = weight
        return self

    def del_edge(self, u, v):
        super(self.__class__, self).del_edge(u, v)
        del self.__edge_attrib[u][v]
        if not self.is_directed:
            del self.__edge_attrib[v][u]
        return self

    def del_vertex(self, v):
        if self.has_node(v):
            super(self.__class__, self).del_vertex(v)
            del self.__node_attrib[v]
            if self.is_directed:
                for u in self.in_neighbors(v):
                    del self.__edge_attrib[u][v]
                del self.__edge_attrib[v]
            else:
                for u in self.neighbors(v):
                    del self.__edge_attrib[u][v]
                del self.__edge_attrib[v]
        return self

    def edge_attribute(self, u, v):
        if (u in self.__edge_attrib) and (v in self.__edge_attrib[u]):
            return self.__edge_attrib[u][v]
        raise RuntimeError("Not a valid edge (%s -> %s)" % (str(u), str(v)))

    def node_attribute(self, v):
        if v in self.__node_attrib:
            return self.__node_attrib[v]
        return None

    def __getitem__(self, item):
        if self.is_vertex(item):
            return self.node_attribute(item)
        elif len(item) == 2:
            return self.edge_attribute(item[0], item[1])
        raise RuntimeError("Requested item is neither a vertex nor an edge")

    @property
    def weight_matrix(self):
        vs = self.vertices
        wts = pd.DataFrame(0, index=vs, columns=vs)
        for vi in vs:
            for vj in self.neighbors(vi):
                wts.loc[vi, vj] = self.__edge_attrib[(vi, vj)]
        return wts

    @property
    def g(self):
        ng = nx.DiGraph() if self.is_directed else nx.Graph()
        for v in self.nodes:
            ng.add_node(v, attr=self.node_attribute(v))
        for u, v in self.edges:
            ng.add_edge(u, v, weight=self.edge_attribute(u, v))
        return ng


class GeometricGraph3d(WeightedGraph):
    def __init__(self, def_weight=1.0):
        super(self.__class__, self).__init__(directed=False,
                                             def_weight=def_weight)

    def add_vertex(self, v, attribute):
        if not isinstance(attribute, Coordinate3d):
            raise RuntimeError("Support Coordinate3d as node attribute!")
        super(self.__class__, self).add_node(v, attribute=attribute)

    def add_edge(self, u, v, weight=None):
        assert self.is_vertex(u)
        assert self.is_vertex(v)
        if weight is None:
            weight = distance(self.node_attribute(u), self.node_attribute(v))
        super(GeometricGraph3d, self).add_edge(u, v, weight=weight)

    def weight(self, n_1, n):
        pass

    def attribute(self, n):
        pass

