import logging
import numpy as np
from .graph import *
from copy import deepcopy
from structure_network.geometry import *

__version__ = "1.0"
__all__ = ['PathFilter',
           'GeometricPathFilter',
           'AllPathIterator']


class PathFilter:
    def __init__(self):
        pass

    def erase_memory(self):
        pass

    def __call__(self, g, n, n_1, n_2=None):
        return True


class GeometricPathFilter(PathFilter):
    def __init__(self, **kwargs):
        self.__logger = logging.getLogger(name='temporal_graph.GeometricPathFilter')
        PathFilter.__init__(self)
        self.__forward_only = kwargs.get('forward', True)
        self.__edge_retrace = kwargs.get('edge_retrace', True)
        self.__edge_memory = dict()
        self.__weight_cutoff = kwargs.get('weight_cutoff', None)
        distance_cutoff = kwargs.get('distance_cutoff', None)
        if distance_cutoff is not None:
            if isinstance(distance_cutoff, DistanceCutoff):
                self._distance_cutoff = distance_cutoff
            elif isinstance(distance_cutoff, np.float) or \
                    isinstance(distance_cutoff, np.double) or \
                    isinstance(distance_cutoff, np.int):
                self._distance_cutoff = DistanceCutoff(def_cutoff=distance_cutoff)
            self.__logger.debug('Setting distance cutoff (%f)!!' % self._distance_cutoff.cutoff)
        else:
            self._distance_cutoff = None

    @property
    def is_forward_only(self):
        return self.__forward_only

    @property
    def has_weight_cutoff(self):
        return self.__weight_cutoff is not None

    @property
    def weight_cutoff(self):
        return self.__weight_cutoff

    @property
    def has_distance_cutoff(self):
        return self._distance_cutoff is not None

    @property
    def distance_cutoff(self):
        return self._distance_cutoff

    def erase_memory(self):
        self.__edge_memory = dict()

    def __call__(self, g, n, n_1, n_2=None):
        if not isinstance(g, GeometricGraph3d):
            self.__logger.debug("Not a geometric graph!")
            return False
        if not g.is_vertex(n) or not g.is_vertex(n_1) or not g.has_edge(n_1, n):
            return False
        if (n_2 is not None) and (not g.is_vertex(n_2) or not g.has_edge(n_2, n_1)):
            return False
        if self.has_weight_cutoff and g.weight(n_1, n) < self.weight_cutoff:
            return False
        if self.has_distance_cutoff:
            d = distance(g.attribute(n), g.attribute(n_1))
            if d > self._distance_cutoff():
                return False
        if self.is_forward_only and n_2 is not None:
            v1 = connecting_vector(g.attribute(n_1), g.attribute(n)).unit_vector
            v2 = connecting_vector(g.attribute(n_2), g.attribute(n_1)).unit_vector
            if dotp(v1, v2) > 0:
                return False

        if n_1 not in self.__edge_memory:
            self.__edge_memory[n_1] = dict()
        if n not in self.__edge_memory[n_1]:
            self.__edge_memory[n_1][n] = True
        else:
            return self.__edge_retrace
        return True


def path_distance(path1, path2):
    assert isinstance(path1, list) and isinstance(path2, list)
    return len(set(path1).intersection(set(path2))) * 1.0/len(set(path1).union(set(path2)))


class AllPathIterator:
    def __init__(self,
                 g,
                 visitor,
                 path_separation=0.4):
        assert isinstance(g, WeightedGraph)
        assert isinstance(visitor, PathFilter)
        self.__logger = logging.getLogger(name="temporal_graph.AllPathIterator")
        self.__g = deepcopy(g)
        self.__visitor = visitor
        self.__nodes = {v: False for v in g.vertices}
        self.__stack = []
        self.__paths = []
        self.__stop_vertex = set()
        self.__visited = set()
        self.__retrack = False
        self.__retrack_length = 0
        self.__path_separation = path_separation
        self.__logger.debug('AllPathIterator instantiated!!')

    def is_vertex(self, v):
        return self.__g.is_vertex(v)

    def all_path(self,
                 start_vertex,
                 stop_vertex,
                 max_paths=50,
                 min_path_length=2,
                 max_path_length=None,
                 retrack_fraction=0):
        assert self.__g.is_vertex(start_vertex)
        g_bk = None
        if isinstance(self.__visitor, GeometricPathFilter) and self.__visitor.has_weight_cutoff:
            g_bk = deepcopy(self.__g)
            cutoff = self.__visitor.weight_cutoff
            remove_edges = []
            for u, v in self.__g.edges:
                if self.__g.weight(u, v) < cutoff:
                    remove_edges.append((u, v))
            self.__logger.debug("Number of edges to be deleted (%d)" % len(remove_edges))
            self.__logger.debug("Original graph dimension G(%d, %d)" % (self.__g.order, self.__g.size))
            for u, v in remove_edges:
                self.__g.del_edge(u, v)
            self.__logger.debug("Post edge deletion graph size G(%d, %d)" % (self.__g.order, self.__g.size))

        if isinstance(stop_vertex, list):
            self.__stop_vertex = set(stop_vertex)
        else:
            self.__stop_vertex = {stop_vertex}
        assert start_vertex not in self.__stop_vertex
        for s in self.__stop_vertex:
            assert self.__g.is_vertex(s)
        if max_path_length is None:
            max_path_length = self.__g.order // 2

        self.__paths = []
        has_path = False
        for v in self.__stop_vertex:
            if self.__g.is_connected(start_vertex, v):
                has_path = True
                self.__logger.debug("Found path between (%s) <--> (%s)" % (start_vertex, v))
                break
        if has_path:
            self.__visited = set()
            self.__logger.debug("Evaluating all paths for residue: %s" % start_vertex)
            self.__retrack = False
            self.__retrack_length = 0
            self.__all_paths(start_vertex,
                             max_paths=max_paths,
                             min_path_length=min_path_length,
                             max_path_length=max_path_length,
                             retrack_fraction=retrack_fraction)
            paths = deepcopy(self.__paths)

        if g_bk is not None:
            self.__g = deepcopy(g_bk)
            self.__logger.debug("Restoring graph to size G(%d, %d)" % (self.__g.order, self.__g.size))
        self.__logger.debug("Number of paths found for vertex (%s) is %d" % (start_vertex, len(paths)))
        self.__paths = []
        return paths

    def __all_paths(self,
                    start_vertex,
                    max_paths,
                    min_path_length=2,
                    max_path_length=10,
                    retrack_fraction=0.5):
        assert self.__g.is_vertex(start_vertex)
        if (not self.__nodes[start_vertex]) and (len(self.__paths) < max_paths):
            append = False
            if len(self.__stack) > 0:
                v = start_vertex
                v_1 = self.__stack[-1]
                v_2 = self.__stack[-2] if len(self.__stack) > 1 else None
                if self.__visitor(self.__g, v, v_1, v_2):
                    append = True
            else:
                append = True
            if append:
                self.__stack.append(start_vertex)
                self.__nodes[start_vertex] = True
                if start_vertex in self.__stop_vertex:
                    if len(self.__stack) >= min_path_length:
                        path = deepcopy(self.__stack)
                        reject = True
                        for p in path:
                            if p not in self.__visited:
                                self.__visited.add(p)
                                reject = False
                        if reject:
                            accept = True
                            for p in self.__paths:
                                if path_distance(path, p['path']) < self.__path_separation:
                                    accept = False
                            reject = not accept
                        if not reject:
                            path_dict = {'path': path[:],
                                         'weights': []}
                            for i in range(1, len(self.__stack)):
                                path_dict['weights'].append(self.__g.weight(self.__stack[i-1], self.__stack[i]))
                            self.__paths.append(path_dict)
                            path_len = int( len(path) * retrack_fraction )
                            if path_len > 1:
                                self.__retrack = True
                                self.__retrack_length = path_len
                            self.__logger.debug("Path found between (%s <-> %s) of length %d" % (path[0], path[-1], len(path)))
                elif len(self.__stack) < max_path_length:
                    for v in self.__g.out_neighbors(start_vertex):
                        self.__all_paths(v,
                                         max_paths=max_paths,
                                         min_path_length=min_path_length,
                                         max_path_length=max_path_length,
                                         retrack_fraction=retrack_fraction)
                        if self.__retrack:
                            if self.__retrack_length > 0:
                                self.__retrack_length -= 1
                            else:
                                self.__retrack = False
                                self.__retrack_length = 0
                        if self.__retrack:
                            break
                n = self.__stack.pop(-1)
                self.__nodes[start_vertex] = False


