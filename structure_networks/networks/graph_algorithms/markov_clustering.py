import logging
import numpy as np
import networkx as nx
from ..core import to_nx
import markov_clustering as mcs
from ..core import add_edge_attribute

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['MarkovCluster']


class MarkovCluster:
    def __init__(self, g):
        assert isinstance(g, nx.Graph)
        self.__g = g
        self.__logger = logging.getLogger(self.__class__.__name__)

    @property
    def __len__(self):
        return self.__g.order()

    def order(self):
        return self.__g.order()

    def size(self, **kwargs):
        return self.__g.size(**kwargs)

    def hpo(self, **kwargs):
        min_inflation = kwargs.get('min_inflation', 1.5)
        max_inflation = kwargs.get('max_inflation', 2.5)
        bins = kwargs.get('bins', 10)
        debug = kwargs.get('verbose', False)
        weight = kwargs.get('weight', None)
        assert (max_inflation > min_inflation) and (min_inflation > 0)
        if weight is not None:
            g = add_edge_attribute(self.__g, attribute_name=weight, **kwargs)
        else:
            g = to_nx(self.__g)
        matrix = nx.to_scipy_sparse_matrix(g)
        step = (max_inflation - min_inflation) / (bins + 1)
        inflation_values = np.arange(start=min_inflation, stop=max_inflation, step=step)
        modularity = list()
        for x in inflation_values:
            result = mcs.run_mcl(matrix, inflation=x)
            clusters = mcs.get_clusters(result)
            modularity.append(mcs.modularity(matrix=result, clusters=clusters))
            if debug:
                self.__logger.info("At inflation: {:.3f} number of "
                                   "clusters found {} with modularity "
                                   "score: {}".format(x, len(clusters), modularity[-1]))
        return inflation_values[np.argmax(modularity)]

    def run_mcs(self, **kwargs):
        inflation = kwargs.get('inflation', None)
        modularity = kwargs.get('compute_modularity', False)
        if inflation is None:
            inflation = self.hpo(**kwargs)

        weight = kwargs.get('weight', None)
        if weight is not None:
            g = add_edge_attribute(self.__g, attribute_name=weight, **kwargs)
        else:
            g = to_nx(self.__g)
        matrix = nx.to_scipy_sparse_matrix(g, weight=weight)

        result = mcs.run_mcl(matrix, inflation=inflation)
        clusters = mcs.get_clusters(result)
        if modularity:
            return mcs.modularity(matrix=result, clusters=clusters)
        node_names = list(g.nodes)
        return [[node_names[index] for index in cluster] for cluster in clusters]
