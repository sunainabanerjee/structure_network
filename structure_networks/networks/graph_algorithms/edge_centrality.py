import warnings
import networkx as nx
from enum import Enum
from structure_networks.networks.core import to_nx
from structure_networks.networks.core import add_edge_attribute

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['EdgeCentrality',
           'edge_centrality',
           'edge_pair_ordering',
           'validate_edge_centrality']


class EdgeCentrality(Enum):
    betweenness = 'betweenness'
    load = 'load'
    current_flow = 'current_flow'


def edge_pair_ordering(x, y):
    return (y, x) if x > y else (x, y)


def edge_centrality(graph, **kwargs):
    normalized = kwargs.get('normalized', True)
    weight = kwargs.get('weight', None)
    centrality = kwargs.get('centrality', EdgeCentrality.betweenness)
    seed = kwargs.get('seed', 121)
    solver = kwargs.get('solver', 'cg')
    if weight is not None:
        gc = add_edge_attribute(graph, attribute_name=weight, **kwargs)
    else:
        gc = to_nx(graph)

    edge_stat = {}
    try:
        if centrality == EdgeCentrality.betweenness:
            edge_stat = nx.edge_betweenness_centrality(gc,
                                                       normalized=normalized,
                                                       weight=weight,
                                                       seed=seed)
        elif centrality == EdgeCentrality.current_flow:
            edge_stat = nx.edge_current_flow_betweenness_centrality(gc,
                                                                    normalized=normalized,
                                                                    weight=weight,
                                                                    solver=solver)
        elif centrality == EdgeCentrality.load:
            cutoff = kwargs.get('cutoff', False)
            edge_stat = nx.edge_load_centrality(gc, cutoff=cutoff)
    except Exception as e:
        warnings.warn("Error while computing the edge "
                      "centraltiy [{}],".format(centrality.name) +
                      " defaulting to zero values")
        edge_stat = {(u, v): 0. for u, v in gc.edges}

    return {edge_pair_ordering(*key): value for key, value in edge_stat.items()}


def validate_edge_centrality(centrality_values):
    return not all([abs(v) < 1e-6 for k, v in centrality_values.items()])

