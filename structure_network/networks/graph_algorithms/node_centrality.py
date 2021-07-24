import warnings
from enum import Enum
import networkx as nx
from structure_network.networks.core import to_nx
from structure_network.networks.core import add_edge_attribute
from structure_network.networks.core import add_node_attribute

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NodeCentrality',
           'node_centrality',
           'validate_node_centrality']


class NodeCentrality(Enum):
    betweenness = 'betweenness'
    load = 'load'
    information = 'information'
    communicability = 'communicability'
    percolation = 'percolation'
    eigen = 'eigen'
    current_flow = 'current_flow'
    current_flow_closeness = 'current_flow_closeness'
    harmonic = 'harmonic'
    subgraph = 'subgraph'
    degree = 'degree'
    closeness = 'closeness'
    katz = 'katz'
    second_order = 'second_order'


def node_centrality(graph, **kwargs):
    normalized = kwargs.get('normalized', True)
    weight = kwargs.get('weight', None)
    centrality = kwargs.get('centrality', NodeCentrality.betweenness)
    seed = kwargs.get('seed', 121)
    if weight is not None:
        gc = add_edge_attribute(graph, attribute_name=weight, **kwargs)
    else:
        gc = to_nx(graph)

    node_stat = {}
    try:
        if centrality == NodeCentrality.betweenness:
            node_stat = nx.betweenness_centrality(gc,
                                                  normalized=normalized,
                                                  weight=weight,
                                                  seed=seed)
        elif centrality == NodeCentrality.current_flow:
            node_stat = nx.current_flow_betweenness_centrality(gc,
                                                               normalized=normalized,
                                                               weight=weight)
        elif centrality == NodeCentrality.current_flow_closeness:
            node_stat = nx.current_flow_closeness_centrality(gc,
                                                             weight=weight)
        elif centrality == NodeCentrality.information:
            node_stat = nx.information_centrality(gc,
                                                  weight=weight)
        elif centrality == NodeCentrality.communicability:
            node_stat = nx.communicability_betweenness_centrality(gc)
        elif centrality == NodeCentrality.load:
            cutoff = kwargs.get('cutoff', None)
            node_stat = nx.load_centrality(gc,
                                           normalized=normalized,
                                           cutoff=cutoff,
                                           weight=weight)
        elif centrality == NodeCentrality.harmonic:
            node_stat = nx.harmonic_centrality(gc,
                                               distance=weight)
        elif centrality == NodeCentrality.eigen:
            tol = kwargs.get('tol', 1e-6)
            node_stat = nx.eigenvector_centrality_numpy(gc,
                                                        max_iter=gc.order(),
                                                        tol=tol,
                                                        weight=weight)
        elif centrality == NodeCentrality.katz:
            alpha = kwargs.get('alpha', 0.1)
            beta = kwargs.get('beta', 1.0)
            node_stat = nx.katz_centrality_numpy(gc,
                                                 alpha=alpha,
                                                 beta=beta,
                                                 weight=weight,
                                                 normalized=normalized)
        elif centrality == NodeCentrality.percolation:
            gc = add_node_attribute(gc,
                                    attribute_name='percolation',
                                    default_value=0.1)
            node_stat = nx.percolation_centrality(gc,
                                                  attribute='percolation',
                                                  weight=weight)
        elif centrality == NodeCentrality.subgraph:
            node_stat = nx.subgraph_centrality(gc)
        elif centrality == NodeCentrality.degree:
            node_stat = nx.degree_centrality(gc)
        elif centrality == NodeCentrality.closeness:
            wf_improved = kwargs.get('wf_improved', True)
            node_stat = nx.closeness_centrality(gc,
                                                distance=weight,
                                                wf_improved=wf_improved)
        elif centrality == NodeCentrality.second_order:
            node_stat = nx.second_order_centrality(gc)
    except Exception as e:
        warnings.warn("Error: computing the centrality {}, ".format(centrality.name) +
                      "defaulting the centrality values to zero")
        node_stat = {n: 0. for n in gc.nodes}
    return node_stat


def validate_node_centrality(centrality_values):
    return not all([abs(v) < 1e-6 for k, v in centrality_values.items()])
