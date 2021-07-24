import warnings
import numpy as np
import pandas as pd
import networkx as nx
from ..core import to_nx
from ..core import to_list
from ..core import add_edge_attribute
from ..core import order_disruption_score
from .edge_centrality import edge_pair_ordering
from .node_centrality import validate_node_centrality
from .edge_centrality import validate_edge_centrality
from .node_centrality import NodeCentrality, node_centrality
from .edge_centrality import EdgeCentrality, edge_centrality

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['gomory_hu_cuts',
           'maxflow',
           'cut_size',
           'nodeflow_disruption_score',
           'edgeflow_disruption_score']


def gomory_hu_cuts(graph, **kwargs):
    capacity = kwargs.pop('capacity', 'weight')
    flow_func = kwargs.pop('flow_func', nx.algorithms.flow.shortest_augmenting_path)
    g = add_edge_attribute(graph, attribute_name=capacity, **kwargs)
    g_tree = nx.gomory_hu_tree(g, capacity=capacity, flow_func=flow_func)
    return [(u, v, g_tree[u][v]['weight']) for u, v in g_tree.edges]


def maxflow(graph, source, target, **kwargs):
    capacity = kwargs.get('capacity', 'weight')
    flow_func = kwargs.pop('flow_func', nx.algorithms.flow.shortest_augmenting_path)
    g = add_edge_attribute(graph, attribute_name=capacity, **kwargs)
    assert g.has_node(source) and g.has_node(target) and (source != target)
    flow_value, flow_dict = nx.maximum_flow(g, source, target,
                                            capacity=capacity,
                                            flow_func=flow_func)
    return flow_value, pd.DataFrame(flow_dict)


def cut_size(graph, source, target, **kwargs):
    capacity = kwargs.get('capacity', None)
    if capacity is not None:
        g = add_edge_attribute(graph, attribute_name=capacity, **kwargs)
    else:
        g = to_nx(graph)
    source = to_list(source)
    target = to_list(target)
    assert len(source) * len(target) > 0
    assert len(set(source).intersection(target)) == 0
    return nx.cut_size(G=g, S=source, T=target, weight=capacity)


def nodeflow_disruption_score(source_graph, target_graph, **kwargs):
    p = kwargs.pop('p', 3)
    beta = kwargs.pop('beta', 0)
    weight = kwargs.pop('weight', None)
    scaled = kwargs.pop('scaled', False)
    normalized = kwargs.pop('normalized', False)
    centrality = kwargs.pop('centrality', NodeCentrality.betweenness)
    node_map = kwargs.pop('node_map', lambda x: x)
    contribution = kwargs.pop('contribution', False)
    source_g = to_nx(source_graph)
    target_g = to_nx(target_graph)
    source_nodes = [node_map(n) for n in source_g.nodes]
    target_nodes = [node_map(n) for n in target_g.nodes]
    common_nodes = list(set(source_nodes).intersection(target_nodes))
    assert len(common_nodes) > 1
    centrality_source = node_centrality(source_g,
                                        centrality=centrality,
                                        weight=weight,
                                        normalized=normalized,
                                        **kwargs)
    centrality_target = node_centrality(target_g,
                                        centrality=centrality,
                                        weight=weight,
                                        normalized=normalized,
                                        **kwargs)
    if validate_node_centrality(centrality_source) and validate_edge_centrality(centrality_target):
        return order_disruption_score(centrality_source,
                                      centrality_target,
                                      node_map=node_map,
                                      beta=beta,
                                      scaled=scaled,
                                      p=p,
                                      contribution=contribution)
    warnings.warn('Error: received inappropriate centrality scores!')
    if contribution:
        common_keys = set(centrality_target.keys()).intersection(centrality_source.keys())
        contrib = {n: 0 for n in common_keys}
        return 0., contrib
    return 0.


def edgeflow_disruption_score(source_graph, target_graph, **kwargs):
    p = kwargs.pop('p', 3.)
    weight = kwargs.pop('weight', None)
    scaled = kwargs.pop('scaled', False)
    beta = kwargs.get('beta', 0.)
    normalized = kwargs.pop('normalized', False)
    centrality = kwargs.pop('centrality', EdgeCentrality.betweenness)
    node_map = kwargs.pop('node_map', lambda x: x)
    source_g = to_nx(source_graph)
    target_g = to_nx(target_graph)

    source_edges = [edge_pair_ordering(node_map(u), node_map(v)) for u, v in source_g.edges]
    target_edges = [edge_pair_ordering(node_map(u), node_map(v)) for u, v in target_g.edges]

    common_edges = list(set(source_edges).intersection(target_edges))

    assert len(common_edges) > 1
    centrality_source = edge_centrality(source_g,
                                        centrality=centrality,
                                        normalized=normalized,
                                        weight=weight,
                                        **kwargs)

    centrality_source = {edge_pair_ordering(node_map(key[0]), node_map(key[1])): value
                         for key, value in centrality_source.items()}

    centrality_target = edge_centrality(target_g,
                                        centrality=centrality,
                                        normalized=normalized,
                                        weight=weight,
                                        **kwargs)
    centrality_target = {edge_pair_ordering(node_map(key[0]), node_map(key[1])): value
                         for key, value in centrality_target.items()}

    if validate_edge_centrality(centrality_source) and validate_edge_centrality(centrality_target):
        return order_disruption_score(centrality_source,
                                      centrality_target,
                                      scaled=scaled,
                                      p=p,
                                      beta=beta)
    return 0.


