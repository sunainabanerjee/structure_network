import numpy as np
import pandas as pd
import networkx as nx
from structure_network.networks.core import to_nx
from structure_network.networks.core import to_list
from structure_network.networks.core import add_edge_attribute

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['shortest_path',
           'shortest_distance',
           'all_shortest_paths',
           'shortest_path_distance']


def shortest_path(graph, **kwargs):
    weight = kwargs.get('weight', None)
    if weight is not None:
        g = add_edge_attribute(graph, attribute_name=weight, **kwargs)
    else:
        g = to_nx(graph)
    source = to_list(kwargs.get('source', list(g.nodes)))
    target = to_list(kwargs.get('target', list(g.nodes)))

    if not all([g.has_node(n) for n in source]):
        raise RuntimeError("Error: requested nodes [{}] not in the graph!".format(source))
    if not all([g.has_node(n) for n in target]):
        raise RuntimeError("Error: requested nodes [{}] not in the graph!".format(target))
    path_df = pd.DataFrame(nx.shortest_path(g, weight=weight))
    return path_df.loc[source, target]


def shortest_distance(graph, **kwargs):
    weight = kwargs.get('weight', None)
    method = kwargs.get('method', 'dijkstra')
    if weight is None:
        g = to_nx(graph)
    else:
        g = add_edge_attribute(graph, attribute_name=weight, **kwargs)

    source = to_list(kwargs.get('source', list(g.nodes)))
    target = to_list(kwargs.get('target', list(g.nodes)))
    assert (len(source) > 0) and (len(target) > 0)
    result = pd.DataFrame(dict(nx.shortest_path_length(g,
                                                       weight=weight,
                                                       method=method)))
    return result.loc[source, target]


def all_shortest_paths(graph, src, tgt, **kwargs):
    weight = kwargs.get('weight', None)
    if weight is not None:
        g = add_edge_attribute(graph, attribute_name=weight, **kwargs)
    else:
        g = to_nx(graph)
    if not g.has_node(src):
        raise RuntimeError("Error: source node {} is not a part of the graph".format(src))
    if not g.has_node(tgt):
        raise RuntimeError("Error: target node {} is not a part of the graph".format(tgt))

    min_length = kwargs.get('min_path_length', 2)
    max_length = kwargs.get('max_path_length', g.order())
    shortest_paths = []
    for path in nx.all_shortest_paths(g, source=src, target=tgt, weight=weight):
        if (len(path) >= min_length) and (len(path) <= max_length):
            shortest_paths.append(path)
    return shortest_paths


def shortest_path_distance(graph, **kwargs):
    weight = kwargs.get('weight', None)
    cutoff = kwargs.get('cutoff', np.inf)
    if weight is not None:
        g = add_edge_attribute(graph, attribute_name=weight, **kwargs)
    else:
        g = to_nx(graph)
    source = to_list(kwargs.get('source', list(g.nodes)))
    target = to_list(kwargs.get('target', list(g.nodes)))
    assert (len(source) > 0) and (len(target) > 0)
    assert all([g.has_node(n) for n in source]), "Error: invalid source node list provided!"
    assert all([g.has_node(n) for n in target]), "Error: invalid target node list provided!"
    all_path = dict(nx.all_pairs_dijkstra(g, cutoff=cutoff, weight=weight))
    vertices = list(g.nodes)
    path_weight, path_trace = [], []
    for u in vertices:
        row_weight, row_trace = [], []
        for v in vertices:
            if (u not in all_path) or (v not in all_path[u][0]):
                row_weight.append(np.nan)
                row_trace.append(None)
            else:
                row_weight.append(all_path[u][0][v])
                row_trace.append(all_path[u][1][v])
        path_weight.append(row_weight)
        path_trace.append(row_trace)
    path_weight = pd.DataFrame(path_weight, index=vertices, columns=vertices)
    path_trace = pd.DataFrame(path_trace, index=vertices, columns=vertices)
    return path_trace.loc[source, target], path_weight.loc[source, target]
