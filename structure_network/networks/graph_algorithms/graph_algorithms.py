import numpy as np
import pandas as pd
import networkx as nx
from structure_network.networks.core import to_nx
from .graph_paths import shortest_distance
from structure_network.networks.core import add_edge_attribute

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__version__ = "1.0.0"
__all__ = ['diameter',
           'disorient',
           'node_degree',
           'edge_inversion',
           'node_inversion',
           'is_connected',
           'connected_components']


def diameter(graph, **kwargs):
    distances = shortest_distance(graph, **kwargs)
    return np.nanmax(np.asarray(distances.values))


def node_degree(graph, **kwargs):
    weight = kwargs.get('weight', None)
    if weight is not None:
        g = add_edge_attribute(graph, attribute_name=weight, **kwargs)
    else:
        g = to_nx(graph)
    degrees = {}
    if weight is None:
        degrees = dict(nx.degree(g))
    else:
        for v in g.nodes:
            nbr = list(nx.neighbors(g, v))
            degrees[v] = np.sum([g.nodes[n][weight] for n in nbr])
    return degrees


def disorient(graph):
    g = to_nx(graph)
    if g.is_directed():
        ng = nx.Graph()
        for u, v in g.edges:
            ng.add_edge(u, v)
            ng[u][v].update(g[u][v])
        g = ng
    return g


def is_connected(graph):
    return nx.is_connected(to_nx(graph))


def edge_inversion(graph, **kwargs):
    g = to_nx(graph)
    attribute_name = kwargs.get('attribute', 'weight')
    max_value = kwargs.get('max_value', None)
    min_value = kwargs.get('min_value', None)
    attr_values = nx.get_edge_attributes(g, attribute_name)

    if len(attr_values) != g.size():
        raise RuntimeError("Error: attribute is not defined over all edges!")

    if max_value is None:
        max_value = np.max(list(attr_values.values()))
    if min_value is None:
        min_value = np.min(list(attr_values.values()))

    if max_value < min_value:
        raise RuntimeError("Error: max < min, edge inversion error!")

    gc = g.copy()
    for u, v in attr_values:
        transformed_value = (max_value - attr_values[(u, v)]) + min_value
        gc[u][v][attribute_name] = transformed_value
    return gc


def node_inversion(graph, **kwargs):
    g = to_nx(graph)
    attribute = kwargs.get('attribute')
    max_value = kwargs.get('max_value', None)
    min_value = kwargs.get('min_value', None)
    attr_values = nx.get_node_attributes(g, attribute)
    if len(attr_values) != g.order():
        raise RuntimeError("Error: attribute note defined on all nodes!")

    if min_value is None:
        min_value = np.min(list(attr_values.values()))
    if max_value is None:
        max_value = np.max(list(attr_values.values()))

    if max_value < min_value:
        raise RuntimeError("Error: max < min, node inversion error!")

    gc = g.copy()
    for u in attr_values:
        transformed_value = (max_value - attr_values[u]) + min_value
        gc.nodes[u][attribute] = transformed_value
    return gc


def connected_components(graph):
    return list(nx.connected_components(to_nx(graph)))


