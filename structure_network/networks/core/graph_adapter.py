import logging
from .graph import *
import numpy as np
import networkx as nx

__all__ = ['to_nx',
           'from_nx',
           'add_edge_attribute',
           'add_node_attribute',
           'add_weight']


def to_nx(g):
    logger = logging.getLogger('to_nx')
    if isinstance(g, nx.Graph):
        return g.copy()

    if not isinstance(g, SimpleGraph):
        logger.error('Does not know how to handle the graph instance')
        raise RuntimeError('Invalid Graph Object!!')
    return g.g


def from_nx(g):
    logger = logging.getLogger('from_nx')
    if not isinstance(g, (nx.Graph, nx.DiGraph)):
        logger.error("Only simple undirected and directed graph "
                     "of networkx instances are supported")
        raise RuntimeError("Unsupported graph instance (%s)" % type(g).__class__.__name__)
    directed = g.is_directed()
    sg = SimpleGraph(directed=directed)
    for v in g.nodes:
        sg.add_node(v)
    for u, v in g.edges:
        sg.add_edge(u, v)
    return sg


def default_value(x):
    if isinstance(x, (np.float, np.int, complex)):
        return 1.
    if isinstance(x, dict):
        return {}
    if isinstance(x, (list, tuple)):
        return []
    if isinstance(x, str):
        return ''
    return None


def add_edge_attribute(graph, attribute_name, **kwargs):
    g = to_nx(graph)
    attr_values = nx.get_edge_attributes(g, attribute_name)
    if (len(attr_values) == 0) or np.isreal(list(attr_values.values())[0]):
        def_weight = 1
    else:
        example_value = list(attr_values.values())[0]
        def_weight = default_value(example_value)
    def_weight = kwargs.get('default_value', def_weight)
    def_fn = kwargs.pop('attribute_fun', None)
    for u, v in g.edges:
        if attribute_name not in g[u][v]:
            if def_fn is None:
                g[u][v][attribute_name] = def_weight
            else:
                g[u][v][attribute_name] = def_fn(g, u, v, **kwargs)
    return g


def add_node_attribute(graph, attribute_name, **kwargs):
    g = to_nx(graph)
    attr_values = nx.get_node_attributes(g, attribute_name)
    if len(attr_values) == 0:
        def_value = None
    else:
        def_value = default_value(list(attr_values.values())[0])
    def_value = kwargs.get('default_value', def_value)
    def_fn = kwargs.pop('attribute_fun', None)
    for u in g.nodes:
        if attribute_name not in g.nodes[u]:
            g.nodes[u][attribute_name] = def_value if def_fn is None else def_fn(g, u, **kwargs)
    return g


def add_weight(graph):
    return add_edge_attribute(graph, attribute_name='weight', default_value=1)
