import numpy as np
from structure_network.networks.core import to_nx
from .node_centrality import node_centrality
from .graph_pruning import  group_connectivity_pruning

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['nodes_pair_centrality']


def nodes_pair_centrality(graph, group1, group2, **kwargs):
    """
    Nodes pair centrality explicitly measures centrality of the subgraph
    of union of paths two node sets. In order to achieve the same it prunes
    the graph till two graph all point connectivity is maintained. The analysis
    holds true only for undirected graph.
    """
    go = to_nx(graph)
    gc = group_connectivity_pruning(graph, group1=group1, group2=group2, **kwargs)
    node_stat = node_centrality(gc, **kwargs)
    for n in go.nodes:
        if n not in node_stat:
            node_stat[n] = np.nan
    return node_stat
