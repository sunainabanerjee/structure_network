import numpy as np
import networkx as nx
from structure_network.networks.core import to_nx
from structure_network.networks.core import to_list
from .node_centrality import NodeCentrality
from .node_centrality import node_centrality
from .graph_paths import shortest_path_distance
from structure_network.networks.core import add_edge_attribute

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['group_connectivity_pruning',
           'shortest_path_pruning',
           'all_path_pruning']


def group_connectivity_pruning(graph, source, target, **kwargs):
    """
    Prunes the network maximally without affecting all to all
    distances between two groups of nodes given a graph. The algorithm
    first computes all pair shortest paths between the two given node
    sets. Then removes the all the nodes that composes the shortest path
    and not part of the given node set. It recompute the shortest path
    and checks two groups remains connected. The algorithm iterates till there
    is no path between the node sets. Every node deletion is stored in a node
    queue. Finally an induced graph is produced by adding this node sets. The
    pruning process can be pre-empted by cutoff criteria like maximum path
    inflation. Path inflation criteria ensures worst connected path length of
    the pruned graph over starting configuration remains bounded by the inflation
    rate.
    """
    group1 = list(set(to_list(source)))
    group2 = list(set(to_list(target)))
    weight = kwargs.get('weight', None)
    return_view = kwargs.get('return_view', True)
    inflation = kwargs.get('max_inflation', 3.)
    drop_k = kwargs.get('drop_k', 10)
    drop_f = kwargs.get('drop_f', 1.0)
    assert drop_f > 0.1, "Error: fraction drop must be a positive number"
    min_connectivity = kwargs.get('min_connectivity', 1)
    if weight is not None:
        g = add_edge_attribute(graph, attribute_name=weight, **kwargs)
    else:
        g = to_nx(graph)
    assert all([g.has_node(n) for n in group1])
    assert all([g.has_node(n) for n in group2])
    assert len(set(group1).symmetric_difference(group2)) == len(group1) + len(group2)

    centrality = node_centrality(g, centrality=NodeCentrality.betweenness, weight=weight)

    def identify_prune_nodes(g_copy):
        paths, distances = shortest_path_distance(g_copy, **kwargs)
        distances = distances.loc[group1, group2]
        paths = paths.loc[group1, group2]
        node_q = []
        connections = 0
        for u in group1:
            increment = 0
            for v in group2:
                new_nodes = to_list(paths.loc[u, v])
                if len(new_nodes) > 0:
                    increment += 1
                node_q = node_q + new_nodes
            if increment > 0:
                connections += 1
        node_q = list(set(node_q) - set(group1 + group2))
        if np.isnan(distances.values).all():
            max_sep = np.inf
        else:
            max_sep = np.nanmax(distances.values)
        return sorted(node_q, key=lambda x: centrality[x]), max_sep, connections

    gc = g.copy()
    stop, start_separation = False, None
    node_set = set(group1 + group2)
    while not stop:
        prune_nodes, separation, connectivity = identify_prune_nodes(gc)
        if start_separation is None:
            start_separation = separation
        stop = (len(prune_nodes) == 0) or \
               ((separation / start_separation) > inflation) or \
               (connectivity < min_connectivity)
        if not stop:
            n_drop = min(drop_k, len(prune_nodes), max(1, int(len(prune_nodes)*drop_f)))
            for k in range(n_drop):
                gc.remove_node(prune_nodes[k])
            for n in prune_nodes:
                node_set.add(n)
    if return_view:
        return g.subgraph(node_set)
    return g.subgraph(node_set).copy()


def shortest_path_pruning(graph, source, target, **kwargs):
    group1 = list(set(to_list(source)))
    group2 = list(set(to_list(target)))
    weight = kwargs.get('weight', None)
    return_view = kwargs.get('return_view', True)
    inflation = kwargs.get('max_inflation', 3.)
    if weight is not None:
        g = add_edge_attribute(graph, attribute_name=weight, **kwargs)
    else:
        g = to_nx(graph)
    assert all([g.has_node(n) for n in group1])
    assert all([g.has_node(n) for n in group2])
    assert len(set(group1).symmetric_difference(group2)) == len(group1) + len(group2)
    node_set = set(group1 + group2)
    prune_set = set()
    min_distance = {}
    centrality = node_centrality(g,
                                 centrality=NodeCentrality.betweenness,
                                 normalized=False,
                                 weight=weight)
    for u in group1:
        min_distance[u] = {}
        for v in group2:
            min_distance[u][v] = None
            gc = g.copy(as_view=False)
            stop = False
            while not stop:
                new_prune = set()
                count = 0
                try:
                    path_generator = nx.all_shortest_paths(gc, u, v, weight=weight)
                    all_paths = list(path_generator)
                except nx.NetworkXNoPath:
                    stop = True
                    continue
                for path in all_paths:
                    cost = 0
                    count += 1
                    for ni in range(1, len(path)):
                        cost += gc[path[ni-1]][path[ni]][weight] \
                            if weight else 1
                        if path[ni] not in node_set:
                            new_prune.add(path[ni])
                    if min_distance[u][v] is None:
                        min_distance[u][v] = cost
                    if cost / min_distance[u][v] > inflation:
                        stop = True
                        break
                stop = stop or (len(new_prune) == 0)
                prune_set = prune_set.union(new_prune)
                if not stop:
                    new_prune = sorted(list(new_prune),
                                       key=lambda x: centrality[x])
                    for i in range(min(count, len(new_prune))):
                        gc.remove_node(new_prune[i])

    gc = g.subgraph(nodes=list(prune_set.union(node_set)))
    if return_view:
        return gc
    return gc.copy()


def all_path_pruning(graph, source, target, **kwargs):
    group1 = list(set(to_list(source)))
    group2 = list(set(to_list(target)))
    weight = kwargs.get('weight', None)
    return_view = kwargs.get('return_view', True)
    inflation = kwargs.get('max_inflation', 1.5)
    if weight is not None:
        g = add_edge_attribute(graph, attribute_name=weight, **kwargs)
    else:
        g = to_nx(graph)
    assert all([g.has_node(n) for n in group1])
    assert all([g.has_node(n) for n in group2])
    assert len(set(group1).symmetric_difference(group2)) == len(group1) + len(group2)

    distances, paths = {}, {}
    for n, (dist, path) in nx.all_pairs_dijkstra(g, weight=weight):
        distances[n] = {v: d for v, d in dist.items() if d > 0}
        paths[n] = {v: p for v, p in path.items() if dist[v] > 0}

    selected_nodes = set(group1).union(group2)
    remaining_node = set(g.nodes) - set(selected_nodes)
    for u in group1:
        assert u in distances
        for v in group2:
            assert v in distances[u]
            d = distances[u][v]
            for n in paths[u][v]:
                if n in remaining_node:
                    selected_nodes.add(n)
                    remaining_node.remove(n)
            if len(remaining_node) == 0:
                break
            node_list = list(remaining_node)
            node_len = len(node_list)
            counter, process_counter, accept_counter = 0, 0, 0
            while (process_counter + counter - accept_counter) < node_len:
                r = node_list[counter]
                if r in selected_nodes:
                    counter += 1
                    accept_counter += 1
                    continue
                if (distances[u][r] + distances[r][v]) < inflation * d:
                    for n in paths[u][r]:
                        if n not in selected_nodes:
                            selected_nodes.add(n)
                            remaining_node.remove(n)
                            process_counter += 1
                    for n in paths[r][v]:
                        if n not in selected_nodes:
                            selected_nodes.add(n)
                            remaining_node.remove(n)
                            process_counter += 1
                    accept_counter += 1
                counter += 1
        if len(remaining_node) == 0:
            break
    gc = g.subgraph(nodes=list(selected_nodes))
    if return_view:
        return gc
    return gc.copy()

