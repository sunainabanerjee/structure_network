import pandas as pd
import networkx as nx
from structure_networks.structure import Mutation
from structure_networks.networks.core import to_list
from structure_networks.networks.core import add_edge_attribute
from structure_networks.networks.network_stat import SupportFunction
from structure_networks.networks.graph_algorithms import edge_inversion
from structure_networks.networks.structure_network import ProteinStructureNetwork

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['preprocess_structure_network', 'mutation_profile']


def preprocess_structure_network(protein_network, **kwargs):
    assert isinstance(protein_network, ProteinStructureNetwork)
    invert_edges = kwargs.get('edge_inversion', False)
    edge_attr = kwargs.get('edge_attribute', 'weight')
    force_update = kwargs.get('force_update', False)
    support_functions = to_list(kwargs.get('support_functions', []))
    assert all([isinstance(fn, SupportFunction) for fn in support_functions])
    g = protein_network.g
    r_names = protein_network.residue_names
    r_ids = protein_network.residue_ids
    assert len(r_names) == len(r_ids)
    assert len(r_names) == g.order()
    node_map = {n: r_ids[i]
                for i, n in enumerate(r_names)}
    gc = nx.Graph()

    for n, v in node_map.items():
        gc.add_node(v, **g.nodes[n])

    for u, v in g.edges:
        gc.add_edge(node_map[u], node_map[v], **g[u][v])
        for supp_func in support_functions:
            if force_update or (supp_func.name not in gc[node_map[u]][node_map[v]]):
                gc[node_map[u]][node_map[v]][supp_func.name] = supp_func(u, v)

    if invert_edges and (edge_attr is not None):
        gc = add_edge_attribute(gc, attribute_name=edge_attr)
        gc = edge_inversion(gc, attribute=edge_attr)

    return gc


def mutation_profile(mutations, delta_values):
    assert len(mutations) == len(delta_values)
    assert all([isinstance(m, Mutation) for m in mutations]), \
        "Error: expect mutation instances only!"
    assert len(mutations) > 1
    effects = {'residue': [], 'mutation': [], 'effect': []}
    for i, m in enumerate(mutations):
        effects['residue'].append(m.residue_id)
        effects['mutation'].append(m.mutation)
        effects['effect'].append(delta_values[i])
    effect = pd.DataFrame(effects)
    return pd.pivot_table(effect,
                          values='effect',
                          index=['residue'],
                          columns=['mutation'])

