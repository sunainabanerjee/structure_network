import networkx as nx
from structure_networks.networks.core import to_list
from .graph_extractor import InducedGraphExtractor
from .graph_extractor import energy_compliant_weight_inversion
from structure_networks.networks.core import add_edge_attribute
from structure_networks.networks.graph_algorithms import MarkovCluster
from structure_networks.networks.graph_algorithms import edge_inversion
from structure_networks.networks.graph_algorithms import is_connected
from structure_networks.networks.graph_algorithms import node_degree
from structure_networks.networks.graph_algorithms import NodeCentrality
from structure_networks.networks.graph_algorithms import EdgeCentrality
from structure_networks.networks.graph_algorithms import node_centrality
from structure_networks.networks.graph_algorithms import edge_centrality
from structure_networks.networks.graph_algorithms import all_path_pruning
from structure_networks.networks.graph_algorithms import shortest_path_pruning
from structure_networks.networks.graph_algorithms import group_connectivity_pruning
from structure_networks.networks.structure_network import ProteinStructureNetwork

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__version__ = "1.0"
__all__ = ['StructureNetworkAnalyzer']


class StructureNetworkAnalyzer:
    def __init__(self, structure_network):
        assert isinstance(structure_network, ProteinStructureNetwork)
        self.__structure_network = structure_network

    @property
    def nodes(self):
        return self.__structure_network.residue_names

    @property
    def residue_ids(self):
        return self.__structure_network.residue_ids

    @property
    def residue_names(self):
        return self.__structure_network.residue_names

    @property
    def threshold(self):
        return self.__structure_network.threshold

    @property
    def g(self):
        return self.__structure_network.g

    @property
    def structure(self):
        return self.__structure_network.structure

    @property
    def structure_network(self):
        return self.__structure_network

    @property
    def threshold_type(self):
        return self.__structure_network.threshold_type

    @property
    def params(self):
        return self.__structure_network.params

    def xyz(self, r=None):
        return self.__structure_network.xyz(r=r)

    @property
    def degree(self):
        return node_degree(self.__structure_network.g)

    def neighbors(self, r):
        if self.__structure_network.g.has_node(r):
            return list(nx.neighbors(self.__structure_network.g, r))
        return []

    def spanning_subgraph(self, source, target, **kwargs):
        subgraph_method = kwargs.get('method', InducedGraphExtractor.percolation)
        inflation = kwargs.pop('max_inflation', 1.5)
        weight = kwargs.pop('weight', None)
        weight_inversion = kwargs.pop('weight_inversion', False)

        if weight is not None:
            g = add_edge_attribute(self.g, attribute_name=weight, **kwargs)
        else:
            g = self.g
        if weight_inversion and (weight is not None):
            g = edge_inversion(g, attribute=weight)

        all_residues = self.residue_ids
        source = [x for x in to_list(source) if x in all_residues]
        target = [x for x in to_list(target) if x in all_residues]
        assert len(source) * len(target) > 0
        assert len(set(source).intersection(target)) == 0
        if subgraph_method == InducedGraphExtractor.percolation:
            cg = all_path_pruning(g, source, target, weight=weight, max_inflation=inflation, **kwargs)
        elif subgraph_method == InducedGraphExtractor.shortest_path:
            cg = shortest_path_pruning(g, source, target, weight=weight, max_inflation=inflation, **kwargs)
        else:
            cg = group_connectivity_pruning(g, source, target, weight=weight, max_inflation=inflation, **kwargs)
        return cg

    def node_centrality(self, **kwargs):
        centrality = kwargs.pop('centrality',
                                NodeCentrality.betweenness)
        weight = kwargs.get('weight',
                            'weight')
        inverse_weight = kwargs.pop('inverse_weight',
                                    energy_compliant_weight_inversion(centrality))

        g = self.__structure_network.g.copy()
        if weight is not None:
            g = add_edge_attribute(g, attribute_name=weight, **kwargs)

        if inverse_weight and (weight is not None):
            g = edge_inversion(g, attribute=weight)
        return node_centrality(g, centrality=centrality, **kwargs)

    def edge_centrality(self, **kwargs):
        centrality = kwargs.pop('centrality',
                                EdgeCentrality.betweenness)

        weight = kwargs.get('weight', 'weight')

        inverse_weight = kwargs.pop('inverse_weight',
                                    energy_compliant_weight_inversion(centrality))

        g = self.__structure_network.g.copy()
        if weight is not None:
            g = add_edge_attribute(g, attribute_name=weight, **kwargs)
        if inverse_weight and (weight is not None):
            g = edge_inversion(g, attribute=weight)
        return edge_centrality(g, centrality=centrality, **kwargs)

    def mark_centrality(self, **kwargs):
        centrality_list = self.node_centrality(**kwargs)
        structure = self.__structure_network.structure.copy()
        for res_id, value in centrality_list.items():
            structure.b_factor(resid=res_id, value=float(value))
        return structure

    def shortest_paths(self, u, v, **kwargs):
        res_ids = self.__structure_network.residue_ids
        res_names = self.__structure_network.residue_names
        weight = kwargs.get('weight', 'weight')
        method = kwargs.get('method', 'dijkstra')
        inverse_weight = kwargs.get('inverse_weight', True)
        g = self.__structure_network.g.copy()
        if weight is not None:
            g = add_edge_attribute(g, attribute_name=weight, **kwargs)
        if inverse_weight and (weight is not None):
            g = edge_inversion(g, attribute=weight)
        if u in res_ids:
            u = res_names[res_ids.index(u)]
        if v in res_ids:
            v = res_names[res_ids.index(v)]
        assert (u in res_names) and (v in res_names)
        if nx.has_path(g, source=u, target=v):
            return list(nx.all_shortest_paths(g,
                                              source=u,
                                              target=v,
                                              weight=weight,
                                              method=method))
        return []

    @property
    def connected(self):
        return is_connected(self.__structure_network.g)

    @property
    def structure_diameter(self):
        return self.__structure_network.structure_diameter

    @property
    def network_diameter(self):
        return self.__structure_network.network_diameter

    @property
    def centroid(self):
        return self.__structure_network.centroid

    def mcs(self, **kwargs):
        weight = kwargs.get('weight', 'weight')
        weight_inversion = kwargs.get('weight_inversion', False)
        g = self.__structure_network.g.copy()
        if weight is not None:
            g = add_edge_attribute(g, attribute_name=weight, **kwargs)
        if weight_inversion and (weight is not None):
            g = edge_inversion(g, attribute=weight)
        return MarkovCluster(g).run_mcs(**kwargs)

    def mark_mcs(self, **kwargs):
        return self.mark_group(self.mcs(**kwargs))

    def mark_group(self, residue_group):
        if isinstance(residue_group, list):
            residue_group = {(i+1): list(group_list)
                             for i, group_list in enumerate(residue_group)}
        assert isinstance(residue_group, dict)
        structure = self.structure.copy()
        for cls_id, residues in residue_group.items():
            for r in residues:
                if isinstance(r, str):
                    r = int(r[3:])
                structure.b_factor(resid=r, value=cls_id)
        return structure

