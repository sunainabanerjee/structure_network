import warnings
import networkx as nx
import multiprocessing as mp
from structure_networks.ds import CacheData
from structure_networks.networks.core import to_list
from .graph_extractor import InducedGraphExtractor
from .graph_extractor import energy_compliant_weight_inversion
from structure_networks.networks.core import add_edge_attribute
from structure_networks.networks.graph_algorithms import EdgeCentrality
from structure_networks.networks.graph_algorithms import NodeCentrality
from structure_networks.networks.graph_algorithms import edge_inversion
from structure_networks.networks.graph_algorithms import node_centrality
from structure_networks.networks.graph_algorithms import edge_centrality
from structure_networks.networks.graph_algorithms import all_path_pruning
from structure_networks.networks.graph_algorithms import shortest_path_pruning
from structure_networks.networks.structure_network import ProteinStructureNetwork
from structure_networks.networks.graph_algorithms import group_connectivity_pruning


__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['SubstructureNetworkAnalyzer']


def f_node_centrality(group1: str,
                      group2: str,
                      g: nx.Graph,
                      centrality: NodeCentrality,
                      normalized: bool,
                      weight: str,
                      kw_args: dict):
    """
    Helper function to call node centrality computation for
    multiprocessing purpose
    """
    try:
        g_centrality = node_centrality(g,
                                       centrality=centrality,
                                       normalized=normalized,
                                       weight=weight,
                                       **kw_args)
    except RuntimeError:
        warnings.warn("Error: Failed to compute centrality "
                      "for ({}, {}, {}), defaulting".format(group1,
                                                            group2,
                                                            centrality.name))
        g_centrality = {n: 0 for n in g.nodes}
    return group1, group2, g_centrality


def f_edge_centrality(group1: str,
                      group2: str,
                      g: nx.Graph,
                      centrality: NodeCentrality,
                      normalized: bool,
                      weight: str,
                      kw_args: dict):
    """
    Helper function to call edge centrality computation for
    multiprocessing purpose.
    """
    try:
        e_centrality = edge_centrality(g,
                                       centrality=centrality,
                                       normalized=normalized,
                                       weight=weight,
                                       **kw_args)
    except RuntimeError:
        warnings.warn("Error: Failed to compute edge centrality for"
                      "({}, {}, {}), defaulting".format(group1,
                                                        group2,
                                                        centrality.name))
        e_centrality = {(u, v): 0.0 for u, v in g.edges}
    return group1, group2, e_centrality


def f_subgraph_extract(group1: str,
                       group2: str,
                       extractor: InducedGraphExtractor,
                       g: nx.Graph,
                       g_i: list,
                       g_j: list,
                       weight: str,
                       kw_args: dict):
    cg = nx.Graph()
    if extractor == InducedGraphExtractor.shortest_path:
        cg = shortest_path_pruning(g,
                                   source=g_i,
                                   target=g_j,
                                   weight=weight,
                                   return_view=False,
                                   **kw_args)
    elif extractor == InducedGraphExtractor.connectivity:
        cg = group_connectivity_pruning(g,
                                        source=g_i,
                                        target=g_j,
                                        weight=weight,
                                        return_view=False,
                                        **kw_args)
    elif extractor == InducedGraphExtractor.percolation:
        cg = all_path_pruning(g,
                              source=g_i,
                              target=g_j,
                              weight=weight,
                              return_view=False,
                              **kw_args)
    return group1, group2, cg


class SubstructureNetworkAnalyzer:
    def __init__(self, structure_network, **kwargs):
        assert isinstance(structure_network, ProteinStructureNetwork)
        self.__extractor = kwargs.pop('subgraph_method',
                                      InducedGraphExtractor.percolation)
        self.__inverse_weight = kwargs.pop('inverse_weight', True)
        self.__weight = kwargs.pop('weight', 'weight')
        self.__kwargs = kwargs
        self.__structure_network = structure_network
        self.__groups = []
        self.__group_name = []
        self.__subgraph = {}
        self.__caching = kwargs.get('caching', True)
        self.__cache = CacheData(allowed_keys=[
                                            ['node_centrality',
                                             'edge_centrality'], # type of centrality
                                            [True,
                                             False], # normalized
                                            [True,
                                             False], # inverse_weight
                                            ])

    def register_group(self, group, **kwargs):
        group_name = kwargs.get('name', int(len(self.__groups)))
        if group_name in self.__group_name:
            raise RuntimeError("Error: duplicate "
                               "group name {}".format(group_name))
        group = to_list(group)
        assert len(group) > 0
        r_names = self.__structure_network.residue_names
        r_ids = self.__structure_network.residue_ids
        group = [r_ids[r_names.index(n)] if n in r_names else n for n in group]
        assert all([n in r_ids for n in group])
        ok = True
        for grp in self.__groups:
            ok = ok and (len(set(group).intersection(grp)) == 0)
            if not ok:
                break
        if ok:
            self.__groups.append(group)
            self.__group_name.append(group_name)
        else:
            raise RuntimeError("Error: Not a non-overlapping node group!")
        self.__build()
        return self

    def __len__(self):
        return len(self.__groups)

    @property
    def weight_inversion(self):
        return self.__inverse_weight

    @property
    def registered_groups(self):
        return self.__group_name

    @property
    def residue_ids(self):
        return self.__structure_network.residue_ids

    @property
    def residue_names(self):
        return self.__structure_network.residue_names

    @property
    def g(self):
        return self.__structure_network.g

    @property
    def params(self):
        return self.__structure_network.params

    @property
    def structure(self):
        return self.__structure_network.structure

    @property
    def structure_network(self):
        return self.__structure_network

    @property
    def is_caching(self):
        return self.__caching

    def cache_on(self):
        self.__caching = True
        return self

    def cache_off(self):
        self.__caching = False
        self.reset_cache()
        return self

    def xyz(self, r=None):
        return self.__structure_network.xyz(r=r)

    @property
    def cover(self):
        group_nodes = set()
        for grp in self.__groups:
            group_nodes = group_nodes.union(grp)
        return list(group_nodes)

    def __getitem__(self, item):
        """
        List of residues in the specified
        residue group. In case while invalid
        group name provided returns None.
        """
        if item in self.__group_name:
            return self.__groups[self.__group_name.index(item)]
        return None

    def has_covered(self, r):
        return r in self.cover

    def group_contains(self, r):
        """
        Checks which residue group contains
        particular residue id, if true returns
        the residue group name, other wise returns
        None
        """
        if r in self.residue_names:
            r = self.residue_ids[self.residue_names.index(n)]
        for i, group in enumerate(self.__groups):
            if r in group:
                return self.__group_name[i]
        return None

    def common_residues(self, *args):
        """
        Returns residues common between different
        spanning sub-graphs. The all pair spanning
        sub-graphs are considered given a set of
        group ids.
        """
        groups = [g for g in list(set(args)) if g in self.__group_name]
        assert len(groups) > 1
        exclude_set = set()
        for g in groups:
            exclude_set = exclude_set.union(
                self.__groups[self.__group_name.index(g)])
        common_set = set(self.g.nodes) - exclude_set
        for i, g_i in enumerate(groups):
            for g_j in groups[i+1:]:
                common_set = common_set.intersection(
                    list(self.__subgraph[g_i][g_j].nodes))
        return list(common_set)

    def __build(self):
        """
        Internal function to build all pair subgraphs between
        all groups. Each subgraph spanning uses same algorithm
        and same parameters.
        """
        if len(self.__groups) > 1:
            n = len(self.__groups)
            g = self.__structure_network.g.copy()
            if self.__weight is not None:
                g = add_edge_attribute(g, attribute_name=self.__weight)
            if self.__inverse_weight and (self.__weight is not None):
                g = edge_inversion(g, attribute=self.__weight)

            r_names = self.residue_names
            r_ids = self.residue_ids

            result_obj = []
            pool = mp.Pool(mp.cpu_count())

            for i in range(n):
                name_i = self.__group_name[i]
                g_i = [r_names[r_ids.index(n)] for n in self.__groups[i]]
                for j in range(i+1, n):
                    name_j = self.__group_name[j]
                    if (name_i in self.__subgraph) and (name_j in self.__subgraph[name_i]):
                        continue
                    if name_j not in self.__subgraph:
                        self.__subgraph[name_j] = {}
                    g_j = [r_names[r_ids.index(n)]
                           for n in self.__groups[j]]
                    result_obj.append(pool.apply_async(f_subgraph_extract,
                                                       args=(name_i,
                                                             name_j,
                                                             self.__extractor,
                                                             g,
                                                             g_i,
                                                             g_j,
                                                             self.__weight,
                                                             self.__kwargs)))
            results = [r.get() for r in result_obj]
            pool.close()
            pool.join()
            for g1, g2, cg in results:
                if g1 not in self.__subgraph:
                    self.__subgraph[g1] = {}
                if g2 not in self.__subgraph:
                    self.__subgraph[g2] = {}
                self.__subgraph[g1][g2] = cg
                self.__subgraph[g2][g1] = cg
        return self

    def subgraph(self, g1, g2):
        """
        Returns spanning subgraph between two groups.
        """
        return self.__subgraph.get(g1, {}).get(g2, None)

    def subgraph_on_structure(self, g1, g2, **kwargs):
        mark = kwargs.get('mark', 50.0)
        other = kwargs.get('other', 0.0)
        skip_other = kwargs.get('skip_other', False)
        structure = self.__structure_network.structure.copy()
        sg = self.subgraph(g1, g2)
        if sg is None:
            raise RuntimeError("Error: group ({}, {}) is not registered".format(g1, g2))
        node_list = list(sg.nodes)
        r_names = self.residue_names
        r_ids = self.residue_ids
        for i, r_id in enumerate(r_ids):
            if r_names[i] in node_list:
                structure.b_factor(resid=r_id, value=mark)
            elif not skip_other:
                structure.b_factor(resid=r_id, value=other)
        return structure

    def reset_cache(self):
        self.__cache.clear()
        return self

    def group_node_centrality(self, g1, g2, **kwargs):
        normalized = kwargs.pop('normalized', False)
        weight = kwargs.pop('weight', 'weight')
        centrality_type = kwargs.pop('centrality',
                                     NodeCentrality.betweenness)
        inverse_weight = kwargs.pop('inverse_weight',
                                    energy_compliant_weight_inversion(centrality_type))
        centrality = {}
        if (g1 in self.__subgraph) or (g2 in self.__subgraph[g1]):
            if self.__caching:
                centrality = self.__cache.get(('node_centrality',
                                               normalized,
                                               inverse_weight,
                                               centrality_type))

            if g2 in centrality.get(g1, {}):
                return centrality[g1][g2]
            sg = self.subgraph(g1, g2)
            if inverse_weight != self.__inverse_weight:
                mx = self.__structure_network.max_edge_weight
                mn = self.__structure_network.min_edge_weight
                sg = edge_inversion(sg,
                                    attribute='weight',
                                    max_value=mx,
                                    min_value=mn)

            cent_computed = node_centrality(sg,
                                            centrality=centrality_type,
                                            normalized=normalized,
                                            weight=weight,
                                            **kwargs)
            if g1 not in centrality:
                centrality[g1] = {}
            centrality[g1][g2] = cent_computed.copy()

            if g2 not in centrality:
                centrality[g2] = {}
            centrality[g2][g1] = cent_computed.copy()

            if self.__caching:
                self.__cache.add(('node_centrality',
                                  normalized,
                                  inverse_weight,
                                  centrality_type),
                                 centrality)
            return centrality[g1][g2]
        return centrality

    def node_centrality(self, **kwargs):
        centrality = {}
        normalized = kwargs.pop('normalized', False)
        weight = kwargs.pop('weight', 'weight')
        centrality_type = kwargs.pop('centrality',
                                     NodeCentrality.betweenness)
        inverse_weight = kwargs.pop('inverse_weight',
                                    energy_compliant_weight_inversion(centrality_type))
        if self.__caching:
            centrality = self.__cache.get(('node_centrality',
                                           normalized,
                                           inverse_weight,
                                           centrality_type))

        pool = mp.Pool(mp.cpu_count())
        results_obj = []
        for i in self.__subgraph:
            for j in self.__subgraph[i]:
                if (i < j) and ((i not in centrality) or (j not in centrality[i])):
                    sg = self.subgraph(g1=i, g2=j)
                    if inverse_weight != self.__inverse_weight:
                        mx = self.__structure_network.max_edge_weight
                        mn = self.__structure_network.min_edge_weight
                        sg = edge_inversion(sg,
                                            attribute='weight',
                                            max_value=mx,
                                            min_value=mn)
                    results_obj.append(pool.apply_async(f_node_centrality,
                                                        args=(i,
                                                              j,
                                                              sg,
                                                              centrality_type,
                                                              normalized,
                                                              weight,
                                                              kwargs)))

        results = [r.get() for r in results_obj]
        pool.close()
        pool.join()

        for g1, g2, score in results:
            if g1 not in centrality:
                centrality[g1] = {}
            if g2 not in centrality:
                centrality[g2] = {}
            centrality[g1][g2] = score
            centrality[g2][g1] = score

        if self.__caching:
            self.__cache.add(['node_centrality',
                              normalized,
                              inverse_weight,
                              centrality_type],
                             centrality)
        return centrality

    def edge_centrality(self, **kwargs):
        centrality = {}
        normalized = kwargs.pop('normalized', False)
        weight = kwargs.pop('weight', 'weight')
        centrality_type = kwargs.pop('centrality',
                                     EdgeCentrality.betweenness)
        inverse_weight = kwargs.pop('inverse_weight',
                                    energy_compliant_weight_inversion(centrality_type))

        if self.__caching:
            centrality = self.__cache.get(('edge_centrality',
                                           normalized,
                                           inverse_weight,
                                           centrality_type))

        pool = mp.Pool(mp.cpu_count())
        results_obj = []
        for i in self.__subgraph:
            for j in self.__subgraph[i]:
                if (i < j) and ((i not in centrality) or (j not in centrality[i])):
                    sg = self.subgraph(g1=i, g2=j)
                    if inverse_weight != self.__inverse_weight:
                        mx = self.__structure_network.max_edge_weight
                        mn = self.__structure_network.min_edge_weight
                        sg = edge_inversion(sg,
                                            attribute='weight',
                                            max_value=mx,
                                            min_value=mn)
                    results_obj.append(pool.apply_async(f_edge_centrality,
                                                        args=(i,
                                                              j,
                                                              sg,
                                                              centrality_type,
                                                              normalized,
                                                              weight,
                                                              kwargs)))
        results = [r.get() for r in results_obj]
        pool.close()
        pool.join()

        for g1, g2, score in results:
            if g1 not in centrality:
                centrality[g1] = {}
            if g2 not in centrality:
                centrality[g2] = {}
            centrality[g1][g2] = score
            centrality[g2][g1] = score

        if self.__caching:
            self.__cache.add(('edge_centrality',
                              normalized,
                              inverse_weight,
                              centrality_type),
                             centrality)

        return centrality
