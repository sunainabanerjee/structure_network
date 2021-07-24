import pandas as pd
import networkx as nx
from structure_network.networks.core import to_list
from .network_metric import NetworkMetric, GroupStatistics
from structure_network.networks.graph_algorithms import EdgeCentrality
from structure_network.networks.graph_algorithms import edge_inversion
from structure_network.networks.graph_algorithms import edge_centrality
from structure_network.networks.graph_algorithms import edge_pair_ordering

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['DeltaEdgeCentrality']


class DeltaEdgeCentrality(NetworkMetric):
    def __init__(self, **kwargs):
        self._g = kwargs.get('g', None)
        self._statistics = kwargs.get('statistics', GroupStatistics.average)
        self._weight = kwargs.get('weight', 'weight')
        self._normalized = kwargs.get('normalized', False)
        self._weight_inversion = kwargs.get('weight_inversion', False)
        self._centrality = kwargs.get('centrality', EdgeCentrality.betweenness)
        def_node_group = list(self._g.nodes) if isinstance(self._g, nx.Graph) else []
        self._node_group = to_list(kwargs.get('nodes', def_node_group))
        self._edge_group = []
        self._check_status = False
        self._check_ready = False
        self._wild_score = {}
        if self.is_ready:
            self._update_reference()

    def _update_reference(self):
        if not self.is_ready:
            raise RuntimeError("Error: update_reference called prior to network setup!")
        self._update_edge_group()
        if self._weight_inversion:
            g = edge_inversion(graph=self._g, attribute=self._weight)
        else:
            g = self._g.copy()
        self._wild_score = edge_centrality(graph=g,
                                           centrality=self._centrality,
                                           weight=self._weight,
                                           normalized=self._normalized)
        return self

    @property
    def is_ready(self):
        if not self._check_status:
            self._check_ready = isinstance(self._g, nx.Graph)
            if self._check_ready:
                all_attrib = nx.get_edge_attributes(self._g, name=self._weight)
                self._check_ready = len(all_attrib) == self._g.size()
            if self._check_ready:
                self._check_ready = len(self._node_group) > 0
            self._check_status = True
        return self._check_ready

    def _update_edge_group(self):
        if not self.is_ready:
            raise RuntimeError("Error: edge group can not be updated "
                               "till network is fully configured!")
        self._edge_group = []
        for u, v in self._g.edges:
            if (u in self._node_group) or (v in self._node_group):
                self._edge_group.append(edge_pair_ordering(u, v))
        return self

    @property
    def g(self):
        return self._g

    @property
    def name(self):
        return 'edge-' + self._centrality.name + \
               '-' + self._statistics.name + \
               '-' + str(self._normalized)

    @property
    def target_group(self):
        return self._edge_group

    def set_target(self, group):
        if self._g is None:
            raise RuntimeError("Error: target can not be set with out network setup!")
        group = to_list(group)
        for x in group:
            if not self._g.has_node(x):
                raise RuntimeError("Error: invalid node specification {}!".format(x))
        self._node_group = group
        self._check_status = False
        if self._check_status:
            self._update_reference()
        return self

    def set_base_network(self, g):
        if isinstance(g, nx.Graph):
            self._g = g
            self._check_status = False
            if len(self._node_group) > 0:
                self._node_group = [n for n in self._node_group if self._g.has_node(n)]
                if len(self._node_group) == 0:
                    raise RuntimeError("Error: target node group setting is improper!")
            else:
                self._node_group = list(self._g.nodes)

            if self.is_ready:
                self._update_reference()
            return self
        raise RuntimeError("Error: expects network Graph object!")

    def centrality_change(self, g, **kwargs):
        assert isinstance(g, nx.Graph)
        assert len(set(g.nodes).symmetric_difference(set(self._g.nodes))) == 0, \
            "Error: graphs are different!"
        if self._weight_inversion:
            gc = edge_inversion(g, attribute='weight')
        else:
            gc = g.copy()
        score = edge_centrality(gc,
                                centrality=self._centrality,
                                weight=self._weight,
                                normalized=self._normalized)
        df = pd.DataFrame({'base': self._wild_score,
                           'perturbed': score})
        return df.loc[self._edge_group, 'perturbed'] - df.loc[self._edge_group, 'base']

    def __call__(self, g, **kwargs):
        df = self.centrality_change(g, **kwargs)
        if self._statistics == GroupStatistics.average:
            return float(df.mean())
        if self._statistics == GroupStatistics.median:
            return float(df.median())
        if self._statistics == GroupStatistics.maximum:
            return float(df.max())
        if self._statistics == GroupStatistics.negative_count:
            return int(df[df < 0].count())
        if self._statistics == GroupStatistics.positive_count:
            return int(df[df > 0].count())
        if self._statistics == GroupStatistics.unchanged:
            return int(df[df.abs() < 1e-3].count())
        if self._statistics == GroupStatistics.absmean:
            return float(df.abs().mean())
        raise NotImplementedError("Error: unknown operation requested!")

    def get_params(self):
        return {'centrality': self._centrality,
                'weight_inversion': self._weight_inversion,
                'statistics': self._statistics,
                'normalized': self._normalized,
                'weight': self._weight}

    def set_params(self, **kwargs):
        self._centrality = kwargs.get('centrality', self._centrality)
        self._weight_inversion = kwargs.get('weight_inversion', self._weight_inversion)
        self._normalized = kwargs.get('normalized', self._normalized)
        self._weight = kwargs.get('weight', self._weight)
        self._statistics = kwargs.get('statistics', self._statistics)
        self._check_status = False
        if self.is_ready:
            self._update_reference()
        return self
