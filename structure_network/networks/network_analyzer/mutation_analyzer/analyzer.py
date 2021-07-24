import numpy as np
import pandas as pd
from structure_network.structure import get_amino
from structure_network.structure import do_mutation
from scipy.optimize import linear_sum_assignment
from structure_network.networks.structure_network import CANetwork
from structure_network.networks.graph_algorithms import edge_inversion
from structure_network.networks.graph_algorithms import NodeCentrality
from structure_network.networks.graph_algorithms import EdgeCentrality
from structure_network.structure.amino_acids.mutation import Mutation
from structure_network.structure.amino_acids.mutation import parse_mutation
from structure_network.networks.graph_algorithms import edge_pair_ordering
from structure_network.networks.structure_network import ProteinStructureNetwork
from structure_network.networks.graph_algorithms import nodeflow_disruption_score
from structure_network.networks.graph_algorithms import edgeflow_disruption_score
from structure_network.networks.network_analyzer.structure_analyzer import StructureNetworkAnalyzer

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['MutationAnalyzer']


class MutationAnalyzer:
    def __init__(self, structure_network):
        assert isinstance(structure_network, ProteinStructureNetwork)
        self.__structure_network = StructureNetworkAnalyzer(structure_network)
        self.__mutations = []
        self.__network_analyzer = []

    def __len__(self):
        return len(self.__mutations)

    @property
    def residue_ids(self):
        return self.__structure_network.residue_ids

    @property
    def residue_names(self):
        return self.__structure_network.residue_names

    @property
    def structure(self):
        return self.__structure_network.structure

    @property
    def structure_network(self):
        return self.__structure_network.structure_network

    def mutant_structure_network(self, mutation):
        if isinstance(mutation, str):
            mutation = parse_mutation(mutation)
        if mutation in self.__mutations:
            return self.__network_analyzer[self.__mutations.index(mutation)].structure_network
        raise RuntimeError("Error: mutation {} is not registered".format(mutation))

    @property
    def g(self):
        return self.__structure_network.g

    def xyz(self, r=None):
        return self.__structure_network.xyz(r=r)

    @property
    def centroid(self):
        return self.__structure_network.centroid

    def register_mutation(self, r_id, r_name):
        if r_id not in self.residue_ids:
            return self
        mut_amino = get_amino(r_name)
        base_amino = self.__structure_network.structure.get_amino(r_id)
        if mut_amino != base_amino:
            mutation = Mutation(residue_id=r_id,
                                base_amino=base_amino,
                                mutated_amino=mut_amino)
            if mutation in self.__mutations:
                return self
            m_structure = do_mutation(self.__structure_network.structure,
                                      res_id=r_id,
                                      to_residue=r_name)
            structure_network = CANetwork(structure=m_structure,
                                          **self.__structure_network.params)
            self.__mutations.append(mutation)
            self.__network_analyzer.append(StructureNetworkAnalyzer(structure_network))
        return self

    @property
    def registered_mutations(self):
        return self.__mutations.copy()

    @property
    def mutation_positions(self):
        return list(set([m.residue_id for m in self.__mutations]))

    def node_centrality(self, **kwargs):
        wild = self.__structure_network.node_centrality(**kwargs.copy())
        centrality = kwargs.get('centrality', NodeCentrality.betweenness)
        return_difference = kwargs.pop('return_difference', False)
        wild_rname = self.__structure_network.residue_names
        wild_hdr = 'wild_{}'.format(centrality.name)
        scores = {wild_hdr: np.array([wild[n] for n in wild_rname])}

        for i, m in enumerate(self.__mutations):
            mut = self.__network_analyzer[i].node_centrality(**kwargs.copy())
            mut_rname = self.__network_analyzer[i].residue_names
            mut_hdr = '{}_{}'.format(m.name, centrality.name)
            scores[mut_hdr] = np.array([mut[n] for n in mut_rname])

        if return_difference:
            for m in scores:
                if m != wild_hdr:
                    scores[m] = scores[m] - scores[wild_hdr]
            del scores[wild_hdr]
        return pd.DataFrame(scores, index=self.residue_ids)

    def edge_centrality(self, **kwargs):
        centrality = kwargs.get('centrality', EdgeCentrality.betweenness)
        return_difference = kwargs.pop('return_difference', False)
        wild = self.__structure_network.edge_centrality(**kwargs.copy())

        edge_maps = [edge_pair_ordering(*edge) for edge in self.__structure_network.g.edges]
        edge_ids = [edge_pair_ordering(int(u[3:]), int(v[3:])) for u, v in edge_maps]
        wild_hdr = 'wild_{}'.format(centrality.name)
        scores = {wild_hdr: {}}
        for i, edge in enumerate(edge_maps):
            scores[wild_hdr][edge_ids[i]] = wild[edge]

        for i, m in enumerate(self.__mutations):
            mut = self.__network_analyzer[i].edge_centrality(**kwargs.copy())
            mut_edges = [edge_pair_ordering(*edge) for edge in self.__network_analyzer[i].g.edges]
            mut_edge_ids = [edge_pair_ordering(int(u[3:]), int(v[3:])) for u, v in mut_edges]
            mut_hdr = '{}_{}'.format(m.name, centrality.name)
            scores[mut_hdr] = {}
            for j, edge in enumerate(mut_edges):
                scores[mut_hdr][mut_edge_ids[j]] = mut[edge]
        scores = pd.DataFrame(scores)
        if return_difference:
            for m in scores.columns:
                if m != wild_hdr:
                    scores[m] = scores[m].values - scores[wild_hdr].values
            scores = scores.drop(columns=[wild_hdr])
        return scores

    def mcs(self, **kwargs):
        return_difference = kwargs.pop('return_difference', False)
        wild_hdr = 'base'
        clusters = {wild_hdr: self.__structure_network.mcs(**kwargs)}
        for i, mutant in enumerate(self.__mutations):
            clusters[mutant.name] = self.__network_analyzer[i].mcs(**kwargs)
        id_clusters = {}
        for entry in clusters:
            id_clusters[entry] = [[int(name[3:]) for name in group]
                                  for group in clusters[entry]]
        wild_clusters = id_clusters[wild_hdr]
        cluster_map = {wild_hdr: {i: cls for i, cls in enumerate(wild_clusters)} }
        for i, mutant in enumerate(self.__mutations):
            mut_clusters = id_clusters[mutant.name]
            cost = np.zeros((len(wild_clusters), len(mut_clusters)))
            for w_i, w_cls in enumerate(wild_clusters):
                for m_i, m_cls in enumerate(mut_clusters):
                    cost[w_i, m_i] = len(set(w_cls).intersection(m_cls))
            w_c, m_c = linear_sum_assignment(cost, maximize=True)
            cluster_map[mutant.name] = {w_i: mut_clusters[m_c[i]]
                                        for i, w_i in enumerate(w_c)}
        if return_difference:
            for mutant in self.__mutations:
                for index in cluster_map[wild_hdr]:
                    w_cls = cluster_map[wild_hdr][index]
                    m_cls = cluster_map[mutant.name].get(index, [])
                    cluster_map[mutant.name][index] = list(set(w_cls).symmetric_difference(m_cls))
            del cluster_map[wild_hdr]
        return pd.DataFrame(cluster_map)

    def node_disruption_score(self, **kwargs):
        centrality = kwargs.pop('centrality', NodeCentrality.betweenness)
        scaled = kwargs.pop('scaled', False)
        normalized = kwargs.pop('normalized', False)
        node_map = kwargs.pop('node_map', lambda x: int(x[3:]))
        weight = kwargs.pop('weight', 'weight')
        contribution = kwargs.get('contribution', False)
        weight_inversion = kwargs.get('weight_inversion', weight is not None)
        base_g = self.__structure_network.g
        if weight_inversion:
            base_g = edge_inversion(base_g, attribute=weight)
        scores, contribs = {}, {}
        for i, m in enumerate(self.__mutations):
            mut_g = self.__network_analyzer[i].g
            if weight_inversion:
                mut_g = edge_inversion(mut_g)
            result = nodeflow_disruption_score(base_g,
                                               mut_g,
                                               node_map=node_map,
                                               scaled=scaled,
                                               normalized=normalized,
                                               centrality=centrality,
                                               weight=weight,
                                               **kwargs)
            if contribution:
                score, contrib = result
                contribs[m] = contrib
            else:
                score = result
            scores[m] = score

        if contribution:
            return pd.DataFrame({'score': scores}), pd.DataFrame(contribs)
        return pd.DataFrame({'score': scores})

    def edge_disruption_score(self, **kwargs):
        centrality = kwargs.pop('centrality', EdgeCentrality.betweenness)
        scaled = kwargs.pop('scaled', False)
        normalized = kwargs.pop('normalized', False)
        node_map = kwargs.pop('node_map', lambda x: int(x[3:]))
        weight = kwargs.pop('weight', 'weight')
        contribution = kwargs.get('contribution', False)
        weight_inversion = kwargs.get('weight_inversion', weight is not None)
        base_g = self.__structure_network.g
        if weight_inversion:
            base_g = edge_inversion(base_g, attribute=weight)
        scores, contribs = {}, {}
        for i, m in enumerate(self.__mutations):
            mut_g = self.__network_analyzer[i].g
            if weight_inversion:
                mut_g = edge_inversion(mut_g)
            result = edgeflow_disruption_score(base_g,
                                               mut_g,
                                               centrality=centrality,
                                               node_map=node_map,
                                               scaled=scaled,
                                               normalized=normalized,
                                               weight=weight,
                                               **kwargs)
            if contribution:
                score, contrib = result
                contribs[m] = contrib
            else:
                score = result
            scores[m] = score

        if contribution:
            return pd.DataFrame({'score': scores}), pd.DataFrame(contribs)
        return pd.DataFrame({'score': scores})

