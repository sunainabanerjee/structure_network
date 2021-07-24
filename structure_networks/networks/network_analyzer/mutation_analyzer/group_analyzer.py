import pandas as pd
import multiprocessing as mp
from structure_networks.structure import get_amino
from structure_networks.structure import do_mutation
from structure_networks.structure.amino_acids import Mutation
from structure_networks.networks.structure_network import CANetwork
from structure_networks.networks.core import order_disruption_score
from structure_networks.networks.graph_algorithms import NodeCentrality
from structure_networks.networks.graph_algorithms import EdgeCentrality
from structure_networks.networks.graph_algorithms import edge_inversion
from structure_networks.structure.amino_acids.mutation import parse_mutation
from structure_networks.networks.graph_algorithms import nodeflow_disruption_score
from structure_networks.networks.graph_algorithms import edgeflow_disruption_score
from structure_networks.networks.graph_algorithms import validate_node_centrality
from structure_networks.networks.network_analyzer.structure_analyzer import SubstructureNetworkAnalyzer

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['MutationGroupAnalyzer']


def f_nds(mutant,
          base_g,
          wild_g,
          centrality,
          scaled,
          normalized,
          weight,
          node_map,
          kwargs):
    score = nodeflow_disruption_score(base_g,
                                      wild_g,
                                      centrality=centrality,
                                      scaled=scaled,
                                      normalized=normalized,
                                      weight=weight,
                                      node_map=node_map,
                                      **kwargs)
    return mutant, score


def f_ods(mutant,
          base_centrality,
          wild_centrality,
          scaled,
          node_map,
          kwargs):
    if (not validate_node_centrality(base_centrality)) or \
            (not validate_node_centrality(wild_centrality)):
        return mutant, 0.
    score = order_disruption_score(base_centrality,
                                   wild_centrality,
                                   scaled=scaled,
                                   node_map=node_map,
                                   **kwargs)
    return mutant, score


def node_mapper(x):
    return int(x[3:])


class MutationGroupAnalyzer:
    def __init__(self, structure_network, **kwargs):
        self.__kwargs = kwargs.copy()
        self.__group_analyzer = SubstructureNetworkAnalyzer(structure_network, **kwargs)
        self.__mutations = []
        self.__mutant_group_analyzer = []

    def register_group(self, group, **kwargs):
        self.__group_analyzer.register_group(group, **kwargs)
        for analyzer in self.__mutant_group_analyzer:
            analyzer.register_group(group, **kwargs)
        return self

    @property
    def registered_groups(self):
        return self.__group_analyzer.registered_groups

    @property
    def registered_mutations(self):
        return self.__mutations

    def __getitem__(self, item):
        return self.__group_analyzer[item]

    @property
    def cover(self):
        return self.__group_analyzer.cover

    @property
    def residue_ids(self):
        return self.__group_analyzer.residue_ids

    @property
    def residue_names(self):
        return self.__group_analyzer.residue_names

    def reset_cache(self):
        self.__group_analyzer.reset_cache()
        for r in self.__mutant_group_analyzer:
            r.reset_cache()
        return self

    @property
    def structure(self):
        return self.__group_analyzer.structure

    @property
    def g(self):
        return self.__group_analyzer.g

    @property
    def params(self):
        return self.__group_analyzer.params

    def xyz(self, r=None):
        return self.__group_analyzer.xyz(r=r)

    def has_covered(self, r):
        return self.__group_analyzer.has_covered(r)

    def group_contains(self, r):
        return self.__group_analyzer.group_contains(r)

    def common_residues(self, *args, **kwargs):
        mutant = kwargs.get('mutant', None)
        if mutant is None:
            return self.__group_analyzer.common_residues(*args)
        if isinstance(mutant, str):
            mutant = parse_mutation(mutant)
        if mutant not in self.__mutations:
            raise RuntimeError("Error: mutation {} not registered!".format(mutant))
        index = self.__mutations.index(mutant)
        return self.__mutant_group_analyzer[index].common_residues(*args)

    def register_mutation(self, r_id, r_name):
        assert r_id in self.residue_ids
        mut_amino = get_amino(r_name)
        base_amino = self.__group_analyzer.structure.get_amino(r_id)
        if mut_amino != base_amino:
            mutation = Mutation(residue_id=r_id,
                                base_amino=base_amino,
                                mutated_amino=mut_amino)
            if mutation in self.__mutations:
                return self
            m_structure = do_mutation(self.__group_analyzer.structure,
                                      res_id=r_id,
                                      to_residue=r_name)
            structure_network = CANetwork(structure=m_structure,
                                          **self.__group_analyzer.params)
            self.__mutations.append(mutation)
            analyzer = SubstructureNetworkAnalyzer(structure_network, **self.__kwargs.copy())
            for group in self.registered_groups:
                analyzer.register_group(group=self.__group_analyzer[group], name=group)
            self.__mutant_group_analyzer.append(analyzer)
        return self

    def __len__(self):
        return len(self.__mutant_group_analyzer)

    @property
    def n_groups(self):
        return len(self.__group_analyzer)

    def subgraph(self, g1, g2, mutant=None):
        if mutant is None:
            return self.__group_analyzer.subgraph(g1, g2)
        if isinstance(mutant, str):
            mutant = parse_mutation(mutant)
        if mutant not in self.__mutations:
            raise RuntimeError("Error: mutation {} not registered".format(mutant))
        return self.__mutant_group_analyzer[self.__mutations.index(mutant)].subgraph(g1, g2)

    def mutant_g(self, mutant):
        if isinstance(mutant, str):
            mutant = parse_mutation(mutant)
        if mutant not in self.__mutations:
            raise RuntimeError("Error: mutation {} not registered".format(mutant))
        return self.__mutant_group_analyzer[self.__mutations.index(mutant)].g

    def node_centrality(self, **kwargs):
        return_difference = kwargs.pop('return_difference', False)
        centrality = {'base': self.__group_analyzer.node_centrality(**kwargs.copy())}

        for i, mutant in enumerate(self.__mutations):
            centrality[mutant.name] = self.__mutant_group_analyzer[i].node_centrality(**kwargs.copy())

        if return_difference:
            centrality_difference = {}
            for mutant in self.__mutations:
                data = {}
                df = centrality[mutant.name]
                for i in self.registered_groups:
                    data[i] = {}
                    for j in self.registered_groups:
                        if i == j:
                            continue
                        data[i][j] = {}
                        mut = df[i][j]
                        wild = centrality['base'][i][j]
                        for k in set(mut.keys()).intersection(wild.keys()):
                            data[i][j][k] = mut[k] - wild[k]
                centrality_difference[mutant.name] = data.copy()
            centrality = centrality_difference
        return centrality

    def edge_centrality(self, **kwargs):
        return_difference = kwargs.pop('return_difference', False)
        centrality = {'base': self.__group_analyzer.edge_centrality(**kwargs.copy())}
        for i, mutant in enumerate(self.__mutations):
            centrality[mutant.name] = self.__mutant_group_analyzer[i].edge_centrality(**kwargs.copy())
        if return_difference:
            centrality_difference = {}
            for mutant in self.__mutations:
                data = {}
                df = centrality[mutant.name]
                for i in self.registered_groups:
                    data[i] = {}
                    for j in self.registered_groups:
                        if i == j:
                            continue
                        data[i][j] = {}
                        mut = df[i][j]
                        wild = centrality['base'][i][j]
                        for k in set(mut.keys()).intersection(wild.keys()):
                            data[i][j][k] = mut[k] - wild[k]
                centrality_difference[mutant.name] = data.copy()
            centrality = centrality_difference
        return centrality

    def node_disruption_score(self, group1, group2, **kwargs):
        assert group1 in self.registered_groups
        assert group2 in self.registered_groups
        assert group1 != group2
        centrality = kwargs.pop('centrality', NodeCentrality.betweenness)
        scaled = kwargs.pop('scaled', False)
        normalized = kwargs.pop('normalized', False)
        node_map = kwargs.pop('node_map', node_mapper)
        weight = kwargs.pop('weight', 'weight')
        use_cache = kwargs.pop('use_cache', False)
        contribution = kwargs.get('contribution', False)
        p = kwargs.pop('p', 3)
        inverse_weight = kwargs.pop('inverse_weight',
                                    self.__group_analyzer.weight_inversion)

        scores, contributions = {}, {}
        pool = mp.Pool(mp.cpu_count())
        result_obj = []
        kwargs['p'] = p

        if use_cache:
            base_centrality = self.__group_analyzer.group_node_centrality(group1,
                                                                          group2,
                                                                          centrality=centrality,
                                                                          normalized=normalized,
                                                                          weight=weight,
                                                                          inverse_weight=inverse_weight,
                                                                          **kwargs)
            for i, m in enumerate(self.__mutations):
                wild_centrality = self.__mutant_group_analyzer[i].group_node_centrality(group1,
                                                                                        group2,
                                                                                        centrality=centrality,
                                                                                        normalized=normalized,
                                                                                        weight=weight,
                                                                                        inverse_weight=inverse_weight,
                                                                                        **kwargs)
                result_obj.append(pool.apply_async(f_ods,
                                                   args=(m,
                                                         base_centrality.copy(),
                                                         wild_centrality,
                                                         scaled,
                                                         node_map,
                                                         kwargs)))
        else:
            base_g = self.__group_analyzer.subgraph(group1, group2)
            if inverse_weight != self.__group_analyzer.weight_inversion:
                mx = self.__group_analyzer.structure_network.max_edge_weight
                mn = self.__group_analyzer.structure_network.min_edge_weight
                base_g = edge_inversion(base_g,
                                        attribute=weight,
                                        max_value=mx,
                                        min_value=mn)

            for i, m in enumerate(self.__mutations):
                wild_g = self.__mutant_group_analyzer[i].subgraph(group1, group2)
                if inverse_weight != self.__mutant_group_analyzer[i].weight_inversion:
                    mx = self.__mutant_group_analyzer[i].structure_network.max_edge_weight
                    mn = self.__mutant_group_analyzer[i].structure_network.min_edge_weight
                    wild_g = edge_inversion(wild_g, attribute=weight, max_value=mx, mn=mn)

                result_obj.append(pool.apply_async(f_nds,
                                                   args=(m,
                                                         base_g,
                                                         wild_g,
                                                         centrality,
                                                         scaled,
                                                         normalized,
                                                         weight,
                                                         node_map,
                                                         kwargs)))
        results = [r.get() for r in result_obj]
        pool.close()
        pool.join()

        for m, score in results:
            if contribution:
                score, contrib = score
                contributions[m] = contrib
            scores[m] = score

        if contribution:
            return pd.DataFrame({'score': scores}), pd.DataFrame(contributions)
        return pd.DataFrame({'score': scores})

    def edge_disruption_score(self, group1, group2, **kwargs):
        assert group1 in self.registered_groups
        assert group2 in self.registered_groups
        assert group1 != group2
        centrality = kwargs.pop('centrality', EdgeCentrality.betweenness)
        scaled = kwargs.pop('scaled', False)
        normalized = kwargs.pop('normalized', False)
        node_map = kwargs.pop('node_map', lambda x: int(x[3:]))
        weight = kwargs.pop('weight', 'weight')
        contribution = kwargs.get('contribution', False)
        inverse_weight = kwargs.pop('inverse_weight',
                                    self.__group_analyzer.weight_inversion)

        base_g = self.__group_analyzer.subgraph(group1, group2)

        if inverse_weight != self.__group_analyzer.weight_inversion:
            mx = self.__group_analyzer.structure_network.max_edge_weight
            mn = self.__group_analyzer.structure_network.min_edge_weight
            base_g = edge_inversion(base_g,
                                    attribute=weight,
                                    max_value=mx,
                                    min_value=mn)

        scores, contributions = {}, {}
        for i, m in enumerate(self.__mutations):
            wild_g = self.__mutant_group_analyzer[i].subgraph(group1, group2)
            if inverse_weight != self.__mutant_group_analyzer[i].weight_inversion:
                mx = self.__mutant_group_analyzer[i].structure_network.max_edge_weight
                mn = self.__mutant_group_analyzer[i].structure_network.min_edge_weight
                wild_g = edge_inversion(wild_g,
                                        attribute=weight,
                                        max_value=mx,
                                        min_value=mn)

            score = edgeflow_disruption_score(base_g,
                                              wild_g,
                                              centrality=centrality,
                                              scaled=scaled,
                                              normalized=normalized,
                                              weight=weight,
                                              node_map=node_map,
                                              **kwargs)
            if contribution:
                score, contrib = score
                contributions[m] = contrib
            scores[m] = score
        if contribution:
            return pd.DataFrame({'score': scores}), pd.DataFrame(contributions)
        return pd.DataFrame({'score': scores})



