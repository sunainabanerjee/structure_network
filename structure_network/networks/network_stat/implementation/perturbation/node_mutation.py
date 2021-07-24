import networkx as nx
from structure_network.networks.core import to_list
from structure_network.structure.amino_acids import Mutation
from structure_network.networks.network_stat.interfaces import Lookup
from structure_network.networks.network_stat.interfaces import NetworkGroupSelector
from structure_network.networks.network_stat.interfaces import NetworkUpdater
from structure_network.networks.network_stat.interfaces import NetworkPerturbation

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NodeMutation']


class NodeMutation(NetworkPerturbation):
    def __init__(self, **kwargs):
        self._updaters = to_list(kwargs.get('network_updaters', []))
        self._lookup = kwargs.get('lookup', None)
        self._selector = kwargs.get('node_selector', None)
        assert all([isinstance(updater, NetworkUpdater) for updater in self._updaters])

    def set_lookup(self, lookup):
        if isinstance(lookup, Lookup):
            self._lookup = lookup
        return self

    def set_selector(self, selector):
        if isinstance(selector, NetworkGroupSelector) and (selector.n_pair == 1):
            self._selector = selector
        return self

    def register_updater(self, updater):
        assert isinstance(updater, NetworkUpdater)
        self._updaters.append(updater)
        return self

    def perturb(self, g, **kwargs):
        node_selector = kwargs.get('node_selector', self._selector)
        k_select = kwargs.get('k', 1)
        assert isinstance(node_selector, NetworkGroupSelector)
        assert isinstance(g, nx.Graph)
        assert node_selector.n_pair == 1
        assert isinstance(self._lookup, Lookup)
        node_attr = kwargs.get('node_attribute', 'name')
        g_new = g.copy(as_view=False)
        nodes = []
        mutation_list = []
        self._lookup.set_network_context(g)
        for _ in range(k_select):
            try:
                u = node_selector.next()
            except StopIteration:
                node_selector.reset()
                u = node_selector.next()
            g_new.nodes[u][node_attr] = self._lookup(u, id=u)
            mutation_list.append(Mutation(residue_id=u,
                                          base_amino=g.nodes[u][node_attr],
                                          mutated_amino=g_new.nodes[u][node_attr]))
            nodes.append(u)
        kwargs['nodes'] = list(nodes)
        for updater in self._updaters:
            g_new = updater(g_new, **kwargs)
        return (g_new, mutation_list[0]) if k_select == 1 else (g_new, mutation_list)

