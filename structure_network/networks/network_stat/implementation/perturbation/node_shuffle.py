import networkx as nx
from structure_network.networks.core import to_list
from structure_network.networks.network_stat.interfaces import NetworkUpdater
from structure_network.networks.network_stat.interfaces import NetworkPerturbation
from structure_network.networks.network_stat.interfaces import NetworkGroupSelector

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NodeShuffle']


class NodeShuffle(NetworkPerturbation):
    def __init__(self, network_updaters, **kwargs):
        network_updaters = to_list(network_updaters)
        assert all([isinstance(updater, NetworkUpdater) for updater in network_updaters])
        self._updaters = to_list(network_updaters)
        self._selector = kwargs.get('selector', None)

    def register_updater(self, updater):
        assert isinstance(updater, NetworkUpdater)
        self._updaters.append(updater)
        return self

    def set_selector(self, selector):
        assert isinstance(selector, NetworkGroupSelector)
        assert selector.n_pair == 2
        self._selector = selector
        return self

    def perturb(self, g, **kwargs):
        node_selector = kwargs.get('node_selector', self._selector)
        node_attr = kwargs.get('node_attribute', 'name')
        k_select = kwargs.get('k', 1)
        assert isinstance(g, nx.Graph)
        assert isinstance(node_selector, NetworkGroupSelector)
        assert node_selector.n_pair == 2
        all_attr = nx.get_node_attributes(g, name=node_attr)
        if len(all_attr) == 0:
            raise RuntimeError("Error: node attribute {} "
                               "not defined".format(node_attr))
        nodes = set()
        g_new = g.copy(as_view=False)
        for _ in range(k_select):
            try:
                u, v = node_selector.next()
            except StopIteration:
                node_selector.reset()
                u, v = node_selector.next()
            if (v in all_attr) and (u in all_attr):
                u_attr = g_new.nodes[u][node_attr]
                v_attr = g_new.nodes[v][node_attr]
                g_new.nodes[v][node_attr] = u_attr
                g_new.nodes[u][node_attr] = v_attr
                nodes.add(u)
                nodes.add(v)
        kwargs['nodes'] = list(nodes)
        for updater in self._updaters:
            g_new = updater(g_new, nodes **kwargs)
        return g_new


