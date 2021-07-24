from enum import Enum
from structure_network.networks.graph_algorithms import NodeCentrality
from structure_network.networks.graph_algorithms import EdgeCentrality

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['InducedGraphExtractor',
           'energy_compliant_weight_inversion']


class InducedGraphExtractor(Enum):
    shortest_path = 'shortest_path'
    connectivity = 'connectivity'
    percolation = 'percolation'


def energy_compliant_weight_inversion(centrality):
    if centrality == NodeCentrality.betweenness:
        return True
    if centrality == NodeCentrality.load:
        return True
    if centrality == NodeCentrality.information:
        return True
    if centrality == NodeCentrality.communicability:
        return True
    if centrality == NodeCentrality.percolation:
        return True
    if centrality == NodeCentrality.eigen:
        return False
    if centrality == NodeCentrality.current_flow:
        return True
    if centrality == NodeCentrality.current_flow_closeness:
        return True
    if centrality == NodeCentrality.harmonic:
        return True
    if centrality == NodeCentrality.subgraph:
        return True
    if centrality == NodeCentrality.degree:
        return False
    if centrality == NodeCentrality.closeness:
        return True
    if centrality == NodeCentrality.katz:
        return False
    if centrality == NodeCentrality.second_order:
        return True

    if centrality == EdgeCentrality.betweenness:
        return True
    if centrality == EdgeCentrality.load:
        return True
    if centrality == EdgeCentrality.current_flow:
        return True
    raise RuntimeError("Error: unknown centrality {}".format(centrality))
