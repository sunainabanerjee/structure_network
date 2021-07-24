import numpy as np
import pandas as pd
import networkx as nx
from structure_networks.structure import CaTrace
from structure_networks.structure import PDBStructure
from ..graph_algorithms import is_connected
from ..graph_algorithms import NodeCentrality
from ..graph_algorithms import node_centrality
from structure_networks.networks.structure_network import CANetwork
from structure_networks.structure import valid_amino_acids, pdb_to_catrace
from structure_networks.networks.structure_network import ProteinStructureNetwork
from structure_networks.networks.structure_network import ContactTypes, ThresholdTypes

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['distance_cutoff_for_connectivity',
           'strongly_connected_component',
           'mutational_centrality_change']


def distance_cutoff_for_connectivity(structure, **kwargs):
    min_cutoff = kwargs.get('min_distance_cutoff', 2.5)
    increment = kwargs.get('distance_increment', 0.1)
    contact_type = kwargs.get('contact_type', ContactTypes.all_contacts)
    assert (increment > 0) and (min_cutoff > 0)
    assert isinstance(structure, (PDBStructure, CaTrace))
    distance_cutoff = float(min_cutoff)
    pdb_network = CANetwork(structure,
                            distance_threshold=distance_cutoff,
                            threshold_type=ThresholdTypes.distance,
                            contact_type=contact_type)
    while not is_connected(pdb_network.g):
        distance_cutoff = distance_cutoff + increment
        pdb_network.set_distance_threshold(d=distance_cutoff)
    return distance_cutoff


def strongly_connected_component(pdb_network):
    assert isinstance(pdb_network, ProteinStructureNetwork)
    return list(nx.kosaraju_strongly_connected_components(pdb_network.g.to_directed()))


def mutational_centrality_change(structure, residue, **kwargs):
    if isinstance(structure, PDBStructure):
        structure = pdb_to_catrace(structure)
    assert isinstance(structure, CaTrace)
    assert residue in structure.residue_ids
    centrality = kwargs.pop('centrality', NodeCentrality.betweenness)
    contact_type = kwargs.pop('contact_type', ContactTypes.all_contacts)
    distance_threshold = kwargs.pop('distance_threshold', distance_cutoff_for_connectivity(structure))
    energy_threshold = kwargs.pop('energy_threshold', 0.)

    base_amino = structure.get_amino(resid=residue).name(one_letter_code=True)
    all_amino = valid_amino_acids(one_letter=True)
    all_residue = structure.residue_ids
    centrality_measure = np.zeros((len(all_amino), len(all_residue)))

    for i, aa in enumerate(all_amino):
        structure.set_amino(resid=residue, aa_type=aa)
        pn = CANetwork(structure,
                       distance_threshold=distance_threshold,
                       energy_threshold=energy_threshold,
                       threshold_type=ThresholdTypes.distance_energy,
                       contact_type=contact_type)
        scores = node_centrality(pn.g, centrality=centrality)
        centrality_measure[i, :] = np.array([scores[r] if r in scores else 0 for r in all_residue])
    base_idx = all_amino.index(base_amino)
    centrality_measure = centrality_measure - centrality_measure[base_idx, :]
    return pd.DataFrame(centrality_measure, columns=all_residue, index=all_amino)

