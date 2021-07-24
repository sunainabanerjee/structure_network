import threading
import numpy as np
import networkx as nx
from structure_networks.structure import get_amino
from structure_networks.structure import valid_amino_acids
from structure_networks.structure import SubstitutionMatrix
from structure_networks.networks.network_stat.interfaces import Lookup

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['MutationSaturationError',
           'AminoSubstitutionLookup']


class MutationSaturationError(Exception):
    def __init__(self, *args):
        super(MutationSaturationError, self).__init__(*args)


class AminoSubstitutionLookup(Lookup):
    def __init__(self, **kwargs):
        self._g = kwargs.get('g', None)
        self._lookup_name = kwargs.get('node_attr', 'name')
        self._selection_bias = kwargs.get('selection_bias', 'UNIFORM')
        self._stateful = kwargs.get('stateful', False)
        self._states = {}
        self._aminos = [get_amino(aa) for aa in valid_amino_acids()]
        self._selection_bias = SubstitutionMatrix(name=self._selection_bias,
                                                  normalized=True)
        self._lock = threading.Lock()
        self._check_ready = False
        self._check_status = False # Caching the result

    @property
    def is_ready(self):
        if not self._check_status:
            self._check_ready = isinstance(self._g, nx.Graph) and \
                                (len(nx.get_node_attributes(self._g, self._lookup_name)) == self._g.order())
            self._check_status = True
        return self._check_ready

    def set_network_context(self, g):
        self._lock.acquire()
        try:
            if isinstance(g, nx.Graph):
                self._g = g
                self._check_status = False
        finally:
            self._lock.release()
        return self

    @property
    def selection_bias(self):
        return self._selection_bias.name

    @property
    def is_random(self):
        return True

    @property
    def is_stateful(self):
        return self._stateful

    def reset(self):
        self._lock.acquire()
        try:
            self._states.clear()
        finally:
            self._lock.release()
        return self

    def __call__(self, x, **kwargs):
        if not self.is_ready:
            raise RuntimeError("Error: network context is not setup!")
        assert x in self._g.nodes
        state_id = None
        if self._stateful:
            self._lock.acquire()
        try:
            if self._stateful:
                state_id = kwargs.get('id', x)
                if state_id not in self._states:
                    self._states[state_id] = set()
            aa_base = get_amino(self._g.nodes[x][self._lookup_name])
            amino_index = self._aminos.index(aa_base)
            select_index = amino_index
            p = np.array([self._selection_bias[aa_base, aa] for aa in self._aminos])
            p[amino_index] = 0.
            select_continue = True
            while select_continue and (np.sum(p) > 0):
                p = p / np.sum(p)
                select_index = np.random.choice(len(self._aminos), p=p)
                select_continue = (select_index == amino_index)
                if self._stateful:
                    select_continue = select_continue and \
                                      (select_index not in self._states[state_id])
                if select_continue:
                    p[select_index] = 0.
            if select_continue:
                raise MutationSaturationError("Error: no available mutation "
                                              "at state {}".format(state_id))
            if self._stateful:
                self._states[state_id].add(select_index)
        finally:
            if self._stateful:
                self._lock.release()
        return self._aminos[select_index].name()
