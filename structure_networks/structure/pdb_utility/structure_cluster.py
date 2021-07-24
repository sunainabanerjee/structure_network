import logging
import numpy as np
from scipy.cluster.hierarchy import linkage
from structure_networks.structure import CaTrace
from .align import catraces_rmsd

__version__ = "1.0"
__all__ = ['cluster_catraces']


def cluster_catraces(lst, method='single'):
    logger = logging.getLogger('structural_dynamics.cluster_catraces')
    assert isinstance(lst, list)
    trace_len = None
    for item in lst:
        assert isinstance(item, CaTrace)
        if trace_len is None:
            trace_len = len(item)
        assert trace_len == len(item)
    logger.debug("Length of each calpha trace is (%d)" % trace_len)
    n = len(lst)
    logger.debug("Number of snapshots is (%d)" % n)
    scores = np.zeros((n, n))
    for i in range(n-1):
        logger.debug("Evaluating distance from structure (%d)" % i)
        for j in range(i+1,n):
            s = catraces_rmsd(lst[i], lst[j],return_structure=False)
            scores[i, j] = s
            scores[j, i] = s
    linkage_matrix = linkage(scores, method=method)
    return linkage_matrix
