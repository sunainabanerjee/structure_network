import numpy as np
import pandas as pd
import networkx as nx
import concurrent.futures
from scipy.stats import kstest
from .utility import mutation_profile
from scipy.stats import gaussian_kde
from .network_metric import NetworkMetric
from ..interfaces import NetworkPerturbation

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NetworkMCMC']


class NetworkMCMC:
    def __init__(self, **kwargs):
        self._g = kwargs.get('g', None)
        self._metric = kwargs.get('metric', None)
        self._network_search = kwargs.get('search', None)

    def register_g(self, g):
        if isinstance(g, nx.Graph):
            self._g = g
        return self

    def register_metric(self, metric):
        if isinstance(metric, NetworkMetric):
            self._metric = metric
        return self

    def register_search_engine(self, search_engine):
        if isinstance(search_engine, NetworkPerturbation):
            self._network_search = search_engine
        return self

    def mcmc(self, **kwargs):
        g = kwargs.get('g', self._g)
        search_engine = kwargs.get('search_engine', self._network_search)
        metric = kwargs.get('metric', self._metric)
        batch_size = kwargs.get('batch_size', 20)
        max_iter = kwargs.get('max_iteration', 50000)
        ks_sample = kwargs.get('ks_sample', 100)
        ks_threshold = kwargs.get('ks_threshold', 0.05)
        p_threads = kwargs.get('p_threads', 3)
        min_sample = kwargs.get('min_sample', 10*batch_size)
        boot = kwargs.get('bootstrap', 10)

        # forcing single point mutation only
        k = min(kwargs.pop('k', 1), 1)
        kwargs['k'] = k

        assert min_sample > batch_size
        assert batch_size > 10
        assert ks_sample > batch_size

        if not isinstance(self._g, nx.Graph):
            raise RuntimeError("Error: no graph specified for search!")

        if not isinstance(search_engine, NetworkPerturbation):
            raise RuntimeError("Error: network perturbation engine for network "
                               "configuration search not registered!")

        def wrapper_function(perturb_fn, metric_fn, graph, **kw_args):
            g_new, mutation = perturb_fn(graph, **kw_args)
            m = metric_fn(g_new, **kw_args)
            return m, mutation

        converged, unstable = False, False
        samples = []
        mutations = []
        trend = []
        while (not converged) and (not unstable):
            if len(samples) > 0:
                print("Sample Size: {} and Average Value: {:.4f}".format(len(samples),
                                                                         np.mean(samples)))
            batch_id = 0
            while batch_id < batch_size:
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future_result = []
                    for _ in range(p_threads):
                        future_result.append(executor.submit(wrapper_function,
                                                             search_engine.perturb,
                                                             metric,
                                                             g,
                                                             **kwargs))
                    for f in future_result:
                        metric_value, mutation_seen = f.result()
                        samples.append(metric_value)
                        mutations.append(mutation_seen)
                        batch_id += 1
            if len(samples) > min_sample:
                ks_scores = []
                half_size = len(samples) // 2
                size = min(max(batch_size, ks_sample), half_size)
                for _ in range(boot):
                    x1 = np.random.choice(samples[:half_size], size)
                    x2 = np.random.choice(samples[half_size:], size)
                    ks_scores.append(kstest(x1, x2).statistic)
                trend.append(np.mean(ks_scores))
                converged = trend[-1] < ks_threshold
                if len(trend) > 10:
                    t_half = len(trend) // 2
                    unstable = (np.mean(trend[:t_half]) > np.mean(trend[t_half:])) or \
                               (len(samples) > max_iter)

        sampling_history = mutation_profile(mutations, samples)
        return converged, gaussian_kde(samples), sampling_history
