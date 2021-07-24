import numpy as np
import pandas as pd

__all__ = ['to_list',
           'order_disruption_score']


def to_list(x):
    """
    Convert an object to list, check string explicitly
    ensuring string does not convert to character list,
    any other object with iterator is explicitly convert
    to list. A None value converted to empty list.
    """
    if isinstance(x, str):
        return [x]
    if hasattr(x, '__iter__'):
        return list(x)
    return [x] if x is not None else []


def order_disruption_score(source_map, target_map, **kwargs):
    """
    Computes the disruption score between two key-value pair list.
    Uses value to rank the keys. In the common keys between two list
    rank disruption is scored.
    """
    node_map = kwargs.get('node_map', lambda x: x)
    scaled = kwargs.get('scaled', False)
    beta = kwargs.get('beta', 0)
    p = kwargs.get('p', 3)
    bias = kwargs.get('bias', 0.1)
    contribution = kwargs.get('contribution', False)

    assert isinstance(source_map, dict)
    assert isinstance(target_map, dict)
    assert (beta <= 1.) and (beta >= -1.)
    source_lookup = {node_map(k): v for k, v in source_map.items()}
    target_lookup = {node_map(k): v for k, v in target_map.items()}
    common_nodes = list(set([node_map(n)
                             for n in source_map]).intersection([node_map(n)
                                                                 for n in target_map]))
    if len(common_nodes) == 0:
        raise RuntimeError("Error: no common nodes between two score map!")
    source_rank = sorted(common_nodes,
                         key=lambda x: source_lookup[x])
    target_rank = sorted(common_nodes,
                         key=lambda x: target_lookup[x])

    df = {'source': {}, 'target': {}}
    for i, pair in enumerate(zip(source_rank, target_rank)):
        df['source'][pair[0]] = source_lookup[pair[0]] if scaled else i + 1
        df['target'][pair[1]] = target_lookup[pair[1]] if scaled else i + 1
    df = pd.DataFrame(df)
    n = len(common_nodes)
    df['weight'] = beta * (df['source'].values / n) + 0.5 * (1 - beta)
    df['score'] = df['weight'].values * (df['target'].values - df['source'].values) ** p
    score = np.mean(df['score'].values.reshape(-1))
    if contribution:
        if abs(score) < 1e-3:
            df['score'] = df['score'] + bias
        denominator = np.sum(df['score'].values.reshape(-1))
        contrib = df['score'].values/denominator
        contrib = {index: contrib[i] for i, index in enumerate(df.index)}
        return score, contrib
    return score
