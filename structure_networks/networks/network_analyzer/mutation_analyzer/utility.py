import numpy as np
from .analyzer import MutationAnalyzer
from .group_analyzer import MutationGroupAnalyzer

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['node_centrality_changes',
           'edge_centrality_changes']


def node_centrality_changes(analyzer, **kwargs):
    assert isinstance(analyzer, (MutationAnalyzer, MutationGroupAnalyzer)), \
        "Error: unsupported analyzer instance {}".format(type(analyzer))
    if isinstance(analyzer, MutationAnalyzer):
        return node_centrality_changes_mutants(analyzer, **kwargs)
    else:
        return node_centrality_changes_group(analyzer, **kwargs)


def edge_centrality_changes(analyzer, **kwargs):
    assert isinstance(analyzer, (MutationAnalyzer, MutationGroupAnalyzer)), \
        "Error: unsupported analyzer instance {}".format(type(analyzer))
    if isinstance(analyzer, MutationAnalyzer):
        return edge_centrality_changes_mutants(analyzer, **kwargs)
    else:
        return edge_centrality_changes_group(analyzer, **kwargs)


def node_centrality_changes_group(analyzer, **kwargs):
    assert isinstance(analyzer, MutationGroupAnalyzer), \
        "Error: not a group mutation analyzer instance!"
    assert (len(analyzer) > 0) or (analyzer.n_groups == 0), \
        "Error: no mutation or group registered for analysis!"
    if 'return_difference' in kwargs:
        kwargs.pop('return_difference')
    q = min(kwargs.pop('q', 0.9), 0.99)
    report_all = kwargs.get('report_all', True)
    if q < 0.5:
        q = 1. - q
    result = analyzer.node_centrality(return_difference=True, **kwargs)
    processed_result = {}
    for mutant in result:
        data = result[mutant]
        values = []
        for i in data:
            for j in data[i]:
                if data[i][j] is not np.nan:
                    for k, v in data[i][j].items():
                        if abs(v) > 0:
                            values.append(v)
        values = np.array(values)
        higher_cutoff, lower_cutoff = None, None
        positive_index = np.where(values > 0)[0]
        if len(positive_index) > 10:
            higher_cutoff = np.quantile(values[positive_index], q=q)
        negative_index = np.where(values < 0)[0]
        if len(negative_index) > 10:
            lower_cutoff = np.quantile(values[negative_index], q=1-q)

        if report_all:
            higher_cutoff = 0
            lower_cutoff = 0
        screened_data = {}

        for i in data:
            screened_data[i] = {}
            for j in data[i]:
                screened_data[i][j] = {}
                if data[i][j] is not np.nan:
                    if higher_cutoff is not None:
                        screened_data[i][j]['positive'] = sorted([k for k, v in data[i][j].items()
                                                                  if v > higher_cutoff],
                                                                 key=lambda x: data[i][j][x],
                                                                 reverse=True)
                    else:
                        screened_data[i][j]['positive'] = []
                    if lower_cutoff is not None:
                        screened_data[i][j]['negative'] = sorted([k for k, v in data[i][j].items()
                                                                  if v < lower_cutoff],
                                                                 key=lambda x: data[i][j][x],
                                                                 reverse=False)
                    else:
                        screened_data[i][j]['negative'] = []
        processed_result[mutant] = screened_data
    return processed_result


def node_centrality_changes_mutants(analyzer, **kwargs):
    assert isinstance(analyzer, MutationAnalyzer), "Error: not a mutational analyzer instance!"
    assert len(analyzer) > 0, "Error: no mutation configured for analysis!"
    if 'return_difference' in kwargs:
        kwargs.pop('return_difference')
    q = min(kwargs.pop('q', 0.9), 0.99)
    report_all = kwargs.get('report_all', False)
    if q < 0.5:
        q = 1. - q
    result = analyzer.node_centrality(return_difference=True, **kwargs)
    differences = result.values.reshape(-1)
    residue_indices = np.asarray(result.index)
    effects = {c: {} for c in result.columns}

    index = np.where(differences > 0)[0]
    if len(index) > (0 if report_all else 10):
        upper_cutoff = np.quantile(differences[index], q=q) if not report_all else 0
        for c in effects:
            data = np.asarray(result[c].values)
            pos_index = np.where(data > upper_cutoff)[0]
            if len(pos_index) == 0:
                effects[c]['positive'] = []
            else:
                effects[c]['positive'] = sorted(residue_indices[pos_index].tolist(),
                                                key=lambda x: result.loc[x, c],
                                                reverse=True)

    index = np.where(differences < 0)[0]
    if len(index) > (0 if report_all else 10):
        lower_cutoff = np.quantile(differences[index], q=(1-q)) if not report_all else 0
        for c in effects:
            data = np.asarray(result[c].values)
            neg_index = np.where(data < lower_cutoff)[0]
            if len(neg_index) > 0:
                effects[c]['negative'] = sorted(residue_indices[neg_index].tolist(),
                                                key=lambda x: result.loc[x, c],
                                                reverse=False)
            else:
                effects[c]['negative'] = []
    return effects


def edge_centrality_changes_mutants(analyzer, **kwargs):
    assert isinstance(analyzer, MutationAnalyzer), "Error: not a mutational analyzer instance!"
    assert len(analyzer) > 0, "Error: no mutation configured!"
    if 'return_difference' in kwargs:
        kwargs.pop('return_difference')
    report_all = kwargs.get('report_all', False)
    q = min(kwargs.pop('q', 0.9), 0.99)
    if q < 0.5:
        q = 1. - q
    result = analyzer.edge_centrality(return_difference=True, **kwargs)
    differences = result.values.reshape(-1)
    differences = np.where(np.isnan(differences), 0, differences)
    edge_indices = np.asarray(result.index)
    effects = {c: {} for c in result.columns}

    index = np.where(differences > 0)[0]
    if len(index) > (0 if report_all else 10):
        upper_cutoff = np.quantile(differences[index], q=q) if not report_all else 0
        for c in effects:
            data = np.asarray(result[c].values)
            data = np.where(np.isnan(data), 0, data)
            pos_index = np.where(data > upper_cutoff)[0]
            if len(pos_index) > 0:
                effects[c]['positive'] = sorted(edge_indices[pos_index].tolist(),
                                                key=lambda x: result.loc[x, c],
                                                reverse=True)
            else:
                effects[c]['positive'] = []

    index = np.where(differences < 0)[0]
    if len(index) > (0 if report_all else 10):
        lower_cutoff = np.quantile(differences[index], q=(1-q)) if not report_all else 0
        for c in effects:
            data = np.asarray(result[c].values)
            data = np.where(np.isnan(data), 0, data)
            neg_index = np.where(data < lower_cutoff)[0]
            if len(neg_index) > 0:
                effects[c]['negative'] = sorted(edge_indices[neg_index].tolist(),
                                                key=lambda x: result.loc[x, c])
            else:
                effects[c]['negative'] = []
    return effects


def edge_centrality_changes_group(analyzer, **kwargs):
    assert isinstance(analyzer, MutationGroupAnalyzer), \
        "Error: not a group mutation analyzer instance!"
    assert (len(analyzer) > 0) or (analyzer.n_groups == 0), \
        "Error: no mutation or group registered for analysis!"
    if 'return_difference' in kwargs:
        kwargs.pop('return_difference')
    report_all = kwargs.get('report_all', True)
    q = min(kwargs.pop('q', 0.9), 0.99)
    if q < 0.5:
        q = 1. - q
    result = analyzer.edge_centrality(return_difference=True, **kwargs)
    processed_result = {}
    for mutant in result:
        data = result[mutant]
        values = []
        for i in data:
            for j in data[i]:
                if data[i][j] is not np.nan:
                    for k, v in data[i][j].items():
                        if abs(v) > 0:
                            values.append(v)
        values = np.array(values)
        higher_cutoff, lower_cutoff = None, None
        positive_index = np.where(values > 0)[0]
        if len(positive_index) > 10:
            higher_cutoff = np.quantile(values[positive_index], q=q)
        if report_all:
            higher_cutoff = 0

        negative_index = np.where(values < 0)[0]
        if len(negative_index) > 10:
            lower_cutoff = np.quantile(values[negative_index], q=1-q)
        if report_all:
            lower_cutoff = 0
        screened_data = {}
        for i in data:
            screened_data[i] = {}
            for j in data[i]:
                screened_data[i][j] = {}
                if data[i][j] is not np.nan:
                    if higher_cutoff is not None:
                        filtered_data = [k for k, v in data[i][j].items() if v > higher_cutoff]
                        screened_data[i][j]['positive'] = sorted(filtered_data,
                                                                 key=lambda x: data[i][j][x],
                                                                 reverse=True)
                    else:
                        screened_data[i][j]['positive'] = []
                    if lower_cutoff is not None:
                        filtered_data = [k for k, v in data[i][j].items() if v < lower_cutoff]
                        screened_data[i][j]['negatives'] = sorted(filtered_data,
                                                                  key=lambda x: data[i][j][x],
                                                                  reverse=False)
                    else:
                        screened_data[i][j]['negative'] = []
        processed_result[mutant] = screened_data
    return processed_result

