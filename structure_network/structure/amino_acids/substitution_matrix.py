import os
import numpy as np
import pandas as pd
from .amino_acids import get_amino

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['SubstitutionMatrix']


class SubstitutionMatrix:
    def __init__(self, name, normalized=True):
        assert name in self.supported_models()
        data_file = os.path.join(os.path.dirname(__file__),
                                 'substitution', "{}.csv".format(name))
        df = pd.read_csv(data_file, header=0)
        self._name = name
        self._data = {}
        for i in df.index:
            aa_from = get_amino(df.loc[i, 'from'])
            aa_to = get_amino(df.loc[i, 'to'])
            value = df.loc[i, 'score']
            if aa_from not in self._data:
                self._data[aa_from] = {}
            self._data[aa_from][aa_to] = float(value)
        self._data = pd.DataFrame(self._data)
        if normalized:
            self._data = np.exp(self._data).div(np.exp(self._data).sum(axis=1), axis=0)

    @staticmethod
    def supported_models():
        data_dir = os.path.join(os.path.dirname(__file__), 'substitution')
        assert os.path.isdir(data_dir)
        return [f[:-4] for f in os.listdir(data_dir) if f.endswith('.csv')]

    @property
    def name(self):
        return self._name

    def __getitem__(self, item):
        from_aa, to_aa = get_amino(item[0]), get_amino(item[1])
        return self._data.loc[from_aa, to_aa]

    @property
    def from_aa(self):
        return list(self._data.index)

    @property
    def to_aa(self):
        return list(self._data.columns)

    @property
    def shape(self):
        return self._data.shape


