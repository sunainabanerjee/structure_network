from enum import Enum
from abc import ABC, abstractmethod

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NetworkMetric', 'GroupStatistics']


class GroupStatistics(Enum):
    average = 'mean'
    maximum = 'max'
    median = 'median'
    absmean = 'absolute_mean'
    positive_count = 'count_positives'
    negative_count = 'count_negatives'
    unchanged = 'unchanged'


class NetworkMetric(ABC):
    @abstractmethod
    def set_base_network(self, g):
        pass

    @property
    @abstractmethod
    def g(self):
        pass

    @property
    @abstractmethod
    def name(self):
        pass

    @property
    @abstractmethod
    def target_group(self):
        pass

    @abstractmethod
    def set_target(self, group):
        pass

    @property
    @abstractmethod
    def is_ready(self):
        pass

    @abstractmethod
    def __call__(self, g, **kwargs):
        pass

    @abstractmethod
    def get_params(self):
        pass

    @abstractmethod
    def set_params(self, **kwargs):
        pass
