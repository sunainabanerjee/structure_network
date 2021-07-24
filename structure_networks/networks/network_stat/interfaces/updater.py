from enum import Enum
from abc import ABC, abstractmethod

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__version__ = "1.0.0"
__all__ = ['NetworkUpdater', 'UpdaterTypes']


class UpdaterTypes(Enum):
    edge_attribute = 'edge_attribute'
    node_attribute = 'node_attribute'


class NetworkUpdater(ABC):
    @property
    @abstractmethod
    def type(self):
        pass

    @property
    @abstractmethod
    def attribute(self):
        pass

    @abstractmethod
    def __call__(self, g, **kwargs):
        pass
