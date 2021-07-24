from abc import ABC, abstractmethod

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NetworkGroupSelector']


class NetworkGroupSelector(ABC):
    @property
    @abstractmethod
    def n_pair(self):
        pass

    def __iter__(self):
        return self.reset()

    def __next__(self):
        return self.next()

    @abstractmethod
    def next(self):
        pass

    @abstractmethod
    def reset(self):
        pass
