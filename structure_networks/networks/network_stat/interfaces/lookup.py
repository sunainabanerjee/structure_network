from abc import ABC, abstractmethod

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['Lookup', 'PairLookup']


class Lookup(ABC):
    @abstractmethod
    def set_network_context(self, g):
        pass

    @property
    @abstractmethod
    def is_ready(self):
        pass

    @property
    @abstractmethod
    def is_random(self):
        pass

    @property
    @abstractmethod
    def is_stateful(self):
        pass

    @abstractmethod
    def reset(self):
        pass

    @abstractmethod
    def __call__(self, x, **kwargs):
        pass


class PairLookup(ABC):
    @property
    @abstractmethod
    def is_ready(self):
        pass

    @abstractmethod
    def set_network_context(self, g):
        pass

    @property
    @abstractmethod
    def is_random(self):
        pass

    @property
    @abstractmethod
    def is_stateful(self):
        pass

    @abstractmethod
    def reset(self):
        pass

    @abstractmethod
    def __call__(self, x, y, **kwargs):
        pass

