from abc import ABC, abstractmethod

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['NetworkPerturbation']


class NetworkPerturbation(ABC):
    @abstractmethod
    def set_selector(self, selector):
        pass

    @abstractmethod
    def register_updater(self, updater):
        pass

    @abstractmethod
    def perturb(self, g, **kwargs):
        pass

