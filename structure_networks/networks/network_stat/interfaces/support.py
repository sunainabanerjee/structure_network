from abc import ABC,abstractmethod

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['SupportFunction']


class SupportFunction(ABC):
    @property
    @abstractmethod
    def attribute(self):
        pass

    @property
    @abstractmethod
    def name(self):
        pass

    @abstractmethod
    def __call__(self, x, y):
        pass
