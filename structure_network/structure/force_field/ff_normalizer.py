import logging
import numpy as np

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__all__ = ['FFNormalizer']


class FFNormalizer:
    def __init__(self, cutoff=0, potential=True):
        self.__cutoff = cutoff
        self.__potential = potential
        self.__logger = logging.getLogger(name="temporal_network.FFNormalizer")

    def __call__(self, potential):
        if self.__potential:
            value = -1 * potential
        else:
            value = potential
        self.__logger.debug("Value before normalization (%f)" % potential)
        return np.clip(value, a_min=self.__cutoff, a_max=None)

