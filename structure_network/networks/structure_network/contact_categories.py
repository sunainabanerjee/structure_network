from enum import Enum

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['ThresholdTypes',
           'ContactTypes']


class ThresholdTypes(Enum):
    distance = 'distance'
    energy = 'energy'
    distance_energy = 'all'


class ContactTypes(Enum):
    backbone_only = 'bb'
    non_bonded_only = 'nb'
    all_contacts = 'all'
