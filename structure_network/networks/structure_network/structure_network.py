from abc import ABC, abstractmethod

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['ProteinStructureNetwork']


class ProteinStructureNetwork(ABC):
    @property
    @abstractmethod
    def size(self):
        pass

    @property
    @abstractmethod
    def residue_ids(self):
        """
        Returns unique residue ids in the structure, in ascending order.
        """
        pass

    @abstractmethod
    def xyz(self, r=None):
        pass

    @abstractmethod
    def get_amino(self, r):
        pass

    @property
    @abstractmethod
    def max_edge_weight(self):
        pass

    @property
    @abstractmethod
    def min_edge_weight(self):
        pass

    @property
    @abstractmethod
    def structure(self):
        pass

    @property
    @abstractmethod
    def g(self):
        """
        Returns networkx graph instance
        """
        pass

    @property
    @abstractmethod
    def residue_names(self):
        """
        Returns name of the residue nodes in the graph. The order
        of residue_names must match the order of residue_ids.
        """
        pass

    @property
    @abstractmethod
    def threshold(self):
        """
        All threshold cutoffs used for the generating the
        network.
        """
        pass

    @property
    @abstractmethod
    def threshold_type(self):
        """
        Threshold type used for building the contact network.
        """
        pass

    @property
    @abstractmethod
    def contact_type(self):
        """
        Specific contact pattern to be used for building  the network.
        """
        pass

    @property
    @abstractmethod
    def centroid(self):
        """
        Mass center of the molecule. x, y, z coordinates
        """
        pass

    @property
    @abstractmethod
    def radius_of_gyration(self):
        """
        Radius of gyration for the protein. It is weighted distance
        of each mass center from the centroid of the molecule. Centroid
        is the mass center of the molecule.
        """
        pass

    @property
    @abstractmethod
    def network_diameter(self):
        """
        In the derived structure network maximum shortest path length.
        The expected analysis is unweighted, and only suggesting of the
        protein structure topology.
        """
        pass

    @property
    @abstractmethod
    def structure_diameter(self):
        """
        Maximum distance between any entity pair compose the structure.
        """
        pass

    @abstractmethod
    def copy(self):
        """
        Deep copy of itself.
        """
        pass

    @property
    @abstractmethod
    def params(self):
        """
        Returns dictionary of key value pairs of internal parameters
        """
        pass
