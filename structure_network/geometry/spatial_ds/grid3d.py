import numpy as np
from copy import deepcopy
from structure_network.geometry.coordinate import *


__all__ = ['Grid3D', 'DistanceCutoff']


class DistanceCutoff:
    def __init__(self, def_cutoff=8.0):
        self.__cutoff = def_cutoff

    @property
    def cutoff(self):
        return self.__cutoff

    def __ge__(self, other):
        if isinstance(other, np.float) or isinstance(other, np.int) or isinstance(other, np.double):
            return self.__cutoff > other
        elif isinstance(other, DistanceCutoff):
            return self.__cutoff > other.cutoff
        raise Exception('Invalid comparison!!')

    def __eq__(self, other):
        if isinstance(other, np.float) or isinstance(other, np.int) or isinstance(other, np.double):
            return np.abs(self.__cutoff - other) < 1e-3
        elif isinstance(other, DistanceCutoff):
            return np.abs(self.__cutoff - other.cutoff) < 1e-3
        return Exception("Invalid equality comparison!!")

    def __float__(self):
        return self.__cutoff

    def __call__(self, context_from=None, context_to=None):
        return self.__cutoff


class Grid3D:
    def __init__(self, min_coord, max_coord, spacing=1.):
        assert isinstance(min_coord, Coordinate3d)
        assert isinstance(max_coord, Coordinate3d)
        assert spacing > 0.1
        assert (max_coord.x - min_coord.x) > spacing
        assert (max_coord.y - min_coord.y) > spacing
        assert (max_coord.z - min_coord.z) > spacing
        self.__max_coordinate = deepcopy(max_coord)
        self.__min_coordinate = deepcopy(min_coord)
        self.__spacing = spacing
        self.__nx = int((self.__max_coordinate.x - self.__min_coordinate.x + spacing) // spacing)
        self.__ny = int((self.__max_coordinate.y - self.__min_coordinate.y + spacing) // spacing)
        self.__nz = int((self.__max_coordinate.z - self.__min_coordinate.z + spacing) // spacing)
        self.__registered_id = dict()
        self.__grid = dict()
        self.__occupied_cells = 0

    def volume(self):
        return (self.__max_coordinate.x - self.__min_coordinate.x) * \
               (self.__max_coordinate.y - self.__min_coordinate.y) * \
               (self.__max_coordinate.z - self.__min_coordinate.z)

    def cell_counts(self):
        return self.__nx * self.__ny * self.__nz

    @property
    def occupied_cells(self):
        return self.__occupied_cells

    def __x_bound(self, coord):
        return (coord.x <= self.__max_coordinate.x) and \
               (self.__min_coordinate.x <= coord.x)

    def __y_bound(self, coord):
        return (coord.y <= self.__max_coordinate.y) and \
               (self.__min_coordinate.y <= coord.y)

    def __z_bound(self, coord):
        return (coord.z <= self.__max_coordinate.z) and \
               (self.__min_coordinate.z <= coord.z)

    def inside(self, coord):
        if isinstance(coord, Coordinate3d):
            return self.__x_bound(coord) and \
                   self.__y_bound(coord) and \
                   self.__z_bound(coord)
        return False

    def __grid_cell(self, coord):
        assert self.inside(coord)
        xi = int((coord.x - self.__min_coordinate.x) // self.__spacing)
        yi = int((coord.y - self.__min_coordinate.y) // self.__spacing)
        zi = int((coord.z - self.__min_coordinate.z) // self.__spacing)
        return xi, yi, zi

    def register_obj(self, id, coord):
        if self.inside(coord) and id not in self.__registered_id:
            xi, yi, zi = self.__grid_cell(coord)
            self.__registered_id[id] = coord
            if xi not in self.__grid:
                self.__grid[xi] = dict()
            if yi not in self.__grid[xi]:
                self.__grid[xi][yi] = dict()
            if zi not in self.__grid[xi][yi]:
                self.__grid[xi][yi][zi] = []
                self.__occupied_cells += 1
            self.__grid[xi][yi][zi].append(id)
            return True
        return False

    def is_registered(self, id):
        return id in self.__registered_id

    @property
    def register_count(self):
        return len(self.__registered_id)

    def neighbor_cells(self, id, radius):
        r = self.__spacing if radius < self.__spacing else radius
        nbr_cells = []
        if self.is_registered(id):
            xi, yi, zi = self.__grid_cell(self.__registered_id[id])
            steps = int(r // self.__spacing)
            x_min = xi - steps if xi > steps else 0
            x_max = xi + steps if xi + steps < self.__nx else self.__nx
            y_min = yi - steps if yi > steps else 0
            y_max = yi + steps if yi + steps < self.__ny else self.__ny
            z_min = zi - steps if zi > steps else 0
            z_max = zi + steps if zi + steps < self.__nz else self.__nz
            for x in range(x_min, x_max):
                for y in range(y_min, y_max):
                    for z in range(z_min, z_max):
                        nbr_cells.append( (x,y,z) )
            return nbr_cells

    def is_occupied(self, xi, yi, zi):
        _x, _y, _z = int(xi), int(yi), int(zi)
        if (_x in self.__grid) and (_y in self.__grid[_x]) and (_z in self.__grid[_x][_y]):
            return len(self.__grid[_x][_y][_z]) > 0
        return False

    def neighbors(self, id, radius):
        nbr_ids = []
        if self.is_registered(id):
            nbr_cells = self.neighbor_cells(id, radius)
            crd = self.__registered_id[id]
            for x, y, z in nbr_cells:
                if self.is_occupied(x, y, z):
                    nbr_ids = nbr_ids + self.__grid[x][y][z]
            nbr_ids = [n for n in list(set(nbr_ids)) if n != id and distance(crd, self.__registered_id[n]) <= radius]
        return nbr_ids

