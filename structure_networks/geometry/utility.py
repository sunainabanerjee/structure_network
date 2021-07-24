import math
import numpy as np
from .coordinate import *
from .vector_algebra import *

__version__ = "1.0"
__all__ = ['cuboid_volume',
           'sphere_volume',
           'sphere_intersection_volume',
           'equivalent_sphere_radius',
           'angle',
           'dihedral',
           'rotation_matrix',
           'reconstruct_coordinate']


def cuboid_volume(sx, sy, sz):
    assert (sx >= 0) and (sy >= 0) and (sz >= 0)
    return sx * sy * sz


def sphere_volume(r):
    assert r >= 0
    return 4.0 * math.pi * (r**3) / 3.0


def sphere_intersection_volume(r1, center1, r2, center2):
    assert isinstance(center1, Coordinate3d) and r1 > 0
    assert isinstance(center2, Coordinate3d) and r2 > 0
    d = distance(center1, center2)
    if d > (r1 + r2):
        return 0
    if d < abs(r1 - r2):
        return sphere_volume(min(r1, r2))
    f = (d**2 + 2*d*(r1 + r2) - 3*(r1**2 + r2**2) + 6*r1*r2) / (12*d)
    return math.pi * ((r1 + r2 - d)**2) * f


def equivalent_sphere_radius(vol):
    assert vol >= 0.
    return ((vol * 3)/(4 * math.pi))**(1.0/3.0)


def angle(c1, c2, c3):
    assert isinstance(c1, Coordinate3d) or len(c1) == 3
    assert isinstance(c2, Coordinate3d) or len(c2) == 3
    assert isinstance(c3, Coordinate3d) or len(c3) == 3
    return math.acos(dotp(connecting_vector(c1, c2).unit_vector,
                          connecting_vector(c2, c3).unit_vector))


def dihedral(c1, c2, c3, c4):
    assert isinstance(c1, Coordinate3d) or len(c1) == 3
    assert isinstance(c2, Coordinate3d) or len(c2) == 3
    assert isinstance(c3, Coordinate3d) or len(c3) == 3
    assert isinstance(c4, Coordinate3d) or len(c4) == 3
    b1 = connecting_vector(c2, c1)
    b2 = connecting_vector(c2, c3).unit_vector
    b3 = connecting_vector(c3, c4)
    v = b1 - b2 * dotp(b1, b2)
    w = b3 - b2 * dotp(b3, b2)
    x = dotp(v, w)
    y = dotp(crossp(b2, v), w)
    return math.atan2(y, x)


def rotation_matrix(v, theta):
    if not isinstance(v, Vector3d) and len(v) == 3:
        v = Vector3d(v[0], v[1], v[2])
    assert isinstance(v, Vector3d)
    u = v.unit_vector
    m = np.zeros((len(u), len(u)), dtype=np.double)
    ct, st = np.cos(theta), np.sin(theta)
    m[0, 0] = ct + u.x*u.x*(1. - ct)
    m[0, 1] = u.x*u.y*(1. - ct) - u.z*st
    m[0, 2] = u.x*u.z*(1. - ct) + u.y*st

    m[1, 0] = u.x*u.y*(1. - ct) + u.z*st
    m[1, 1] = ct + u.y*u.y*(1. - ct)
    m[1, 2] = u.y*u.z*(1. - ct) - u.x*st

    m[2, 0] = u.z*u.x*(1. - ct) - u.y*st
    m[2, 1] = u.z*u.y*(1. - ct) + u.x*st
    m[2, 2] = ct + u.z*u.z*(1. - ct)
    return m


def reconstruct_coordinate(c1, c2, c3, dist, angle, dihed):
    assert isinstance(c1, Coordinate3d) or (len(c1) == 3)
    assert isinstance(c2, Coordinate3d) or (len(c2) == 3)
    assert isinstance(c3, Coordinate3d) or (len(c3) == 3)
    assert dist > 0
    v23 = connecting_vector(c2, c3)
    n321 = crossp(connecting_vector(c2, c3), connecting_vector(c1, c2))
    dihed_rot = rotation_matrix(v23, dihed)
    bond_rot = rotation_matrix(n321, -angle)
    m = np.matmul(dihed_rot, bond_rot)
    u = connecting_vector(c2, c3).unit_vector.toarray()
    v = np.matmul(m, u)
    return Coordinate3d(c3[0] + v[0]*dist,
                        c3[1] + v[1]*dist,
                        c3[2] + v[2]*dist)

