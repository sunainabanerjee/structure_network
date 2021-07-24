import numpy as np
from copy import deepcopy
from .coordinate import Coordinate3d

__version__ = "1.0"
__all__ = ['rmsd',
           'kabsch_coordinate',
           'kabsch_rmsd',
           'rotation_matrix3d']


def rmsd(V, W):
    v = np.array(V, dtype=np.float)
    w = np.array(W, dtype=np.float)
    assert v.shape == w.shape
    if len(v.shape) == 1:
        return np.sqrt(np.mean(np.square(v - w)))
    else:
        return np.sqrt(np.mean(np.square(v - w).reshape(v.shape[0], -1), axis=1))


def kabsch_rmsd(x, y,
                return_matrix=False,
                return_coordinate=False):
    assert isinstance(x, np.ndarray)
    assert isinstance(y, np.ndarray)
    assert (len(x.shape) == len(y.shape)) and len(x.shape) == 2
    assert (x.shape[0] == y.shape[0]) and (x.shape[1] == y.shape[1]) and (x.shape[0] > x.shape[1])
    xc = x - x.mean(axis=0)
    yc = y - y.mean(axis=0)
    C = np.dot(np.transpose(xc), yc)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = np.dot(V, W)
    xca = np.matmul(xc, U)
    r = rmsd(xca, yc)
    if (not return_matrix) and (not return_coordinate):
        return r
    elif (not return_matrix) and return_coordinate:
        return r, xca + y.mean(axis=0)
    elif return_matrix and (not return_coordinate):
        return r, U
    else:
        return r, U, xca + y.mean(axis=0)


def kabsch_coordinate(x1, x2, return_matrix=False, return_coordinate=False):
    if isinstance(x1, list):
        x1a = np.array([[x1[i][0], x1[i][1], x1[i][2]] for i in range(len(x1))])
    else:
        x1a = deepcopy(x1)
    if isinstance(x2, list):
        x2a = np.array([[x2[i][0], x2[i][1], x2[i][2]] for i in range(len(x2))])
    else:
        x2a = deepcopy(x2)

    assert isinstance(x1a, np.ndarray) and isinstance(x2a, np.ndarray)
    assert x1a.shape == x2a.shape

    if not return_coordinate:
        return kabsch_rmsd(x1a, x2a,
                           return_matrix=return_matrix,
                           return_coordinate=return_coordinate)
    else:
        d, m, x1a = kabsch_rmsd(x1a, x2a,
                                return_matrix=return_matrix,
                                return_coordinate=return_coordinate)
        return d, m, [Coordinate3d(x1a[i, 0], x1a[i, 1], x1a[i, 2]) for i in range(x1a.shape[0])]


def rotation_matrix3d(axis_vector, rotation_angle):
    assert len(axis_vector) == 3
    r = np.zeros([3, 3])
    n = np.sqrt(axis_vector[0] ** 2 + axis_vector[1] ** 2 + axis_vector[2] ** 2)
    n = 1. if n < 1e-6 else n
    ux, uy, uz = axis_vector[0] / n, axis_vector[1] / n, axis_vector[2] / n
    ct, st = np.cos(rotation_angle), np.sin(rotation_angle)

    r[0, 0] = ct + ux * ux * (1 - ct)
    r[0, 1] = ux * uy * (1 - ct) - uz * st
    r[0, 2] = ux * uz * (1 - ct) + uy * st

    r[1, 0] = ux * uy * (1 - ct) + uz * st
    r[1, 1] = ct + uy * uy * (1 - ct)
    r[1, 2] = uy * uz * (1 - ct) - ux * st

    r[2, 0] = uz * ux * (1 - ct) - uy * st
    r[2, 1] = uz * uy * (1 - ct) + ux * st
    r[2, 2] = ct + uz * uz * (1 - ct)
    return r

