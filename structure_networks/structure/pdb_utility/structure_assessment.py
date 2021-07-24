import logging
import numpy as np
from structure_networks.geometry import *
from structure_networks.structure import *

__version__ = "1.0"
__all__ = ['CADistanceReport', 'fix_ca_trace']


class CADistanceReport:
    def __init__(self,
                 allowed_lb=0.2,
                 allowed_ub=0.2):
        assert (allowed_lb > 0) and (allowed_lb < 1.0)
        assert (allowed_ub > 0) and (allowed_ub < 1.0)
        self._chain = None
        self._valid_sequence = []
        self._violated_sequence = []
        self._gapped_sequence = []
        self._min_separation = 3.8 - allowed_lb
        self._max_separation = 3.8 + allowed_ub

    def assess(self, structure):
        assert isinstance(structure, CaTrace)
        self._chain = structure.chain
        residues = structure.residue_ids
        self._violated_sequence.clear()
        self._valid_sequence.clear()
        self._gapped_sequence.clear()
        for i in range(1, len(residues)):
            if (residues[i] - residues[i-1]) > 1:
                self._gapped_sequence.append((residues[i-1], residues[i]))
            crd_0 = Coordinate3d(*structure.xyz(residues[i-1]))
            crd_1 = Coordinate3d(*structure.xyz(residues[i]))
            d = distance(crd_0, crd_1)
            pair = (residues[i-1], residues[i], d)
            if (d > self._max_separation) or (d < self._min_separation):
                self._violated_sequence.append(pair)
            else:
                self._valid_sequence.append(pair)

    @property
    def chain(self):
        return self._chain

    @property
    def n_violations(self):
        return len(self._violated_sequence)

    @property
    def n_valid(self):
        return len(self._valid_sequence)

    @property
    def n_gap(self):
        return len(self._gapped_sequence)

    @property
    def has_sequence_gap(self):
        return (len(self._gapped_sequence) > 0) and (self._chain is not None)

    def get_valid(self, idx):
        assert (idx >= 0) and (idx < self.n_valid)
        return self._valid_sequence[idx]

    def get_violation(self, idx):
        assert (idx >= 0) and (idx < self.n_violations)
        return self._violated_sequence[idx]

    def get_gap(self, idx):
        assert (idx >= 0) and (idx < self.n_gap)
        return self._gapped_sequence[idx]


def shrink_edge(xyz, pos, dist_lb, dist_ub):
    assert isinstance(xyz, list)
    n = len(xyz)
    assert (pos < n - 1)
    assert (dist_lb > 0) and (dist_ub > 0) and (dist_lb < dist_ub)
    dfound = distance(Coordinate3d(*xyz[pos]), Coordinate3d(*xyz[pos+1]))
    if (dfound < dist_lb) or (dfound > dist_ub):
        if dfound < dist_lb:
            d = dfound - dist_lb
        else:
            d = dfound - dist_ub
        u = connecting_vector(Coordinate3d(*xyz[pos]), Coordinate3d(*xyz[pos+1])).unit_vector
        shift_vector = Vector3d(u.x * d, u.y * d, u.z * d)
        for i in range(pos+1):
            xyz[i] = (xyz[i][0] + shift_vector.x,
                      xyz[i][1] + shift_vector.y,
                      xyz[i][2] + shift_vector.z)
    return xyz


def shrink_loop(xyz,
                start_pos,
                end_pos,
                dist):
    assert isinstance(xyz, list)
    n = len(xyz)
    logger = logging.getLogger('structural_dynamics.shrink_loop')
    logger.debug('Start (%d) End (%d)' % (start_pos, end_pos))
    assert (start_pos < end_pos) and (end_pos - start_pos > 1)
    assert (start_pos > 0) and (end_pos < n)
    for i in range(start_pos, end_pos-1):
        p1, p2, p3 = i, i+1, i+2
        pt1, pt2, pt3 = Coordinate3d(*xyz[p1]), Coordinate3d(*xyz[p2]), Coordinate3d(*xyz[p3])
        if distance(pt1, pt3) >= 2 * dist:
            xyz[p2] = ((pt1.x+pt3.x)/2.0, (pt1.y + pt3.y)/2.0, (pt1.z+pt3.z)/2.0)
        else:
            u = (connecting_vector(pt2, pt1) + connecting_vector(pt2, pt3)).unit_vector
            b = 2*((pt1.x - pt2.x)*u.x + (pt1.y - pt2.y)*u.y + (pt1.z - pt2.z)*u.z)
            c = (pt1.x - pt2.x)**2 + (pt1.y - pt2.y)**2 + (pt1.z - pt2.z)**2 - dist**2
            d = (b - np.sqrt(b**2 - 4*c))/2.0
            xyz[p2] = (pt2.x + d*u.x, pt2.y + d*u.y, pt2.z + d*u.z)
    return xyz


def find_centroid(xyz):
    assert isinstance(xyz, list)
    cx, cy, cz = 0, 0, 0
    n = len(xyz)
    for x, y, z in xyz:
        cx += x
        cy += y
        cz += z
    cx, cy, cz = cx/n, cy/n, cz/n
    return cx, cy, cz


def fix_ca_trace(caTrace,
                 allowed_lb=0.2,
                 allowed_ub=0.2):
    logger = logging.getLogger('structural_dynamics.fix_ca_trace')
    assert isinstance(caTrace, CaTrace)
    report = CADistanceReport(allowed_lb=allowed_lb,
                              allowed_ub=allowed_ub)
    report.assess(caTrace)
    assert report.n_gap == 0
    assert report.n_violations < (len(caTrace)//2)
    violation_counts = {r: 0 for r in caTrace.residue_ids}
    for i in range(report.n_violations):
        ri, rj, _ = report.get_violation(i)
        violation_counts[ri] += 1
        violation_counts[rj] += 1
    resids = caTrace.residue_ids
    xyz = [(caTrace.xyz(r)) for r in resids]
    i, start_res = 1, True
    while i < len(resids):
        r = resids[i]
        if (start_res is True) and (violation_counts[r] == 0):
            start_res = False
        if (violation_counts[r] > 0) and (start_res is True):
            logger.debug("Fixing end violation @ (%d)" % r)
            xyz = shrink_edge(xyz, i, 3.8 - allowed_lb, 3.8 + allowed_ub)
        elif (violation_counts[r] > 0) and (start_res is False):
            l = 0
            while (violation_counts[resids[i+l+1]] > 0) and (i+l+1 < len(resids)):
                l = l + 1
            if l > 1:
                logger.debug("Loop violation detected @ (%d - %d)" % (r, r + l))
                xyz = shrink_loop(xyz, i, i+l, 3.8)
                i += l
        i += 1
    return CaTrace(caTrace.name,
                   [{'resid': r,
                     'resname': caTrace.get_amino(r).name(one_letter_code=False),
                     'x': xyz[i][0],
                     'y': xyz[i][1],
                     'z': xyz[i][2]} for i, r in enumerate(resids)],
                   caTrace.chain)


