import os
import logging
from structure_network.structure.pdb_processor.pdb_structure import CaTrace

__all__ = ['CompressedTrajectoryStore']
__version__ = "1.0"


class CompressedTrajectoryStore:
    @staticmethod
    def save(trajectory, filename):
        if os.path.isfile(filename):
            logging.warning("[%s] file will be overwritten" % filename)

        if isinstance(trajectory, CaTrace):
            trajectory = [trajectory]
        assert isinstance(trajectory, list)
        assert len(trajectory) > 0
        assert all([isinstance(snp, CaTrace) and len(trajectory[0]) == len(snp) for snp in trajectory])

        with open(filename, 'w') as fp:
            residue_ids = trajectory[0].residue_ids
            name = trajectory[0].name
            n = len(trajectory[0])
            chain = trajectory[0].chain
            id_str = ";".join(['%d' % r for r in residue_ids])
            name_str = ";".join(['%s' % trajectory[0].get_amino(r) for r in residue_ids])
            fp.write("%s;%d;%s\n" % (name, n, chain))
            fp.write("%s\n" % id_str)
            fp.write("%s\n" % name_str)
            for i, snp in enumerate(trajectory):
                data = [(*snp.xyz(r), snp.b_factor(r)) for r in residue_ids]
                data_str = ";".join(['%.3f,%.3f,%.3f,%.2f' % (x, y, z, b) for x, y, z, b in data])
                fp.write("%s\n" % data_str)

    @staticmethod
    def load(filename):
        assert os.path.isfile(filename)
        trajectory = list()
        with open(filename, "r") as fp:
            counter = 0
            for line in fp:
                line = line.strip()
                if counter == 0:
                    name, n, chain = line.split(";")
                    n = int(n)
                elif counter == 1:
                    residue_ids = [int(r) for r in line.split(";")]
                    assert n == len(residue_ids)
                elif counter == 2:
                    residue_names = line.split(";")
                    assert n == len(residue_names)
                else:
                    data =[ entry.split(",") for entry in line.split(";") ]
                    assert n == len(data)
                    entry = [{'x': float(entry[0]),
                              'y': float(entry[1]),
                              'z': float(entry[2]),
                              'bfactor': float(entry[3]),
                              'resname': residue_names[i],
                              'resid': residue_ids[i]} for i, entry in enumerate(data)]
                    trajectory.append(CaTrace(name=name, entry=entry, chainId=chain))
                counter = counter + 1
        return trajectory

