import os


__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__all__ = ['parse_amino_charmff', 'parse_nb_charmff']


def parse_amino_charmff(filename):
    assert os.path.isfile(filename)
    atoms, bonds, donors, acceptors = dict(), dict(), dict(), dict()
    residue_name, atomname, atom_type, charge = None, None, None, None

    with open(filename,"r") as f:
        lines = f.readlines()

    for line in lines:
        if "!" in line:
            line = line[:line.index('!')].strip()
        if line.startswith('RESI'):
            fields = line.split()
            if len(fields) == 3:
                residue_name = fields[1]

        if line.startswith('ATOM'):
            assert residue_name is not None
            fields = line.split()
            if len(fields) == 4:
                atomname, atom_type, charge = fields[1], fields[2], float(fields[3])
                if residue_name not in atoms:
                    atoms[residue_name] = dict()
                atoms[residue_name][atomname] = (atom_type, charge)

        if line.startswith('BOND'):
            fields = line.split()
            assert len(fields) % 2 == 1
            n = len(fields) // 2
            if residue_name not in bonds:
                bonds[residue_name] = []
            for i in range(n):
                bonds[residue_name].append((fields[2*i+1], fields[2*i+2]))

        if line.startswith('DONOR'):
            fields = line.split()
            assert len(fields) > 1
            if residue_name not in donors:
                donors[residue_name] = list()
            for i in range(1, len(fields)):
                donors[residue_name].append(fields[i])

        if line.startswith('ACCEPTOR'):
            fields = line.split()
            assert len(fields) > 1
            if residue_name not in acceptors:
                acceptors[residue_name] = list()
            for i in range(1, len(fields)):
                acceptors[residue_name].append(fields[i])
    return atoms, bonds, donors, acceptors


def parse_nb_charmff(filename):
    assert os.path.isfile(filename)
    with open(filename, "r") as f:
        lines = f.readlines()

    read = False
    parameters = dict()

    for line in lines:
        if "!" in line:
            line = line[:line.index("!")].strip()
        if len(line) < 2:
            continue

        if read is True:
            fields = line.split()
            epsilon14, rmin14 = 0, 0
            if len(fields) != 4 and len(fields) != 7:
                continue
            atom_type, epsilon, rmin = fields[0], float(fields[2]), float(fields[3])
            if len(fields) == 7:
                epsilon14, rmin14 = float(fields[5]), float(fields[6])
            parameters[atom_type] = (epsilon, rmin, epsilon14, rmin14)

        if line.startswith('NONBONDED'):
            read = True
    return parameters