__all__ = ['atoms',
           'bonds',
           'hb_donor',
           'hb_acceptor']

atoms = {
    'A': ["C", "CA", "N", "O", "CB"],
    'C': ["C", "CA", "N", "O", "CB", "SG"],
    'D': ["C", "CA", "N", "O", "CB", "CG", "OD1", "OD2"],
    'E': ["C", "CA", "N", "O", "CB", "CG", "CD", "OE1", "OE2"],
    'F': ["C", "CA", "N", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    'G': ["C", "CA", "N", "O"],
    'H': ["C", "CA", "N", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    'I': ["C", "CA", "N", "O", "CB", "CG1", "CG2", "CD1"],
    'K': ["C", "CA", "N", "O", "CB", "CG", "CD", "CE", "NZ"],
    'L': ["C", "CA", "N", "O", "CB", "CG", "CD1", "CD2"],
    'M': ["C", "CA", "N", "O", "CB", "CG", "SD", "CE"],
    'N': ["C", "CA", "N", "O", "CB", "CG", "OD1", "ND2"],
    'P': ["C", "CA", "N", "O", "CB", "CG", "CD"],
    'Q': ["C", "CA", "N", "O", "CB", "CG", "CD", "NE2", "OE1"],
    'R': ["C", "CA", "N", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    'S': ["C", "CA", "N", "O", "CB", "OG"],
    'T': ["C", "CA", "N", "O", "CB", "CG2", "OG1"],
    'V': ["C", "CA", "N", "O", "CB", "CG1", "CG2"],
    'W': ["C", "CA", "N", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    'Y': ["C", "CA", "N", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]
    }

bonds = {
    'A': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'),  ('C:0', 'N:1')],
    'C': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'),  ('C:0', 'N:1'), ('CB:0', 'SG:0'),
          ('SG:0', 'SG:*')],
    'D': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'),  ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'OD1:0'), ('CG:0', 'OD2:0')],
    'E': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'),  ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'CD:0'), ('CD:0', 'OE1:0'), ('CD:0', 'OE2:0')],
    'F': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'),  ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'CD1:0'), ('CD1:0', 'CE1:0'), ('CE1:0', 'CZ:0'),
          ('CG:0', 'CD2:0'), ('CD2:0', 'CE2:0'), ('CE2:0', 'CZ:0')],
    'G': [('N:0', 'CA:0'), ('CA:0', 'C:0'), ('C:0', 'O:0'),
          ('C:0', 'N:1')],
    'H': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'ND1:0'), ('ND1:0', 'CE1:0'), ('CG:0', 'CD2:0'),
          ('CD2:0', 'NE2:0'), ('CE1:0', 'NE2:0')],
    'I': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG1:0'),
          ('CB:0', 'CG2:0'), ('CG1:0', 'CD1:0')],
    'K': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'CD:0'), ('CD:0', 'CE:0'), ('CE:0', 'NZ:0')],
    'L': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'CD1:0'), ('CG:0', 'CD2:0')],
    'M': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'SD:0'), ('SD:0', 'CE:0')],
    'N': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'OD1:0'), ('CG:0', 'ND2:0')],
    'P': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'CD:0'), ('CD:0', 'N:0')],
    'Q': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'CD:0'), ('CD:0', 'OE1:0'), ('CD:0', 'NE2:0')],
    'R': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'CD:0'), ('CD:0', 'NE:0'), ('NE:0', 'CZ:0'),
          ('CZ:0', 'NH1:0'), ('CZ:0', 'NH2:0')],
    'S': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'OG:0')],
    'T': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'OG1:0'),
          ('CB:0', 'CG2:0')],
    'V': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG1:0'),
          ('CB:0', 'CG2:0')],
    'W': [('N:0', 'CA:0'),    ('CA:0', 'CB:0'),   ('CA:0', 'C:0'),
          ('C:0', 'O:0'),     ('C:0', 'N:1'),     ('CB:0', 'CG:0'),
          ('CG:0', 'CD1:0'),  ('CG:0', 'CD2:0'),  ('CD2:0', 'CE2:0'),
          ('CD1:0', 'NE1:0'), ('CE2:0', 'CZ2:0'), ('CD2:0', 'CE3:0'),
          ('CZ2:0', 'CH2:0'), ('CH2:0', 'CZ3:0'), ('CE3:0', 'CZ3:0')],
    'Y': [('N:0', 'CA:0'), ('CA:0', 'CB:0'), ('CA:0', 'C:0'),
          ('C:0', 'O:0'), ('C:0', 'N:1'), ('CB:0', 'CG:0'),
          ('CG:0', 'CD1:0'), ('CG:0', 'CD2:0'), ('CD1:0', 'CE1:0'),
          ('CD2:0', 'CE2:0'), ('CE2:0', 'CZ:0'), ('CE1:0', 'CZ:0'),
          ('CZ:0', 'OH:0')]
}

hb_donor = {
    'A': ['N'],
    'C': ['N', 'SG'],
    'D': ['N'],
    'E': ['N'],
    'F': ['N'],
    'G': ['N'],
    'H': ['N', 'ND1'],
    'I': ['N'],
    'K': ['N'],
    'L': ['N'],
    'M': ['N'],
    'N': ['N', 'ND2'],
    'P': [],
    'Q': ['N', 'NE2'],
    'R': ['N', 'NE', 'NH1', 'NH2'],
    'S': ['N', 'OG'],
    'T': ['N', 'OG1'],
    'V': ['N'],
    'W': ['N', 'NE1'],
    'Y': ['N', 'OH']
}

hb_acceptor = {
    'A': ['O'],
    'C': ['O'],
    'D': ['O', 'OD1', 'OD2'],
    'E': ['O', 'OE1', 'OE2'],
    'F': ['O'],
    'G': ['O'],
    'H': ['O', 'NE2'],
    'I': ['O'],
    'K': ['O'],
    'L': ['O'],
    'M': ['O'],
    'N': ['O', 'OD1'],
    'P': ['O'],
    'Q': ['O', 'OE1'],
    'R': ['O'],
    'T': ['O', 'OG1'],
    'V': ['O'],
    'W': ['O'],
    'Y': ['O', 'OH']
}