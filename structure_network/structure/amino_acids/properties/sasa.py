__all__ = ['aa_sasa_free', 'aa_sasa_folded']


# Residue accessible surface area in tripeptide (Chothia, 1976)
aa_sasa_free = {"A": 115.000, "C": 135.000, "D": 150.000,
                "E": 190.000, "F": 210.000, "G":  75.000,
                "H": 195.000, "I": 175.000, "K": 200.000,
                "L": 170.000, "M": 185.000, "N": 160.000,
                "P": 145.000, "Q": 180.000, "R": 225.000,
                "S": 115.000, "T": 140.000, "V": 155.000,
                "W": 255.000, "Y": 230.000}


# Residue accessible surface area in folded protein (Chothia, 1976)
aa_sasa_folded = {"A": 25.000, "C": 19.000, "D": 50.000,
                  "E": 49.000, "F": 24.000, "G": 23.000,
                  "H": 43.000, "I": 18.000, "K": 97.000,
                  "L": 23.000, "M": 31.000, "N": 63.000,
                  "P": 50.000, "Q": 71.000, "R": 90.000,
                  "S": 44.000, "T": 47.000, "V": 18.000,
                  "W": 32.000, "Y": 60.000}
