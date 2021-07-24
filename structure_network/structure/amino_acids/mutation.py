import re
from .amino_acids import get_amino

__author__ = ['sumanta mukherjee', 'sunaina banerjee']
__email__ = ["mukherjee.sumanta.iisc@gmail.com", "ms.sunaina.b@gmail.com"]
__all__ = ['Mutation', 'parse_mutation']


class Mutation:
    def __init__(self,
                 residue_id,
                 base_amino,
                 mutated_amino):
        base_amino = get_amino(base_amino)
        mutated_amino = get_amino(mutated_amino)
        assert base_amino != mutated_amino, "Error: not a valid mutation!"
        self._residue_id = residue_id
        self._base = base_amino
        self._mutant = mutated_amino

    @property
    def residue_id(self):
        return self._residue_id

    @property
    def base(self):
        return self._base

    @property
    def mutation(self):
        return self._mutant

    @property
    def name(self):
        return '{}{}{}'.format(self._base.name(one_letter_code=True),
                               self._residue_id,
                               self._mutant.name(one_letter_code=True))

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def __format__(self, **kwargs):
        return self.name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __lt__(self, other):
        return self.name < other.name


def parse_mutation(mutation_expr):
    mutation_expr = mutation_expr.upper()
    if (len(mutation_expr) == 0) or (not re.match('^[A-Z][0-9]+[A-Z]$', mutation_expr)):
        raise RuntimeError("Error: not a valid mutation expression {}".format(mutation_expr))
    base_amino = get_amino(mutation_expr[0])
    mut_amino = get_amino(mutation_expr[-1])
    residue_id = int(mutation_expr[1:-1])
    return Mutation(residue_id=residue_id, base_amino=base_amino, mutated_amino=mut_amino)
