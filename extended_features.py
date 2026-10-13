"""Sequence-independent substitution descriptors; no labels or fitted statistics."""
import math
from functools import lru_cache
from Bio.Align import substitution_matrices
from analysis import AMINO_ACID_PROPERTIES
from prepare_dataset import MODEL_FEATURES

EXTENDED_FEATURES = MODEL_FEATURES + [
    'absolute_hydrophobicity_change', 'absolute_weight_change', 'blosum62_score',
]


@lru_cache(maxsize=1)
def blosum62():
    return substitution_matrices.load('BLOSUM62')


def extended_values(row):
    ref, alt = row['ref'], row['alt']
    if ref not in AMINO_ACID_PROPERTIES or alt not in AMINO_ACID_PROPERTIES or ref == alt:
        raise ValueError('Expected a canonical amino-acid substitution')
    original, new = AMINO_ACID_PROPERTIES[ref], AMINO_ACID_PROPERTIES[alt]
    delta_h = new['hydrophobicity'] - original['hydrophobicity']
    delta_w = new['weight'] - original['weight']
    if not all(math.isclose(float(row[name]), value, abs_tol=1e-8) for name, value in
               [('hydrophobicity_change', delta_h), ('weight_change', delta_w)]):
        raise ValueError('Stored features disagree with amino-acid substitution')
    return [float(row[f]) for f in MODEL_FEATURES] + [abs(delta_h), abs(delta_w), float(blosum62()[ref, alt])]
