"""Fixed local sequence descriptors from the exact validated MANE protein."""
import math
from collections import Counter
from functools import lru_cache
from pathlib import Path
from Bio import SeqIO
from analysis import AMINO_ACID_PROPERTIES
from validate_dataset import checked_protein

CONTEXT_FEATURES = ['neighbor_mean_hydrophobicity', 'neighbor_glycine_proline_fraction', 'neighbor_entropy']


def context_values(protein, position, radius=5):
    if not 1 <= position <= len(protein):
        raise ValueError('Context position outside protein')
    center = position - 1
    neighbors = protein[max(0, center - radius):center] + protein[center + 1:center + 1 + radius]
    if not neighbors or any(aa not in AMINO_ACID_PROPERTIES for aa in neighbors):
        raise ValueError('Unsupported or missing local amino-acid context')
    counts = Counter(neighbors)
    size = len(neighbors)
    entropy = -sum((n / size) * math.log2(n / size) for n in counts.values())
    return [sum(AMINO_ACID_PROPERTIES[aa]['hydrophobicity'] for aa in neighbors) / size,
            sum(aa in 'GP' for aa in neighbors) / size, entropy]


@lru_cache(maxsize=1000)
def cached_protein(gene, transcript):
    import re
    if not re.fullmatch(r'NM_\d+\.\d+', transcript):
        raise ValueError('Invalid transcript accession')
    record = SeqIO.read(Path('data/sequence_cache') / f'{transcript}.gb', 'genbank')
    return checked_protein(record, gene, transcript)


def row_context(row):
    protein = cached_protein(row['gene'], row['transcript'])
    position = int(row['position'])
    if not 1 <= position <= len(protein) or protein[position - 1] != row['ref']:
        raise ValueError('Sequence context reference differs from validated row')
    return context_values(str(protein), position)
