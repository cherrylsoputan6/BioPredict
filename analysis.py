from utils import (
    calculate_gc,
    nucleotide_counts,
    reverse_complement,
    translate
)


def analyze_sequence(sequence):

    gc = calculate_gc(sequence)

    counts = nucleotide_counts(sequence)

    protein = translate(sequence)

    return {
        "length": len(sequence),
        "gc": gc,
        "counts": counts,
        "protein": protein,
        "reverse": reverse_complement(sequence)
    }