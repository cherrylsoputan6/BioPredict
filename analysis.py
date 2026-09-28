from utils import (
    calculate_gc,
    nucleotide_counts,
    reverse_complement
)


def analyze_sequence(sequence):
    gc = calculate_gc(sequence)
    counts = nucleotide_counts(sequence)

    return {
        "length": len(sequence),
        "gc": gc,
        "counts": counts,
        "reverse": reverse_complement(sequence)
    }


def extract_cds(record):
    for feature in record.features:

        if feature.type == "CDS":
            return feature.extract(record.seq)

    return None


def translate_cds(cds_sequence):
    if cds_sequence is None:
        return None

    return cds_sequence.translate(to_stop=True)