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

def extract_cds(record):

    for feature in record.features:

        if feature.type == "CDS":
            cds_sequence = feature.extract(record.seq)

            return cds_sequence

    return None

def translate_cds(cds_sequence):

    if cds_sequence is None:
        return None

    protein = cds_sequence.translate(to_stop=True)

    return protein