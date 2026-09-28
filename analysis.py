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

def parse_variant(variant):

    variant = variant.upper().strip()

    original_amino_acid = variant[0]
    new_amino_acid = variant[-1]
    position = int(variant[1:-1])

    return {
        "original": original_amino_acid,
        "position": position,
        "new": new_amino_acid
    }

def validate_variant(protein, variant_info):

    position = variant_info["position"]
    expected_amino_acid = variant_info["original"]

    if position < 1 or position > len(protein):
        return False

    actual_amino_acid = protein[position - 1]

    return actual_amino_acid == expected_amino_acid