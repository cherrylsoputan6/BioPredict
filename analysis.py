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

AMINO_ACID_PROPERTIES = {
    "A": {"hydrophobicity": 1.8, "weight": 89.09},
    "R": {"hydrophobicity": -4.5, "weight": 174.20},
    "N": {"hydrophobicity": -3.5, "weight": 132.12},
    "D": {"hydrophobicity": -3.5, "weight": 133.10},
    "C": {"hydrophobicity": 2.5, "weight": 121.16},
    "Q": {"hydrophobicity": -3.5, "weight": 146.15},
    "E": {"hydrophobicity": -3.5, "weight": 147.13},
    "G": {"hydrophobicity": -0.4, "weight": 75.07},
    "H": {"hydrophobicity": -3.2, "weight": 155.16},
    "I": {"hydrophobicity": 4.5, "weight": 131.17},
    "L": {"hydrophobicity": 3.8, "weight": 131.17},
    "K": {"hydrophobicity": -3.9, "weight": 146.19},
    "M": {"hydrophobicity": 1.9, "weight": 149.21},
    "F": {"hydrophobicity": 2.8, "weight": 165.19},
    "P": {"hydrophobicity": -1.6, "weight": 115.13},
    "S": {"hydrophobicity": -0.8, "weight": 105.09},
    "T": {"hydrophobicity": -0.7, "weight": 119.12},
    "W": {"hydrophobicity": -0.9, "weight": 204.23},
    "Y": {"hydrophobicity": -1.3, "weight": 181.19},
    "V": {"hydrophobicity": 4.2, "weight": 117.15}
}

def extract_variant_features(protein, variant_info):

    original = variant_info["original"]
    new = variant_info["new"]
    position = variant_info["position"]

    original_properties = AMINO_ACID_PROPERTIES[original]
    new_properties = AMINO_ACID_PROPERTIES[new]

    relative_position = position / len(protein)

    hydrophobicity_change = (
        new_properties["hydrophobicity"]
        - original_properties["hydrophobicity"]
    )

    weight_change = (
        new_properties["weight"]
        - original_properties["weight"]
    )

    return {
        "position": position,
        "protein_length": len(protein),
        "relative_position": relative_position,
        "original_hydrophobicity": original_properties["hydrophobicity"],
        "new_hydrophobicity": new_properties["hydrophobicity"],
        "hydrophobicity_change": hydrophobicity_change,
        "original_weight": original_properties["weight"],
        "new_weight": new_properties["weight"],
        "weight_change": weight_change
    }