from Bio.Seq import Seq

def calculate_gc(sequence):
    return (
        sequence.count("G") +
        sequence.count("C")
    ) / len(sequence)


def nucleotide_counts(sequence):
    return {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C"),
    }


def reverse_complement(sequence):
    return Seq(sequence).reverse_complement()


def translate(sequence):
    trimmed = sequence[:len(sequence) - (len(sequence) % 3)]
    return Seq(trimmed).translate()