from Bio import Entrez
from Bio import SeqIO
from utils import (
    calculate_gc,
    nucleotide_counts,
    reverse_complement,
    translate,
)
Entrez.email = "cherrylsoputan@gmail.com"

gene_id = "NM_000546"

handle = Entrez.efetch(
    db="nucleotide",
    id=gene_id,
    rettype="fasta",
    retmode="text"
)

record = SeqIO.read(handle, "fasta")

sequence = str(record.seq)

gc_content = calculate_gc(sequence)

print()
print("Gene:")
print(record.description)

print()
print(f"Length: {len(sequence)} bp")

print()
print(f"GC Content: {gc_content * 100:.2f}%")

counts = nucleotide_counts(sequence)

print("A:", counts["A"])
print("T:", counts["T"])
print("G:", counts["G"])
print("C:", counts["C"])

print()
print("First 100 bases:")
print(sequence[:100])

print()
print("First 60 amino acids:")

protein = translate(sequence)

print(protein[:60])

print()
print("Protein length:")
print(len(protein), "amino acids")

