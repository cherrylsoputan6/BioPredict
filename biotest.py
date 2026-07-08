from Bio.Seq import Seq

dna = Seq("ATGGCC")

print("DNA:", dna)
print("Protein:", dna.translate())
print("Reverse complement:", dna.reverse_complement()) 