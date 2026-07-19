from Bio.Seq import Seq

sequence = input("Enter a DNA sequence: ")

dna = Seq(sequence)

gc_content = (
    sequence.count("G")
    + sequence.count("C")
) / len(sequence)

a_count = sequence.count("A")
t_count = sequence.count("T")
g_count = sequence.count("G")
c_count = sequence.count("C")

print()
print("DNA Analysis Report")
print("===================")
print("Nucleotide Counts:")
print("A:", a_count)
print("T:", t_count)
print("G:", g_count)
print("C:", c_count)

print()
print("Length:", len(sequence), "bp")
print("GC content:", round(gc_content * 100, 2), "%")
print("Reverse complement:")
print(dna.reverse_complement())
print("Protein translation:")
print(dna.translate())

if gc_content > 0.60:
    print("\nInterpretation:")
    print("High GC content.")
elif gc_content < 0.40:
    print("\nInterpretation:")
    print("Low GC content.")
else:
    print("\nInterpretation:")
    print("Moderate GC content.")