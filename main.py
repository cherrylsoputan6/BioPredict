from Bio import Entrez
from Bio import SeqIO

Entrez.email = "cherrylsoputan@gmail.com"

def search_gene(gene_name):

    handle = Entrez.esearch(
        db="nucleotide",
        term=f"{gene_name}[Gene] AND Homo sapiens[Organism] AND mRNA",
        retmax=5
    )

    results = Entrez.read(handle)

    return results["IdList"]

gene = input("Enter gene name: ")

ids = search_gene(gene)

for gene_id in ids:

    handle = Entrez.efetch(
        db="nucleotide",
        id=gene_id,
        rettype="fasta",
        retmode="text"
    )

    record = SeqIO.read(handle, "fasta")

    print(record.description)
    print("-" * 50)