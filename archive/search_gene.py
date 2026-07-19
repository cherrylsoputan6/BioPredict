from Bio import Entrez

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

if ids:
    print("First ID:", ids[0])
else:
    print("No genes found.")