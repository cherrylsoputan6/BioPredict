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

def download_sequence(gene_id):

    handle = Entrez.efetch(
        db="nucleotide",
        id=gene_id,
        rettype="fasta",
        retmode="text"
    )

    record = SeqIO.read(handle, "fasta")

    return record

def get_description(gene_id):

    handle = Entrez.efetch(
        db="nucleotide",
        id=gene_id,
        rettype="fasta",
        retmode="text"
    )

    record = SeqIO.read(handle, "fasta")

    return record.description