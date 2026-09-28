from Bio import Entrez
from Bio import SeqIO

Entrez.email = "cherrylsoputan@gmail.com"

def search_gene(gene_name):

    handle = Entrez.esearch(
        db="gene",
        term=f"{gene_name}[Gene Name] AND Homo sapiens[Organism]",
        retmax=1
    )

    results = Entrez.read(handle)
    handle.close()

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

def search_gene_summary(gene_name):

    handle = Entrez.esearch(
        db="nucleotide",
        term=f"{gene_name}[Gene] AND Homo sapiens[Organism] AND refseq[filter] AND biomol_mrna[PROP]",
        retmax=5
    )

    results = Entrez.read(handle)
    handle.close()

    summaries = []

    for gene_id in results["IdList"]:

        summary_handle = Entrez.esummary(
            db="nucleotide",
            id=gene_id
        )

        summary = Entrez.read(summary_handle)[0]
        summary_handle.close()

        summaries.append(summary)

    return summaries

def get_transcript_ids(gene_id):

    handle = Entrez.elink(
        dbfrom="gene",
        db="nucleotide",
        id=gene_id
    )

    results = Entrez.read(handle)
    handle.close()

    transcript_ids = []

    for link_set in results[0]["LinkSetDb"]:
        for link in link_set["Link"]:
            transcript_ids.append(link["Id"])

    return list(dict.fromkeys(transcript_ids))

def find_refseq_transcripts(nucleotide_ids):

    transcripts = []
    batch_size = 100

    for i in range(0, len(nucleotide_ids), batch_size):

        batch = nucleotide_ids[i:i + batch_size]

        handle = Entrez.esummary(
            db="nucleotide",
            id=",".join(batch)
        )

        summaries = Entrez.read(handle)
        handle.close()

        for summary in summaries:

            accession = str(summary["AccessionVersion"])

            if accession.startswith("NM_"):
                transcripts.append(accession)

    return list(dict.fromkeys(transcripts))

def find_mane_transcript(gene_name):
        handle = Entrez.esearch(
            db="nucleotide",
            term=f"{gene_name}[Gene] AND Homo sapiens[Organism] AND MANE Select[keyword]",
            retmax=1
        )

        results = Entrez.read(handle)
        handle.close()

        if not results["IdList"]:
            return None

        nucleotide_id = results["IdList"][0]

        handle = Entrez.esummary(
            db="nucleotide",
            id=nucleotide_id
        )

        summary = Entrez.read(handle)[0]
        handle.close()

        return str(summary["AccessionVersion"])

def download_genbank_record(accession):

    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="gb",
        retmode="text"
    )

    record = SeqIO.read(handle, "genbank")
    handle.close()

    return record