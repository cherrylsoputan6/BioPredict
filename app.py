from ncbi import (
    search_gene,
    find_mane_transcript,
    download_genbank_record
)

from analysis import extract_cds, translate_cds


gene = input("Enter gene: ")

ids = search_gene(gene)

if ids:
    gene_id = ids[0]

    print()
    print("Gene:", gene.upper())
    print("NCBI Gene ID:", gene_id)

    transcript = find_mane_transcript(gene)

    if transcript:
        print("MANE Select transcript:", transcript)

        print("Downloading GenBank record...")
        record = download_genbank_record(transcript)

        cds = extract_cds(record)
        protein = translate_cds(cds)

        if cds:
            print()
            print("Transcript length:", len(record.seq), "bp")
            print("CDS length:", len(cds), "bp")
            print("Protein length:", len(protein), "amino acids")
            print()
            print("Protein:")
            print(protein)

        else:
            print("No CDS annotation found.")

    else:
        print("No MANE Select transcript found.")

else:
    print("Gene not found.")