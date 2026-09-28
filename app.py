from ncbi import (
    search_gene,
    find_mane_transcript,
    download_genbank_record
)

from analysis import (
    analyze_sequence,
    extract_cds,
    translate_cds
)


gene = input("Enter gene: ")

ids = search_gene(gene)

if ids:
    gene_id = ids[0]

    transcript = find_mane_transcript(gene)

    if transcript:
        record = download_genbank_record(transcript)

        sequence = str(record.seq)

        results = analyze_sequence(sequence)

        cds = extract_cds(record)
        protein = translate_cds(cds)

        print()
        print("BIOPREDICT ANALYSIS REPORT")
        print("--------------------------")
        print()
        print("Gene:", gene.upper())
        print("NCBI Gene ID:", gene_id)
        print("MANE Select Transcript:", transcript)

        print()
        print("SEQUENCE ANALYSIS")
        print("Transcript Length:", results["length"], "bp")
        print("GC Content:", f'{results["gc"] * 100:.2f}%')

        print()
        print("Nucleotide Counts")
        print("A:", results["counts"]["A"])
        print("T:", results["counts"]["T"])
        print("G:", results["counts"]["G"])
        print("C:", results["counts"]["C"])

        if cds:
            print()
            print("CODING SEQUENCE")
            print("CDS Length:", len(cds), "bp")

            print()
            print("PROTEIN")
            print("Protein Length:", len(protein), "amino acids")
            print("Protein Sequence:")
            print(protein)

        else:
            print()
            print("No CDS annotation found.")

    else:
        print("No MANE Select transcript found.")

else:
    print("Gene not found.")