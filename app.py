from ncbi import (
    search_gene,
    find_mane_transcript,
    download_sequence
)

from analysis import analyze_sequence


gene = input("Enter gene: ")

ids = search_gene(gene)

if ids:
    gene_id = ids[0]

    print()
    print("Gene:", gene.upper())
    print("NCBI Gene ID:", gene_id)

    print("Finding MANE Select transcript...")
    transcript = find_mane_transcript(gene)

    if transcript:
        print("MANE Select transcript:", transcript)

        print("Downloading sequence...")
        record = download_sequence(transcript)

        sequence = str(record.seq)

        print("Analyzing sequence...")
        results = analyze_sequence(sequence)

        print()
        print("DNA Analysis Report")
        print("-------------------")
        print("Transcript:", transcript)
        print("Length:", results["length"], "bp")
        print("GC content:", f'{results["gc"] * 100:.2f}%')
        print("A:", results["counts"]["A"])
        print("T:", results["counts"]["T"])
        print("G:", results["counts"]["G"])
        print("C:", results["counts"]["C"])

    else:
        print("No MANE Select transcript found.")

else:
    print("Gene not found.")