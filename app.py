from ncbi import (
    search_gene,
    download_sequence,
    get_description
)

from analysis import analyze_sequence

gene = input("Enter gene: ")

ids = search_gene(gene)

print("\nSearch Results:")
print("-" * 50)

for gene_id in ids:
    print(get_description(gene_id))