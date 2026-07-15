from ncbi import search_gene

gene = input("Enter gene: ")

ids = search_gene(gene)

print(ids)