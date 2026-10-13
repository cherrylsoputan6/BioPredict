"""Conservative grouping by symbol, stable NCBI GeneID and exact protein identity."""
from functools import lru_cache
from pathlib import Path
import hashlib
from Bio import SeqIO
from validate_dataset import matches_gene, checked_protein


def connected_groups(symbols, identities):
    parent = {s: s for s in symbols}
    def root(s):
        while parent[s] != s:
            parent[s] = parent[parent[s]]
            s = parent[s]
        return s
    owners = {}
    for symbol, keys in zip(symbols, identities):
        for key in keys:
            if key in owners:
                a, b = root(symbol), root(owners[key])
                parent[max(a, b)] = min(a, b)
            else:
                owners[key] = symbol
    return [root(s) for s in symbols]


@lru_cache(maxsize=1000)
def record_identities(gene, transcript):
    record = SeqIO.read(Path('data/sequence_cache') / f'{transcript}.gb', 'genbank')
    protein = str(checked_protein(record, gene, transcript))
    gene_ids = {xref for f in record.features if f.type in ('gene', 'CDS') and matches_gene(f, gene)
                for xref in f.qualifiers.get('db_xref', []) if xref.startswith('GeneID:')}
    if len(gene_ids) != 1:
        raise ValueError('Expected one stable NCBI GeneID for grouping')
    return list(gene_ids) + ['protein:' + hashlib.sha256(protein.encode()).hexdigest()]


def biological_groups(rows):
    return connected_groups([r['gene'] for r in rows],
                            [record_identities(r['gene'], r['transcript']) for r in rows])
