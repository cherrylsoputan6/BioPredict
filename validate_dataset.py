"""Validate a small, ordered ClinVar candidate batch; cache MANE GenBank records."""
import argparse
import csv
import json
import re
import socket
from collections import Counter
from itertools import islice
from pathlib import Path

from Bio import SeqIO, Entrez
from Bio.SeqFeature import ExactPosition
from analysis import AMINO_ACID_PROPERTIES, extract_variant_features, validate_variant
from ncbi import find_mane_transcript, download_genbank_record

FEATURES = ['protein_length', 'relative_position', 'hydrophobicity_change', 'weight_change']


def matches_gene(feature, gene):
    names = set(feature.qualifiers.get('gene', []))
    for value in feature.qualifiers.get('gene_synonym', []):
        names.update(name.strip() for name in value.split(';'))
    return gene in names


def checked_protein(record, gene, accession):
    if record.id != accession:
        raise ValueError('Downloaded accession version differs from candidate')
    if record.annotations.get('organism') != 'Homo sapiens':
        raise ValueError('Record is not human')
    cds = [f for f in record.features if f.type == 'CDS' and matches_gene(f, gene)]
    if len(cds) != 1:
        raise ValueError('Expected one CDS for the candidate gene')
    feature = cds[0]
    if not all(isinstance(p.start, ExactPosition) and isinstance(p.end, ExactPosition)
               for p in feature.location.parts):
        raise ValueError('Partial CDS')
    if feature.qualifiers.get('codon_start', ['1']) != ['1'] or 'transl_except' in feature.qualifiers:
        raise ValueError('Unsupported CDS translation')
    table = int(feature.qualifiers.get('transl_table', ['1'])[0])
    protein = str(feature.extract(record.seq).translate(table=table, cds=True))
    if feature.qualifiers.get('translation') != [protein]:
        raise ValueError('Translation does not match annotation')
    return protein


class SequenceCache:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.proteins = {}
        self.mane = {}

    def get(self, gene, accession):
        if not re.fullmatch(r'[A-Za-z0-9_-]+', gene) or not re.fullmatch(r'NM_\d+\.\d+', accession):
            raise ValueError('Invalid gene or accession')
        if gene not in self.mane:
            print(f'Looking up MANE transcript for {gene}...', flush=True)
            self.mane[gene] = find_mane_transcript(gene)
        if self.mane[gene] != accession:
            return None
        key = (gene, accession)
        if key not in self.proteins:
            path = self.directory / f'{accession}.gb'
            if path.exists():
                record = SeqIO.read(path, 'genbank')
            else:
                print(f'Downloading {accession}...', flush=True)
                record = download_genbank_record(accession)
                protein = checked_protein(record, gene, accession)
                SeqIO.write(record, path, 'genbank')
                self.proteins[key] = protein
            self.proteins[key] = checked_protein(record, gene, accession)
        return self.proteins[key]


def validate_batch(source, output, cache, limit=100):
    if limit < 1:
        raise ValueError('Limit must be positive')
    source, output = Path(source), Path(output)
    if output.exists():
        raise ValueError('Output exists; choose a new filename')
    counts = Counter()
    with source.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        required = {'gene', 'transcript', 'variant', 'ref', 'alt', 'position', 'label'}
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Missing candidate CSV columns')
        if set(FEATURES + ['validation_status', 'validation_detail']) & set(reader.fieldnames):
            raise ValueError('Use the original candidate CSV as input')
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x', encoding='utf-8', newline='') as dest:
            writer = csv.DictWriter(dest, fieldnames=reader.fieldnames + FEATURES +
                                    ['validation_status', 'validation_detail'])
            writer.writeheader()
            for row in islice(reader, limit):
                counts['processed'] += 1
                print(f'Checking {counts["processed"]}/{limit}: {row["gene"]} {row["variant"]}', flush=True)
                status, detail = 'invalid_candidate', ''
                try:
                    position = int(row['position'])
                    ref, alt = row['ref'], row['alt']
                    if (position < 1 or ref not in AMINO_ACID_PROPERTIES or
                            alt not in AMINO_ACID_PROPERTIES or ref == alt or
                            row['variant'] != f'{ref}{position}{alt}' or row['label'] not in ('0', '1')):
                        raise ValueError('Invalid substitution or label')
                    status = 'sequence_error'
                    protein = cache.get(row['gene'], row['transcript'])
                    if protein is None:
                        status = 'not_current_mane'
                    else:
                        variant = {'original': ref, 'new': alt, 'position': position}
                        if validate_variant(protein, variant):
                            values = extract_variant_features(protein, variant)
                            row.update({name: values[name] for name in FEATURES})
                            status = 'valid'
                        else:
                            status = 'reference_mismatch'
                except (ValueError, OSError, RuntimeError) as error:
                    detail = str(error)
                row.update(validation_status=status, validation_detail=detail)
                writer.writerow(row)
                dest.flush()
                counts[status] += 1
                print(f'  {status}' + (f': {detail}' if detail else ''), flush=True)
    return dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--limit', type=int, default=100)
    parser.add_argument('--timeout', type=float, default=30, help='Socket timeout in seconds')
    parser.add_argument('--cache', type=Path, default=Path('data/sequence_cache'))
    args = parser.parse_args()
    if args.limit < 1 or not 0 < args.timeout <= 300 or args.output.exists():
        parser.error('Use a positive limit, timeout between 0 and 300 seconds, and a new output filename')
    previous_timeout, previous_tries = socket.getdefaulttimeout(), Entrez.max_tries
    socket.setdefaulttimeout(args.timeout)
    Entrez.max_tries = 1
    try:
        print(json.dumps(validate_batch(args.source, args.output, SequenceCache(args.cache), args.limit), indent=2))
    finally:
        socket.setdefaulttimeout(previous_timeout)
        Entrez.max_tries = previous_tries


if __name__ == '__main__':
    main()
