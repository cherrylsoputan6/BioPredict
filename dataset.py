"""Build labeled missense candidates from the official human ClinVar summary.

Candidates still need validation against their exact transcript before training.
This module uses only Python's standard library and makes no network requests.
"""

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

AMINO_ACIDS = dict(zip(
    'Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val'.split(),
    'ARNDCQEGHILKMFPSTWYV',
))
LABELS = {
    'Benign': 0, 'Likely benign': 0, 'Benign/Likely benign': 0,
    'Pathogenic': 1, 'Likely pathogenic': 1, 'Pathogenic/Likely pathogenic': 1,
}
NAME_PATTERN = re.compile(
    r'^(NM_\d+\.\d+)\(([^()]+)\):c\.[^\s]+ '
    r'\(p\.([A-Z][a-z]{2})([1-9]\d*)([A-Z][a-z]{2})\)$'
)
FIELDS = [
    'variation_id', 'allele_id', 'gene', 'transcript', 'variant',
    'position', 'ref', 'alt', 'label', 'clinical_significance',
    'review_status', 'rcv_accession', 'assembly', 'last_evaluated', 'clinvar_name',
]
REQUIRED = {
    '#AlleleID', 'VariationID', 'GeneSymbol', 'Name', 'Assembly',
    'ClinicalSignificance', 'ReviewStatus', 'RCVaccession', 'LastEvaluated',
}


def candidate_from_row(row):
    """Return (candidate, exclusion reason), keeping only unambiguous changes."""
    if row['Assembly'] != 'GRCh38':
        return None, 'other_assembly'
    label = LABELS.get(row['ClinicalSignificance'])
    if label is None:
        return None, 'excluded_classification'
    match = NAME_PATTERN.fullmatch(row['Name'])
    if not match:
        return None, 'unsupported_name'
    transcript, gene, ref_three, position, alt_three = match.groups()
    ref, alt = AMINO_ACIDS.get(ref_three), AMINO_ACIDS.get(alt_three)
    if ref is None or alt is None or ref == alt:
        return None, 'not_missense'
    if gene != row['GeneSymbol'] or gene in ('', '-') or ';' in gene:
        return None, 'ambiguous_gene'
    if not row['VariationID'].isdigit() or not row['#AlleleID'].isdigit():
        return None, 'missing_identifier'
    return dict(zip(FIELDS, [
        row['VariationID'], row['#AlleleID'], gene, transcript,
        f'{ref}{position}{alt}', int(position), ref, alt, label,
        row['ClinicalSignificance'], row['ReviewStatus'], row['RCVaccession'],
        row['Assembly'], row['LastEvaluated'], row['Name'],
    ])), None


def build_dataset(source, output):
    """Stream .txt or .txt.gz into CSV; report every exclusion and deduplicate IDs."""
    source, output = Path(source), Path(output)
    if source.resolve() == output.resolve():
        raise ValueError('Input and output must be different files.')
    opener = gzip.open if source.suffix == '.gz' else open
    counts, seen = Counter(), set()
    with opener(source, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f'Missing ClinVar columns: {", ".join(sorted(missing))}')
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('w', encoding='utf-8', newline='') as destination:
            writer = csv.DictWriter(destination, fieldnames=FIELDS)
            writer.writeheader()
            for row in reader:
                counts['rows_read'] += 1
                if None in row or any(row.get(key) is None for key in REQUIRED):
                    counts['malformed_row'] += 1
                    continue
                candidate, reason = candidate_from_row(row)
                if reason:
                    counts[reason] += 1
                    continue
                key = candidate['variation_id']
                if key in seen:
                    counts['duplicate_variation_id'] += 1
                    continue
                seen.add(key)
                writer.writerow(candidate)
                counts['retained'] += 1
                counts[f'label_{candidate["label"]}'] += 1
    return dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path, help='ClinVar variant_summary.txt[.gz]')
    parser.add_argument('output', type=Path, help='Candidate CSV destination')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output already exists; choose a new filename.')
    print(json.dumps(build_dataset(args.source, args.output), indent=2))


if __name__ == '__main__':
    main()
