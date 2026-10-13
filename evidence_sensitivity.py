"""Training-only sensitivity analysis excluding records without assertion criteria."""
import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from cross_validate_model import cross_validate


def filter_rows(rows):
    if any(not r.get('review_status', '').strip() for r in rows):
        raise ValueError('Missing review status; inspect records before filtering')
    kept = [r for r in rows if r['review_status'].strip().lower() != 'no assertion criteria provided']
    return kept


def analyze(source, destination):
    source, destination = Path(source), Path(destination)
    if destination.exists():
        raise ValueError('Destination exists; choose a new folder')
    with source.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if not {'review_status', 'label', 'gene'} <= set(fields):
            raise ValueError('Missing metadata columns')
        original = list(reader)
    filtered = filter_rows(original)
    if {r['label'] for r in filtered} != {'0', '1'}:
        raise ValueError('Both classes must remain after filtering')
    destination.mkdir(parents=True)
    filtered_path = destination / 'training_filtered.csv'
    with filtered_path.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(filtered)
    comparisons = {}
    for name, extended in [('original_features', False), ('extended_features', True)]:
        print(f'Checking {name} on filtered training rows...', flush=True)
        result = cross_validate(filtered_path, destination / f'{name}.json', extended=extended)
        comparisons[name] = result['aggregate']
    with source.open('rb') as handle:
        source_hash = hashlib.file_digest(handle, 'sha256').hexdigest()
    summary = dict(original_rows=len(original), retained_rows=len(filtered),
                   excluded_rows=len(original) - len(filtered),
                   labels_before=dict(Counter(r['label'] for r in original)),
                   labels_after=dict(Counter(r['label'] for r in filtered)),
                   retained_genes=len({r['gene'] for r in filtered}),
                   source_sha256=source_hash, results=comparisons,
                   note='Sensitivity analysis only. Filtering changes class and gene composition and CV folds; differences cannot be attributed solely to label quality. No test data used and no saved model changed.')
    (destination / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('training_csv', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(analyze(args.training_csv, args.destination), indent=2))
    except ValueError as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
