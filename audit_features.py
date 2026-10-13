"""Audit training metadata and compare predefined feature groups on identical folds."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from cross_validate_model import cross_validate
from extended_features import EXTENDED_FEATURES
from prepare_dataset import MODEL_FEATURES


def audit(source, output):
    source, output = Path(source), Path(output)
    if output.exists():
        raise ValueError('Output folder exists; choose a new folder')
    with source.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        if not {'review_status', 'clinical_significance', 'label', 'gene'} <= set(reader.fieldnames or []):
            raise ValueError('Missing audit metadata')
        rows = list(reader)
    reviews, classifications = defaultdict(Counter), defaultdict(Counter)
    for row in rows:
        reviews[row['review_status']][row['label']] += 1
        classifications[row['clinical_significance']][row['label']] += 1
    groups = {
        'original': MODEL_FEATURES,
        'original_plus_absolute': EXTENDED_FEATURES[:-1],
        'original_plus_blosum': MODEL_FEATURES + ['blosum62_score'],
        'all_six': EXTENDED_FEATURES,
    }
    output.mkdir(parents=True)
    comparisons = {}
    for name, features in groups.items():
        print(f'Comparing {name}...', flush=True)
        report = cross_validate(source, output / f'{name}.json', extended=True, feature_subset=features)
        comparisons[name] = report['aggregate']['logistic_regression']
    report = dict(training_rows=len(rows), genes=len({r['gene'] for r in rows}),
                  labels=dict(Counter(r['label'] for r in rows)),
                  review_status_by_label={k: dict(v) for k, v in sorted(reviews.items())},
                  clinical_significance_by_label={k: dict(v) for k, v in sorted(classifications.items())},
                  comparisons=comparisons,
                  note='All comparisons use identical training-only folds and fixed model settings. Review status describes evidence, not verified label correctness. No model selected or test evaluation performed.')
    (output / 'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('training_csv', type=Path)
    parser.add_argument('output_directory', type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(audit(args.training_csv, args.output_directory), indent=2))
    except ValueError as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
