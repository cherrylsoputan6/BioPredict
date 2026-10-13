"""Compare six versus eight features on exactly the same training cohort/folds."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from cross_validate_model import cross_validate


def load(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames or [], list(reader)


def identity(row):
    return tuple(row[f] for f in ('variation_id', 'gene', 'transcript', 'variant'))


def matched_cohort(training_rows, ortholog_rows):
    original = {identity(r): r for r in training_rows}
    observed = {identity(r): r for r in ortholog_rows}
    if (len(original) != len(training_rows) or len(observed) != len(ortholog_rows) or
            set(original) != set(observed)):
        raise ValueError('Ortholog data must cover exactly the original training candidate identities')
    for row in ortholog_rows:
        if any(row.get(f) != value for f, value in original[identity(row)].items()):
            raise ValueError('Original training values changed in the ortholog CSV')
    return [observed[identity(row)] for row in training_rows
            if observed[identity(row)].get('ortholog_status') == 'available']


def compare(training_csv, ortholog_csv, output):
    output = Path(output)
    if output.exists():
        raise ValueError('Output directory exists; choose a new folder')
    _, training_rows = load(training_csv)
    fields, ortholog_rows = load(ortholog_csv)
    cohort = matched_cohort(training_rows, ortholog_rows)
    if not cohort:
        raise ValueError('No available conservation rows')
    output.mkdir(parents=True)
    cohort_path = output / 'matched_training.csv'
    with cohort_path.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(cohort)
    reports = {}
    for name, conservation in [('six_features', False), ('with_conservation', True)]:
        print(f'Comparing {name}...', flush=True)
        reports[name] = cross_validate(cohort_path, output / f'{name}.json', extended=True, conservation=conservation)
    # Cross-validation depends only on row order, labels, gene IDs, seed and fold count.
    assert reports['six_features']['source_sha256'] == reports['with_conservation']['source_sha256']
    summary = dict(input_rows=len(training_rows), compared_rows=len(cohort),
                   compared_labels=dict(Counter(r['label'] for r in cohort)),
                   coverage_by_status_and_label=dict(Counter(f'{r["ortholog_status"]}:label_{r["label"]}' for r in ortholog_rows)),
                   identical_source_sha256=reports['six_features']['source_sha256'],
                   features={k: v['features'] for k, v in reports.items()},
                   results={k: v['aggregate']['logistic_regression'] for k, v in reports.items()},
                   note='Training-only, identical eligible rows and folds. Coverage may differ by label; results apply to this covered subset. Rows with connection errors or insufficient alignment coverage are excluded. No test evaluation or model replacement.')
    (output / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('training_csv', type=Path)
    parser.add_argument('ortholog_csv', type=Path)
    parser.add_argument('output_directory', type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(compare(args.training_csv, args.ortholog_csv, args.output_directory), indent=2))
    except ValueError as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
