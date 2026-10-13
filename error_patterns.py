"""Descriptive feature-group error rates from existing training-only predictions."""
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from extended_features import extended_values


def groups(row):
    ref, alt = float(row['ortholog_ref_agreement']), float(row['ortholog_alt_agreement'])
    score = extended_values(row)[-1]
    position = float(row['relative_position'])
    return {
        'reference_agreement': 'all observed species agree' if ref == 1 else 'some or none agree',
        'alternate_agreement': 'alternate observed' if alt > 0 else 'alternate not observed',
        'blosum62': 'negative' if score < 0 else 'zero' if score == 0 else 'positive',
        'relative_position': 'first third' if position <= 1 / 3 else 'middle third' if position <= 2 / 3 else 'last third',
    }


def summarize(rows):
    counts = defaultdict(Counter)
    for row in rows:
        for dimension, group in groups(row).items():
            bucket = counts[(dimension, group)]
            label, prediction = row['label'], row['model_label']
            bucket['rows'] += 1
            if label == '0':
                bucket['benign'] += 1
                bucket['false_positive'] += prediction == '1'
            else:
                bucket['pathogenic'] += 1
                bucket['false_negative'] += prediction == '0'
    output = []
    for (dimension, group), count in sorted(counts.items()):
        output.append(dict(dimension=dimension, group=group, **dict(count),
                           false_positive_rate=count['false_positive'] / count['benign'] if count['benign'] else None,
                           false_negative_rate=count['false_negative'] / count['pathogenic'] if count['pathogenic'] else None))
    return output


def analyze(training, predictions, output):
    output = Path(output)
    if output.exists():
        raise ValueError('Output exists; choose a new filename')
    def load(path):
        with Path(path).open(encoding='utf-8', newline='') as handle:
            return list(csv.DictReader(handle))
    def key(row):
        return tuple(row[f] for f in ('gene', 'transcript', 'variant'))
    rows, predicted = load(training), load(predictions)
    index = {key(r): r for r in predicted}
    if len(index) != len(predicted) or len(rows) != len(predicted) or {key(r) for r in rows} != set(index):
        raise ValueError('Predictions must cover the same unique training variants')
    joined = []
    for row in rows:
        result = index[key(row)]
        if result['label'] != row['label'] or result['model_label'] not in ('0', '1'):
            raise ValueError('Prediction labels differ from the source cohort')
        joined.append(dict(row, model_label=result['model_label']))
    report = dict(rows=len(joined), groups=summarize(joined),
                  scope='Exploratory training-only feature-group rates. Groups overlap across dimensions. Small groups are unstable; associations are not causal. No threshold or feature changes made.')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('training_csv', type=Path)
    parser.add_argument('predictions_csv', type=Path)
    parser.add_argument('output_json', type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze(args.training_csv, args.predictions_csv, args.output_json), indent=2))


if __name__ == '__main__':
    main()
