"""Evaluate the frozen configuration on covered rows of the existing test split."""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
import joblib
import sklearn
from sklearn.dummy import DummyClassifier
from compare_conservation import load, matched_cohort
from extended_features import EXTENDED_FEATURES
from train_model import read_split, fit_baseline, metrics

EXTRA = ['ortholog_ref_agreement', 'ortholog_alt_agreement']


def checksum(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def add_conservation(rows, values, minimum=3):
    result = []
    for row, base in zip(rows, values):
        if row['ortholog_status'] != 'available' or int(row['ortholog_observations']) < minimum:
            raise ValueError('Insufficient ortholog coverage')
        extras = [float(row[f]) for f in EXTRA]
        if not all(math.isfinite(v) and 0 <= v <= 1 for v in extras) or sum(extras) > 1 + 1e-9:
            raise ValueError('Invalid agreement fractions')
        result.append(list(base) + extras)
    return result


def evaluate(training, original_test, ortholog_test, output, config_path):
    output = Path(output)
    if output.exists():
        raise ValueError('Output exists; evaluation will not overwrite a previous run')
    config = json.loads(Path(config_path).read_text(encoding='utf-8'))
    if (config['features'] != EXTENDED_FEATURES + EXTRA or
            config['estimator'] != dict(type='LogisticRegression', C=1.0, solver='lbfgs', max_iter=1000) or
            config['threshold_label_1'] != 0.5 or config['conservation']['minimum_observations'] != 3):
        raise ValueError('Frozen settings differ from the implemented baseline')
    if checksum(training) != config['training_cohort_sha256']:
        raise ValueError('Training cohort differs from the frozen configuration')
    train_rows, x_train, y_train = read_split(training, extended=True)
    if len(train_rows) != config['training_rows']:
        raise ValueError('Training row count differs from frozen configuration')
    _, original = load(original_test)
    fields, observations = load(ortholog_test)
    eligible = matched_cohort(original, observations)
    if {r['gene'] for r in train_rows} & {r['gene'] for r in original}:
        raise ValueError('Training and test genes overlap')
    if not eligible:
        raise ValueError('No eligible test rows')
    output.mkdir(parents=True)
    test_path = output / 'eligible_test.csv'
    with test_path.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(eligible)
    test_rows, x_test, y_test = read_split(test_path, extended=True)
    x_train_eight = add_conservation(train_rows, x_train)
    x_test_eight = add_conservation(test_rows, x_test)
    # Fixed configuration and threshold; test labels enter only the metric calls.
    six = fit_baseline(x_train, y_train)
    eight = fit_baseline(x_train_eight, y_train)
    dummy = DummyClassifier(strategy='most_frequent').fit(x_train, y_train)
    predictions, results = {}, {}
    for name, model, values in [('six_features', six, x_test), ('eight_features', eight, x_test_eight), ('majority_baseline', dummy, x_test)]:
        predicted, scores = model.predict(values), model.predict_proba(values)[:, 1]
        predictions[name] = (predicted, scores)
        results[name] = metrics(y_test, predicted, scores)
    report = dict(train_rows=len(train_rows), test_total=len(original), test_evaluated=len(test_rows),
                  test_labels=dict(Counter(r['label'] for r in test_rows)),
                  coverage_by_status_and_label=dict(Counter(f'{r["ortholog_status"]}:label_{r["label"]}' for r in observations)),
                  results=results, features=config['features'], threshold=0.5, gene_overlap=0,
                  sklearn_version=sklearn.__version__, config_sha256=checksum(config_path),
                  training_sha256=checksum(training), original_test_sha256=checksum(original_test),
                  ortholog_test_sha256=checksum(ortholog_test),
                  confusion_matrix_order='rows: true 0,1; columns: predicted 0,1',
                  limitations='Covered test subset only; differential coverage by class. This split was previously evaluated for the original baseline. No tuning performed; probabilities are uncalibrated.')
    (output / 'metrics.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    joblib.dump(dict(pipeline=eight, features=config['features'], configuration=config, sklearn_version=sklearn.__version__), output / 'model.joblib')
    with (output / 'predictions.csv').open('x', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['gene', 'transcript', 'variant', 'true_label', 'six_predicted', 'eight_predicted', 'eight_score_label_1'])
        for i, row in enumerate(test_rows):
            writer.writerow([row['gene'], row['transcript'], row['variant'], y_test[i],
                             int(predictions['six_features'][0][i]), int(predictions['eight_features'][0][i]),
                             float(predictions['eight_features'][1][i])])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('training_csv', type=Path)
    parser.add_argument('original_test_csv', type=Path)
    parser.add_argument('ortholog_test_csv', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--config', type=Path, default=Path('model_config.json'))
    args = parser.parse_args()
    print(json.dumps(evaluate(args.training_csv, args.original_test_csv, args.ortholog_test_csv, args.output, args.config), indent=2))


if __name__ == '__main__':
    main()
