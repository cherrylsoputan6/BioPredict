"""Fit a fixed pilot baseline and evaluate once on the existing gene-held-out test."""
import argparse
import csv
import hashlib
import json
import math
import platform
from pathlib import Path

import joblib
import sklearn
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, confusion_matrix,
                             precision_score, recall_score, f1_score, roc_auc_score)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from prepare_dataset import MODEL_FEATURES


def read_split(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        required = set(MODEL_FEATURES + ['gene', 'transcript', 'variant', 'label', 'validation_status'])
        if not required <= set(reader.fieldnames or []):
            raise ValueError('Missing split columns')
        rows = list(reader)
    if not rows:
        raise ValueError('Empty split')
    x, y, keys = [], [], set()
    for row in rows:
        values = [float(row[f]) for f in MODEL_FEATURES]
        key = tuple(row[f] for f in ('gene', 'transcript', 'variant'))
        if (row['validation_status'] != 'valid' or row['label'] not in ('0', '1') or
                not all(key) or not all(math.isfinite(v) for v in values) or key in keys):
            raise ValueError('Invalid or duplicate split row')
        keys.add(key)
        x.append(values)
        y.append(int(row['label']))
    if set(y) != {0, 1}:
        raise ValueError('Both labels are required in each split')
    return rows, x, y


def fit_baseline(x, y):
    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('classifier', LogisticRegression(C=1.0, solver='lbfgs', max_iter=1000)),
    ])
    pipeline.fit(x, y)
    return pipeline


def metrics(y, predicted, scores):
    return dict(accuracy=float(accuracy_score(y, predicted)),
                balanced_accuracy=float(balanced_accuracy_score(y, predicted)),
                precision_label_1=float(precision_score(y, predicted, zero_division=0)),
                recall_label_1=float(recall_score(y, predicted, zero_division=0)),
                f1_label_1=float(f1_score(y, predicted, zero_division=0)),
                roc_auc=float(roc_auc_score(y, scores)),
                confusion_matrix=confusion_matrix(y, predicted, labels=[0, 1]).tolist())


def train(split_directory, output):
    split_directory, output = Path(split_directory), Path(output)
    if output.exists():
        raise ValueError('Output folder already exists; choose a new folder')
    train_rows, x_train, y_train = read_split(split_directory / 'train.csv')
    test_rows, x_test, y_test = read_split(split_directory / 'test.csv')
    if {r['gene'] for r in train_rows} & {r['gene'] for r in test_rows}:
        raise ValueError('Gene overlap between training and test sets')
    model = fit_baseline(x_train, y_train)
    dummy = DummyClassifier(strategy='most_frequent').fit(x_train, y_train)
    predicted, scores = model.predict(x_test), model.predict_proba(x_test)[:, 1]
    hashes = {}
    for name in ('train', 'test'):
        with (split_directory / f'{name}.csv').open('rb') as handle:
            hashes[name] = hashlib.file_digest(handle, 'sha256').hexdigest()
    report = dict(
        train_rows=len(train_rows), test_rows=len(test_rows), gene_overlap=0,
        features=MODEL_FEATURES, threshold=0.5,
        logistic_regression=metrics(y_test, predicted, scores),
        majority_baseline=metrics(y_test, dummy.predict(x_test), dummy.predict_proba(x_test)[:, 1]),
        confusion_matrix_order='rows: true 0,1; columns: predicted 0,1',
        scaling_mean=model.named_steps['scaler'].mean_.tolist(),
        coefficients=model.named_steps['classifier'].coef_[0].tolist(),
        intercept=float(model.named_steps['classifier'].intercept_[0]),
        source_sha256=hashes, python_version=platform.python_version(),
        sklearn_version=sklearn.__version__,
        limitation='Small balanced pilot; probabilities are uncalibrated, not clinical risk. No test-based tuning.',
    )
    output.mkdir(parents=True)
    joblib.dump(dict(pipeline=model, features=MODEL_FEATURES, sklearn_version=sklearn.__version__), output / 'model.joblib')
    (output / 'metrics.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    with (output / 'test_predictions.csv').open('x', encoding='utf-8', newline='') as handle:
        fields = ['gene', 'transcript', 'variant', 'true_label', 'predicted_label', 'score_label_1']
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row, label, prediction, score in zip(test_rows, y_test, predicted, scores):
            writer.writerow(dict(gene=row['gene'], transcript=row['transcript'], variant=row['variant'],
                                 true_label=label, predicted_label=int(prediction), score_label_1=float(score)))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('split_directory', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(train(args.split_directory, args.output), indent=2))
    except ValueError as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
