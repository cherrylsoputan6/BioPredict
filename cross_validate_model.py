"""Evaluate the fixed baseline on gene-held-out folds of training data only."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import sklearn
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import StratifiedGroupKFold
from train_model import read_split, fit_baseline, metrics
from prepare_dataset import MODEL_FEATURES


def cross_validate(source, output, folds=5, seed=42, extended=False, feature_subset=None, conservation=False, context=False):
    source, output = Path(source), Path(output)
    if output.exists():
        raise ValueError('Output already exists; choose a new filename')
    rows, x, y = read_split(source, extended=extended)
    feature_names = MODEL_FEATURES
    if extended:
        from extended_features import EXTENDED_FEATURES
        feature_names = EXTENDED_FEATURES
    if conservation:
        extra_names = ['ortholog_ref_agreement', 'ortholog_alt_agreement']
        feature_names = list(feature_names) + extra_names
        for row, values in zip(rows, x):
            if row.get('ortholog_status') != 'available' or int(row.get('ortholog_observations', 0)) < 3:
                raise ValueError('Conservation CV requires available observations from at least three species')
            extras = [float(row[name]) for name in extra_names]
            if not all(np.isfinite(v) and 0 <= v <= 1 for v in extras):
                raise ValueError('Invalid ortholog agreement fractions')
            values.extend(extras)
    if context:
        from sequence_context import CONTEXT_FEATURES, row_context
        feature_names = list(feature_names) + CONTEXT_FEATURES
        for row, values in zip(rows, x):
            values.extend(row_context(row))
    if feature_subset is not None:
        if not feature_subset or len(set(feature_subset)) != len(feature_subset) or not set(feature_subset) <= set(feature_names):
            raise ValueError('Feature subset must contain unique allowed feature names')
        indices = [feature_names.index(name) for name in feature_subset]
        x = [[values[index] for index in indices] for values in x]
        feature_names = list(feature_subset)
    x, y = np.asarray(x), np.asarray(y)
    groups = np.asarray([r['gene'] for r in rows])
    if folds < 2 or len(set(groups)) < folds:
        raise ValueError('Need at least two folds and enough distinct genes')
    for label in (0, 1):
        if len(set(groups[y == label])) < folds:
            raise ValueError('Too few genes for this number of folds in one label class')
    splitter = StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=seed)
    fold_results = []
    predictions = np.empty(len(y), dtype=int)
    scores = np.empty(len(y), dtype=float)
    covered = np.zeros(len(y), dtype=int)
    for number, (train, validation) in enumerate(splitter.split(x, y, groups), 1):
        overlap = set(groups[train]) & set(groups[validation])
        if overlap or any(set(y[index]) != {0, 1} for index in (train, validation)):
            raise ValueError('Fold contains overlapping genes or a missing label class')
        model = fit_baseline(x[train], y[train])
        dummy = DummyClassifier(strategy='most_frequent').fit(x[train], y[train])
        predicted = model.predict(x[validation])
        probability = model.predict_proba(x[validation])[:, 1]
        predictions[validation], scores[validation] = predicted, probability
        covered[validation] += 1
        result = dict(fold=number, training_rows=len(train), validation_rows=len(validation),
                      gene_overlap=0,
                      validation_label_0=int(sum(y[validation] == 0)),
                      validation_label_1=int(sum(y[validation] == 1)),
                      logistic_regression=metrics(y[validation], predicted, probability),
                      majority_baseline=metrics(y[validation], dummy.predict(x[validation]),
                                               dummy.predict_proba(x[validation])[:, 1]))
        fold_results.append(result)
        print(f'Fold {number}/{folds}: accuracy {result["logistic_regression"]["accuracy"]:.3f}', flush=True)
    if not np.all(covered == 1):
        raise ValueError('Each row must be held out exactly once')
    aggregated = {}
    for model_name in ('logistic_regression', 'majority_baseline'):
        aggregated[model_name] = {}
        for measure in ('accuracy', 'balanced_accuracy', 'precision_label_1',
                        'recall_label_1', 'f1_label_1', 'roc_auc'):
            values = [f[model_name][measure] for f in fold_results]
            aggregated[model_name][measure] = dict(mean=float(np.mean(values)),
                                                  std=float(np.std(values, ddof=1)))
    with source.open('rb') as handle:
        checksum = hashlib.file_digest(handle, 'sha256').hexdigest()
    report = dict(rows=len(rows), genes=len(set(groups)), folds=folds, seed=seed,
                  source_sha256=checksum, sklearn_version=sklearn.__version__,
                  features=feature_names, fold_results=fold_results,
                  aggregate=aggregated,
                  out_of_fold=metrics(y, predictions, scores),
                  note='Training-only cross-validation. Fold standard deviations are not confidence intervals. No tuning performed.')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as handle:
        json.dump(report, handle, indent=2)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('training_csv', type=Path)
    parser.add_argument('output_json', type=Path)
    parser.add_argument('--folds', type=int, default=5)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--extended', action='store_true', help='Add absolute changes and BLOSUM62 score')
    parser.add_argument('--conservation', action='store_true', help='Use available ortholog agreement columns')
    parser.add_argument('--context', action='store_true', help='Add fixed five-residue local context from cached MANE protein')
    args = parser.parse_args()
    try:
        report = cross_validate(args.training_csv, args.output_json, args.folds, args.seed, args.extended,
                                conservation=args.conservation, context=args.context)
        print(json.dumps(report['aggregate'], indent=2))
    except ValueError as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
