"""Out-of-fold error and calibration diagnostics; no model or threshold tuning."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import numpy as np
from sklearn.metrics import brier_score_loss, log_loss
from sklearn.model_selection import StratifiedGroupKFold
from train_model import read_split, fit_baseline, metrics
from evaluate_conservation import add_conservation, checksum


def score_bins(labels, scores, count=5):
    bins = []
    for index in range(count):
        low, high = index / count, (index + 1) / count
        selected = [i for i, s in enumerate(scores) if low <= s < high or (index == count - 1 and s == 1)]
        bins.append(dict(lower=low, upper=high, rows=len(selected),
                         mean_score=float(np.mean([scores[i] for i in selected])) if selected else None,
                         observed_label_1_fraction=float(np.mean([labels[i] for i in selected])) if selected else None))
    return bins


def diagnose(source, output, config_path):
    output = Path(output)
    if output.exists():
        raise ValueError('Output exists; choose a new folder')
    config = json.loads(Path(config_path).read_text(encoding='utf-8'))
    if checksum(source) != config['training_cohort_sha256']:
        raise ValueError('Use the frozen matched training cohort')
    rows, six, labels = read_split(source, extended=True)
    values = np.asarray(add_conservation(rows, six))
    y = np.asarray(labels)
    groups = np.asarray([r['gene'] for r in rows])
    scores, reference_scores = np.empty(len(rows)), np.empty(len(rows))
    coverage = np.zeros(len(rows), dtype=int)
    splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
    for number, (train, heldout) in enumerate(splitter.split(values, y, groups), 1):
        if set(groups[train]) & set(groups[heldout]) or set(y[train]) != {0, 1}:
            raise ValueError('Invalid gene-held-out fold')
        model = fit_baseline(values[train], y[train])
        scores[heldout] = model.predict_proba(values[heldout])[:, 1]
        reference_scores[heldout] = np.mean(y[train])
        coverage[heldout] += 1
        print(f'Diagnostic fold {number}/5 complete', flush=True)
    if not np.all(coverage == 1):
        raise ValueError('Every row must be held out once')
    predicted = (scores >= 0.5).astype(int)
    categories = []
    for actual, prediction in zip(y, predicted):
        categories.append('false_positive' if actual == 0 and prediction == 1 else
                          'false_negative' if actual == 1 and prediction == 0 else 'correct')
    by_review = {}
    for row, category in zip(rows, categories):
        status = row['review_status']
        if status not in by_review:
            by_review[status] = Counter()
        by_review[status]['rows'] += 1
        by_review[status][f'label_{row["label"]}'] += 1
        by_review[status][category] += 1
    report = dict(rows=len(rows), gene_overlap_per_fold=0,
                  metrics=metrics(y, predicted, scores),
                  brier_score=float(brier_score_loss(y, scores)),
                  prevalence_baseline_brier=float(brier_score_loss(y, reference_scores)),
                  log_loss=float(log_loss(y, scores, labels=[0, 1])),
                  score_bins=score_bins(y, scores), error_counts=dict(Counter(categories)),
                  errors_by_review_status={k: dict(v) for k, v in sorted(by_review.items())},
                  source_sha256=checksum(source), threshold=0.5,
                  note='Training-only out-of-fold diagnostics. Bin frequencies describe this selected covered dataset, not disease risk. No calibrator fitted, threshold tuned, or test set used.')
    output.mkdir(parents=True)
    (output / 'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    with (output / 'out_of_fold_predictions.csv').open('x', encoding='utf-8', newline='') as handle:
        fields = ['gene', 'transcript', 'variant', 'label', 'review_status', 'model_label', 'score_label_1', 'error_category']
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row, prediction, score, category in zip(rows, predicted, scores, categories):
            writer.writerow(dict({f: row[f] for f in fields[:5]}, model_label=int(prediction),
                                 score_label_1=float(score), error_category=category))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('training_csv', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--config', type=Path, default=Path('model_config.json'))
    args = parser.parse_args()
    print(json.dumps(diagnose(args.training_csv, args.output, args.config), indent=2))


if __name__ == '__main__':
    main()
