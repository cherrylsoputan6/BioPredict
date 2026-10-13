"""Reproducible training-only scientific audit; never loads test rows or replaces models."""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import sklearn
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score, brier_score_loss
from sklearn.model_selection import StratifiedGroupKFold
from train_model import read_split, fit_baseline, metrics
from evaluate_conservation import add_conservation, checksum
from sequence_context import row_context, cached_protein
from training_diagnostics import score_bins
from validation_groups import biological_groups


def paired_gene_bootstrap(y, base, context, genes, repetitions=2000, seed=2026):
    """Conditional paired cluster bootstrap of fixed OOF predictions, not model refits."""
    rng = np.random.default_rng(seed)
    unique = np.unique(genes)
    members = [np.flatnonzero(genes == gene) for gene in unique]
    measures = {
        'accuracy': lambda a, s: accuracy_score(a, s >= 0.5),
        'balanced_accuracy': lambda a, s: balanced_accuracy_score(a, s >= 0.5),
        'roc_auc': roc_auc_score,
        'brier': brier_score_loss,
    }
    draws = {name: [] for name in measures}
    for _ in range(repetitions):
        indices = np.concatenate([members[i] for i in rng.integers(0, len(unique), len(unique))])
        if len(np.unique(y[indices])) < 2:
            continue
        for name, measure in measures.items():
            draws[name].append(float(measure(y[indices], context[indices]) - measure(y[indices], base[indices])))
    return {name: dict(context_minus_baseline=float(measures[name](y, context) - measures[name](y, base)),
                       conditional_95_percentile_interval=np.quantile(values, [0.025, 0.975]).tolist(),
                       bootstrap_draws=len(values)) for name, values in draws.items()}


def validate(source, output, seeds=(42, 17, 123), bootstrap=2000):
    output = Path(output)
    if output.exists():
        raise ValueError('Output folder exists')
    config = json.loads(Path('model_config.json').read_text(encoding='utf-8'))
    if checksum(source) != config['training_cohort_sha256']:
        raise ValueError('Use frozen training cohort; no test input accepted')
    rows, six, labels = read_split(source, extended=True)
    x = np.asarray(add_conservation(rows, six))
    context_x = np.column_stack((x, np.asarray([row_context(r) for r in rows])))
    y = np.asarray(labels)
    symbols = np.asarray([r['gene'] for r in rows])
    genes = np.asarray(biological_groups(rows))
    proteins = [str(cached_protein(r['gene'], r['transcript'])) for r in rows]
    hashes = [hashlib.sha256(p.encode()).hexdigest() for p in proteins]
    sequence_genes = defaultdict(set)
    for gene, digest in zip(symbols, hashes):
        sequence_genes[digest].add(gene)
    duplicate_sequences = {k: sorted(v) for k, v in sequence_genes.items() if len(v) > 1}
    fold_reports, prediction_records, first_scores = [], [], None
    for seed in seeds:
        scores = {name: np.empty(len(y)) for name in ('baseline', 'context')}
        seen = np.zeros(len(y), dtype=int)
        for fold, (train, heldout) in enumerate(StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed).split(x, y, genes), 1):
            if set(genes[train]) & set(genes[heldout]) or any(set(y[i]) != {0, 1} for i in (train, heldout)):
                raise ValueError('Invalid grouped fold')
            report = dict(seed=seed, fold=fold, train_rows=len(train), validation_rows=len(heldout),
                          gene_overlap=0, identical_protein_overlap=len(set(np.asarray(hashes)[train]) & set(np.asarray(hashes)[heldout])),
                          heldout_gene_ids=sorted(set(genes[heldout])))
            for name, features in [('baseline', x), ('context', context_x)]:
                model = fit_baseline(features[train], y[train])
                values = model.predict_proba(features[heldout])[:, 1]
                scores[name][heldout] = values
                report[name] = metrics(y[heldout], values >= 0.5, values)
                report[name]['brier'] = float(brier_score_loss(y[heldout], values))
            report['paired_delta'] = {m: report['context'][m] - report['baseline'][m] for m in ('accuracy', 'balanced_accuracy', 'roc_auc', 'brier')}
            fold_reports.append(report)
            seen[heldout] += 1
            print(f'Seed {seed}, fold {fold}/5 complete', flush=True)
        if not np.all(seen == 1):
            raise ValueError('OOF coverage invalid')
        if first_scores is None:
            first_scores = scores
        for i, row in enumerate(rows):
            prediction_records.append(dict(seed=seed, gene=row['gene'], transcript=row['transcript'], variant=row['variant'],
                                           label=int(y[i]), review_status=row['review_status'],
                                           baseline_score=float(scores['baseline'][i]), context_score=float(scores['context'][i])))
    by_gene, by_review = [], []
    for field, destination in [('gene', by_gene), ('review_status', by_review)]:
        for name in sorted({r[field] for r in rows}):
            indices = np.array([i for i, r in enumerate(rows) if r[field] == name])
            truth, score = y[indices], first_scores['baseline'][indices]
            pred = score >= 0.5
            destination.append(dict(group=name, rows=len(indices), benign=int(sum(truth == 0)), pathogenic=int(sum(truth == 1)),
                                    false_positive=int(sum((truth == 0) & pred)), false_negative=int(sum((truth == 1) & ~pred)),
                                    accuracy=float(accuracy_score(truth, pred)),
                                    note='Too few rows/classes for reliable gene-level AUC' if len(indices) < 10 or len(set(truth)) < 2 else 'Descriptive subgroup'))
    shuffled_y = np.random.default_rng(2026).permutation(y)
    negative_scores = np.empty(len(y))
    for train, heldout in StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42).split(x, shuffled_y, genes):
        negative_scores[heldout] = fit_baseline(x[train], shuffled_y[train]).predict_proba(x[heldout])[:, 1]
    report = dict(training_rows=len(rows), gene_symbols=len(set(symbols)), biological_groups=len(set(genes)),
                  grouping='Connected components sharing a symbol, NCBI GeneID or exact protein sequence',
                  source_sha256=checksum(source), sklearn_version=sklearn.__version__,
                  seeds=list(seeds), folds=fold_reports, identical_sequence_groups=duplicate_sequences,
                  paired_uncertainty=paired_gene_bootstrap(y, first_scores['baseline'], first_scores['context'], genes, bootstrap),
                  calibration=dict(brier=float(brier_score_loss(y, first_scores['baseline'])), bins=score_bins(y, first_scores['baseline'])),
                  shuffled_label_negative_control=metrics(shuffled_y, negative_scores >= .5, negative_scores),
                  performance_by_review_status=by_review, gene_rows_distribution=dict(Counter(Counter(genes).values())),
                  limitations=[
                      'Existing test split has already been evaluated; no test rows loaded by this audit.',
                      'Three repeated five-fold runs are dependent, not 15 independent experiments.',
                      'Gene bootstrap intervals condition on fixed seed-42 OOF predictions; they omit refitting uncertainty and prior adaptive feature selection.',
                      'Most genes have one variant; gene-level scores are not stable generalization estimates.',
                      'Gene separation does not establish protein-family, ancestry, temporal, or clinical-population independence.',
                      'Labels are ClinVar aggregates across conditions; review status does not establish label correctness.',
                      'Balanced gene-diverse sampling, review imbalance and ortholog missingness limit calibration and deployment transfer.',
                      'One shuffled-label run is a diagnostic, not proof that leakage is absent.'
                  ], decision='Retain eight-feature baseline; context remains experimental. Independent reserved cohort required before final selection.')
    output.mkdir(parents=True)
    (output / 'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    for filename, records in [('predictions.csv', prediction_records), ('by_gene.csv', by_gene)]:
        with (output / filename).open('x', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('training_csv', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    result = validate(args.training_csv, args.output)
    print(json.dumps({'paired_uncertainty': result['paired_uncertainty'], 'negative_control': result['shuffled_label_negative_control']}, indent=2))
