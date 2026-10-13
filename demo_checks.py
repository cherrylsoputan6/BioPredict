"""Replay benign development examples and refusal checks from local caches only."""
import argparse
import csv
import json
from pathlib import Path
import joblib
import sklearn
from Bio import SeqIO
from conservation import SPECIES
from predict_variant import FEATURES, parse_missense, feature_vector
from validate_dataset import checked_protein


def cached_payload(gene, cache):
    directory = Path(cache) / 'oct2024'
    full = directory / f'{gene}.json'
    if full.exists():
        return json.loads(full.read_text(encoding='utf-8'))['payload']
    combined = {'data': []}
    for species in sorted(SPECIES):
        record = json.loads((directory / f'{gene}.{species}.json').read_text(encoding='utf-8'))
        combined['data'].extend(record['payload']['data'])
    return combined


def run(cohort, model_path, output):
    output = Path(output)
    if output.exists():
        raise ValueError('Output exists; choose a new filename')
    artifact = joblib.load(model_path)
    if artifact['features'] != FEATURES or artifact['sklearn_version'] != sklearn.__version__:
        raise ValueError('Model environment or feature list differs')
    with Path(cohort).open(encoding='utf-8', newline='') as handle:
        examples = [r for r in csv.DictReader(handle) if r['label'] == '0'][:3]
    results = []
    for row in examples:
        record = SeqIO.read(Path('data/sequence_cache') / f'{row["transcript"]}.gb', 'genbank')
        protein = checked_protein(record, row['gene'], row['transcript'])
        payload = cached_payload(row['gene'], 'data/ortholog_cache')
        variant = parse_missense(row['variant'])
        values, observations = feature_vector(protein, variant, payload)
        score = float(artifact['pipeline'].predict_proba([values])[0, 1])
        results.append(dict(gene=row['gene'], variant=row['variant'], source_label=0,
                            model_label=int(score >= 0.5), score_label_1=score,
                            species=len(observations)))
        # These are behavior checks, not demands that the model agree with every label.
        for name, bad in [('wrong_reference', dict(variant, original='A' if variant['original'] != 'A' else 'R')),
                          ('out_of_range', dict(variant, position=len(protein) + 1))]:
            try:
                feature_vector(protein, bad, payload)
            except ValueError:
                results.append(dict(gene=row['gene'], check=name, result='refused'))
            else:
                raise AssertionError(f'{name} should have been refused')
    if len(examples) != 3:
        raise ValueError('Need three benign development examples')
    report = dict(results=results, scope='Cached development demonstrations only; not independent accuracy evidence. No model tuning or network access.')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort', type=Path, default=Path('data/conservation_comparison/matched_training.csv'))
    parser.add_argument('--model', type=Path, default=Path('data/conservation_model/model.joblib'))
    parser.add_argument('--output', type=Path, default=Path('data/demo_checks.json'))
    args = parser.parse_args()
    print(json.dumps(run(args.cohort, args.model, args.output), indent=2))


if __name__ == '__main__':
    main()
