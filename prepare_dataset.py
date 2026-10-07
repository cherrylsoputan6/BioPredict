"""Prepare pilot train/test CSVs with disjoint genes and raw validated features."""
import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter
from pathlib import Path

MODEL_FEATURES = ['relative_position', 'hydrophobicity_change', 'weight_change']


def prepare(source, destination, seed=42, test_fraction=0.2):
    source, destination = Path(source), Path(destination)
    if destination.exists():
        raise ValueError('Destination already exists; choose a new folder')
    if not 0 < test_fraction < 1:
        raise ValueError('Test fraction must be between zero and one')
    counts, rows, seen = Counter(), [], set()
    with source.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        required = set(MODEL_FEATURES + ['gene', 'transcript', 'variant', 'label',
                                       'validation_status', 'position', 'protein_length'])
        if not required <= set(fields):
            raise ValueError('Missing validated dataset columns')
        for row in reader:
            counts['input_rows'] += 1
            if row['validation_status'] != 'valid':
                counts['excluded_nonvalid'] += 1
                continue
            try:
                values = [float(row[f]) for f in MODEL_FEATURES]
                position, length = int(row['position']), int(row['protein_length'])
                if (None in row or any(row.get(f) is None for f in fields) or
                        not all(math.isfinite(v) for v in values) or
                        row['label'] not in ('0', '1') or
                        not all(row[f] for f in ('gene', 'transcript', 'variant')) or
                        not 1 < position <= length or
                        not math.isclose(values[0], position / length, rel_tol=1e-9)):
                    raise ValueError('Invalid validated features')
            except (ValueError, TypeError, ZeroDivisionError):
                raise ValueError(f'Malformed valid row {counts["input_rows"]}; repair input first') from None
            key = (row['gene'], row['transcript'], row['variant'])
            if key in seen:
                raise ValueError('Duplicate substitution; deduplicate and resolve labels before splitting')
            seen.add(key)
            rows.append(row)
    genes = sorted({row['gene'] for row in rows})
    if len(genes) < 4 or {row['label'] for row in rows} != {'0', '1'}:
        raise ValueError('Need both labels and at least four genes')
    test_count = max(1, min(len(genes) - 1, round(len(genes) * test_fraction)))
    rng = random.Random(seed)
    # Find a group split with both classes; never inspect model performance.
    for attempt in range(1000):
        shuffled = genes.copy()
        rng.shuffle(shuffled)
        test_genes = set(shuffled[:test_count])
        train = [row for row in rows if row['gene'] not in test_genes]
        test = [row for row in rows if row['gene'] in test_genes]
        if all({r['label'] for r in part} == {'0', '1'} for part in (train, test)):
            break
    else:
        raise ValueError('Cannot create gene-disjoint sets with both labels; use a larger sample')
    overlap = {r['gene'] for r in train} & {r['gene'] for r in test}
    assert not overlap
    with source.open('rb') as handle:
        checksum = hashlib.file_digest(handle, 'sha256').hexdigest()
    summary = dict(counts, retained=len(rows), seed=seed, test_fraction=test_fraction,
                   split_attempts=attempt + 1, gene_overlap=len(overlap),
                   source_sha256=checksum, model_features=MODEL_FEATURES,
                   scope='pilot workflow check; no scaling or model fitting performed')
    for name, part in [('train', train), ('test', test)]:
        summary[name] = dict(rows=len(part), genes=len({r['gene'] for r in part}),
                            label_0=sum(r['label'] == '0' for r in part),
                            label_1=sum(r['label'] == '1' for r in part))
    destination.mkdir(parents=True)
    for name, part in [('clean', rows), ('train', train), ('test', test)]:
        with (destination / f'{name}.csv').open('x', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(part)
    (destination / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--test-fraction', type=float, default=0.2)
    args = parser.parse_args()
    try:
        print(json.dumps(prepare(args.source, args.destination, args.seed, args.test_fraction), indent=2))
    except ValueError as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
