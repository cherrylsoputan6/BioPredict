"""Select a reproducible, gene-diverse pilot sample, not a train/test split."""
import argparse
import csv
import hashlib
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

from dataset import FIELDS, AMINO_ACIDS


def select_sample(source, output, per_label=50, seed=42):
    source, output = Path(source), Path(output)
    report = output.with_suffix('.summary.json')
    if per_label < 1:
        raise ValueError('per-label must be positive')
    if output.exists() or report.exists():
        raise ValueError('Output or summary already exists; choose a new filename')
    counts, groups = Counter(), {}
    amino_acids = set(AMINO_ACIDS.values())
    with source.open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        if not set(FIELDS) <= set(reader.fieldnames or []):
            raise ValueError('Use the original candidate CSV with all provenance columns')
        for row in reader:
            counts['rows_read'] += 1
            try:
                position = int(row['position'])
                if (None in row or any(row.get(f) is None for f in FIELDS) or
                        row['label'] not in ('0', '1') or position < 1 or
                        row['ref'] not in amino_acids or row['alt'] not in amino_acids or
                        row['ref'] == row['alt'] or not row['gene'] or not row['transcript'] or
                        not row['variation_id'].isdigit() or not row['allele_id'].isdigit() or
                        row['variant'] != f'{row["ref"]}{position}{row["alt"]}'):
                    raise ValueError('Malformed candidate')
            except (ValueError, TypeError):
                counts['malformed_rows'] += 1
                continue
            if position == 1:
                counts['position_one_excluded'] += 1
                continue
            key = (row['gene'], row['transcript'], row['variant'])
            if key not in groups:
                groups[key] = {'row': {f: row[f] for f in FIELDS}, 'labels': set(), 'count': 0}
            group = groups[key]
            group['labels'].add(row['label'])
            group['count'] += 1
            # Deterministic representative; retain the chosen ClinVar record's provenance.
            if (int(row['variation_id']), int(row['allele_id'])) < (
                    int(group['row']['variation_id']), int(group['row']['allele_id'])):
                group['row'] = {f: row[f] for f in FIELDS}

    pools = {label: defaultdict(list) for label in ('0', '1')}
    for group in groups.values():
        if len(group['labels']) > 1:
            counts['conflicting_substitution_groups'] += 1
            counts['conflicting_rows_excluded'] += group['count']
            continue
        counts['duplicate_rows_removed'] += group['count'] - 1
        row = group['row']
        pools[row['label']][row['gene']].append(row)

    rng, selected = random.Random(seed), []
    for label in ('0', '1'):
        genes = sorted(pools[label])
        counts[f'eligible_label_{label}'] = sum(len(v) for v in pools[label].values())
        if counts[f'eligible_label_{label}'] < per_label:
            raise ValueError(f'Not enough label {label} candidates for {per_label} rows')
        rng.shuffle(genes)
        for gene in genes:
            pools[label][gene].sort(key=lambda r: (r['transcript'], r['variant']))
            rng.shuffle(pools[label][gene])
        picked = []
        # One variant per gene per round, preventing large genes from dominating.
        while len(picked) < per_label:
            for gene in genes:
                if pools[label][gene]:
                    picked.append(pools[label][gene].pop())
                    if len(picked) == per_label:
                        break
        counts[f'selected_label_{label}'] = len(picked)
        counts[f'selected_genes_label_{label}'] = len({r['gene'] for r in picked})
        selected.extend(picked)
    rng.shuffle(selected)
    counts['selected'] = len(selected)
    counts['selected_genes'] = len({r['gene'] for r in selected})
    with source.open('rb') as handle:
        source_hash = hashlib.file_digest(handle, 'sha256').hexdigest()
    summary = dict(counts, seed=seed, per_label=per_label, source_sha256=source_hash,
                   duplicate_key=['gene', 'transcript', 'variant'],
                   sampling='equal labels; one variant per gene per round',
                   scope='pilot only; validate before training')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(selected)
    with report.open('x', encoding='utf-8') as handle:
        json.dump(summary, handle, indent=2)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--per-label', type=int, default=50)
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()
    try:
        print(json.dumps(select_sample(args.source, args.output, args.per_label, args.seed), indent=2))
    except ValueError as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
