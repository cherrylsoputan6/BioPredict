"""Freeze all baseline ClinVar VariationIDs; screen a later clinical snapshot.

Passing this gate establishes record/gene novelty only, not protein-family or
source independence. Sequences and family identity must be checked before use.
"""
import argparse
import csv
import gzip
import hashlib
import json
import sqlite3
from pathlib import Path
from dataset import candidate_from_row, FIELDS


def file_hash(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def novel_candidate(row, previous_ids, development_gene_ids):
    if row['VariationID'] in previous_ids:
        return None, 'previously_present_variant'
    gene_id = row.get('GeneID', '')
    if not gene_id.isdigit() or gene_id == '0':
        return None, 'unresolved_gene_id'
    if gene_id in development_gene_ids:
        return None, 'development_gene'
    if row.get('ReviewStatus') not in ('criteria provided, multiple submitters, no conflicts', 'reviewed by expert panel', 'practice guideline'):
        return None, 'insufficient_label_evidence'
    candidate, reason = candidate_from_row(row)
    if candidate and candidate['position'] == 1:
        return None, 'position_one'
    return (dict(candidate, gene_id=gene_id), None) if candidate else (None, reason)


def freeze(source, development, directory):
    directory = Path(directory)
    if directory.exists():
        raise ValueError('Registry exists; preserve it')
    with Path(development).open(encoding='utf-8', newline='') as handle:
        development_rows = list(csv.DictReader(handle))
    used_variants = {r['variation_id'] for r in development_rows if r.get('variation_id', '').isdigit()}
    used_symbols = {r['gene'] for r in development_rows}
    directory.mkdir(parents=True)
    connection = sqlite3.connect(directory / 'registry.sqlite')
    connection.executescript('CREATE TABLE variants(id TEXT PRIMARY KEY); CREATE TABLE genes(id TEXT PRIMARY KEY);')
    found_ids, batch, gene_ids = set(), [], set()
    number = 0
    opener = gzip.open if str(source).endswith('.gz') else open
    with opener(source, 'rt', encoding='utf-8', newline='') as handle:
        for number, row in enumerate(csv.DictReader(handle, delimiter='\t'), 1):
            identifier = row['VariationID']
            if identifier.isdigit():
                batch.append((identifier,))
            if identifier in used_variants or row['GeneSymbol'] in used_symbols:
                if row['GeneID'].isdigit() and row['GeneID'] != '0':
                    gene_ids.add(row['GeneID'])
                    found_ids.add(identifier)
            if len(batch) >= 50000:
                connection.executemany('INSERT OR IGNORE INTO variants VALUES (?)', batch)
                connection.commit()
                batch.clear()
            if number % 500000 == 0:
                print(f'Indexed {number} source rows', flush=True)
    connection.executemany('INSERT OR IGNORE INTO variants VALUES (?)', batch)
    connection.executemany('INSERT INTO genes VALUES (?)', [(g,) for g in sorted(gene_ids)])
    connection.commit()
    manifest = dict(source_sha256=file_hash(source), development_sha256=file_hash(development),
                    source_rows=number, unique_variation_ids=connection.execute('SELECT count(*) FROM variants').fetchone()[0],
                    excluded_gene_ids=len(gene_ids), unmatched_development_variants=sorted(used_variants - found_ids),
                    status='baseline registry only; no independent cohort collected',
                    required_next_checks=['verify later archived release provenance', 'screen novel IDs and stable genes',
                                          'exclude alias/identical-protein and family overlap', 'lock eligible cohort before scoring'])
    connection.close()
    (directory / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    if manifest['unmatched_development_variants']:
        raise ValueError('Unresolved development identities; do not use this registry until resolved')
    return manifest


def extend_cached_history(directory, candidates):
    """Add historical cohorts reconstructed from preserved candidate order/caches."""
    from validation_groups import record_identities
    directory = Path(directory)
    with Path(candidates).open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle))
    manifest_path = directory / 'manifest.json'
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    connection = sqlite3.connect(directory / 'registry.sqlite')
    gene_ids = set()
    for row in rows:
        if not connection.execute('SELECT 1 FROM variants WHERE id=?', (row['variation_id'],)).fetchone():
            raise ValueError('Historical ID absent from the baseline source registry')
        gene_ids.update(key.split(':', 1)[1] for key in record_identities(row['gene'], row['transcript']) if key.startswith('GeneID:'))
    connection.executemany('INSERT OR IGNORE INTO genes VALUES (?)', [(g,) for g in sorted(gene_ids)])
    connection.commit()
    manifest['excluded_gene_ids'] = connection.execute('SELECT count(*) FROM genes').fetchone()[0]
    connection.close()
    manifest.setdefault('history_extensions', []).append(dict(candidate_sha256=file_hash(candidates), rows=len(rows),
                                                             gene_ids=len(gene_ids), method='Stable GeneIDs from validated cached records; all historical IDs verified in full source'))
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline_summary')
    parser.add_argument('development_candidates')
    parser.add_argument('output_directory')
    args = parser.parse_args()
    print(json.dumps(freeze(args.baseline_summary, args.development_candidates, args.output_directory), indent=2))
