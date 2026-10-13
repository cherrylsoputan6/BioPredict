"""Extract retryable connection failures and merge a completed retry by candidate ID."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from dataset import FIELDS


def read_rows(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        if not set(FIELDS + ['validation_status', 'validation_detail']) <= set(fields):
            raise ValueError('Missing validated CSV columns')
        rows = list(reader)
    if any(None in r or any(r.get(f) is None for f in fields) for r in rows):
        raise ValueError('Malformed CSV row')
    return fields, rows


def retryable(row):
    detail = row['validation_detail'].lower()
    return row['validation_status'] == 'sequence_error' and any(
        token in detail for token in (
            'getaddrinfo failed', 'timed out', 'connection reset',
            'connection refused', 'remote end closed connection',
            'http error 429', 'http error 500', 'http error 502',
            'http error 503', 'http error 504',
        ))


def write_rows(output, fields, rows):
    output = Path(output)
    if output.exists():
        raise ValueError('Output exists; choose a new filename')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def extract(source, output):
    _, rows = read_rows(source)
    selected = [{f: r[f] for f in FIELDS} for r in rows if retryable(r)]
    write_rows(output, FIELDS, selected)
    return dict(input_rows=len(rows), retry_candidates=len(selected),
                retained_without_retry=len(rows) - len(selected))


def key(row):
    return tuple(row[f] for f in ('variation_id', 'gene', 'transcript', 'variant'))


def merge(original, retried, output):
    fields, rows = read_rows(original)
    retry_fields, replacements = read_rows(retried)
    if fields != retry_fields:
        raise ValueError('Validation column layouts differ')
    index = {key(r): r for r in rows}
    if len(index) != len(rows):
        raise ValueError('Duplicate original candidate identities')
    changes = {}
    for row in replacements:
        candidate_key = key(row)
        old = index.get(candidate_key)
        if candidate_key in changes or old is None or not retryable(old):
            raise ValueError('Retry contains a duplicate, unknown, or nonretryable candidate')
        if any(row[f] != old[f] for f in FIELDS):
            raise ValueError('Candidate provenance changed during retry')
        changes[candidate_key] = row
    expected = {key(r) for r in rows if retryable(r)}
    if set(changes) != expected:
        raise ValueError('Retry is incomplete; finish all retry candidates before merging')
    result = [changes.get(key(r), r) for r in rows]
    write_rows(output, fields, result)
    summary = dict(rows=len(result), replaced=len(changes),
                   statuses=dict(Counter(r['validation_status'] for r in result)),
                   valid_labels=dict(Counter(r['label'] for r in result if r['validation_status'] == 'valid')))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    extraction = sub.add_parser('extract')
    extraction.add_argument('source')
    extraction.add_argument('output')
    merging = sub.add_parser('merge')
    merging.add_argument('original')
    merging.add_argument('retried')
    merging.add_argument('output')
    args = parser.parse_args()
    try:
        result = extract(args.source, args.output) if args.action == 'extract' else merge(args.original, args.retried, args.output)
        print(json.dumps(result, indent=2))
    except ValueError as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
