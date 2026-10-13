"""Pilot residue agreement across five species using Ensembl ortholog alignments.

This unweighted agreement descriptor is not a calibrated conservation predictor.
Only one-to-one orthologs whose aligned human protein exactly matches MANE qualify.
"""
import argparse
import csv
import hashlib
import json
import re
import time
from collections import Counter
from datetime import datetime, timezone
from itertools import islice
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

from Bio import SeqIO
from validate_dataset import checked_protein

SPECIES = {'mus_musculus', 'canis_lupus_familiaris', 'bos_taurus', 'gallus_gallus', 'danio_rerio'}
CANONICAL = set('ARNDCQEGHILKMFPSTWYV')
SERVERS = {'current': 'https://rest.ensembl.org',
           'oct2024': 'https://oct2024.rest.ensembl.org'}
OUTPUT_FIELDS = ['ortholog_status', 'ortholog_detail', 'ortholog_observations',
                 'ortholog_ref_agreement', 'ortholog_alt_agreement', 'ortholog_species']


def residue_observations(payload, protein, position):
    if not 1 <= position <= len(protein):
        raise ValueError('Position outside MANE protein')
    observations = {}
    matched = 0
    for entry in payload.get('data', []):
        for homology in entry.get('homologies', []):
            source, target = homology.get('source', {}), homology.get('target', {})
            species = target.get('species')
            if homology.get('type') != 'ortholog_one2one' or species not in SPECIES:
                continue
            human = source.get('align_seq', '').upper()
            other = target.get('align_seq', '').upper()
            if source.get('species') != 'homo_sapiens' or human.replace('-', '') != protein:
                continue
            if len(human) != len(other):
                raise ValueError('Unequal aligned sequence lengths')
            matched += 1
            residue = 0
            for column, amino_acid in enumerate(human):
                if amino_acid != '-':
                    residue += 1
                if amino_acid != '-' and residue == position:
                    observed = other[column]
                    if observed in CANONICAL:
                        if species in observations and observations[species] != observed:
                            raise ValueError('Inconsistent ortholog observations for a species')
                        observations[species] = observed
                    break
    return observations, matched


def fetch_homology(gene, directory, timeout, server='current'):
    if not re.fullmatch(r'[A-Za-z0-9_-]+', gene):
        raise ValueError('Invalid gene symbol')
    directory = Path(directory)
    if server not in SERVERS:
        raise ValueError('Unsupported Ensembl server')
    if server != 'current':
        directory = directory / server
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f'{gene}.json'
    url = f'{SERVERS[server]}/homology/symbol/homo_sapiens/{quote(gene)}?type=orthologues;sequence=protein;aligned=1;content-type=application/json'
    if path.exists():
        saved = json.loads(path.read_text(encoding='utf-8'))
        if saved['url'] != url:
            raise ValueError('Cache URL differs from requested source')
        return saved['payload']
    # Keep legacy full responses, but fetch missing genes one species at a time.
    # Persist each successful species so interrupted runs resume without repeats.
    combined = {'data': []}
    for species in sorted(SPECIES):
        species_path = directory / f'{gene}.{species}.json'
        species_url = url + f';target_species={species}'
        payload = cached_request(species_url, species_path, timeout, f'{gene}: {species}')
        combined['data'].extend(payload['data'])
    return combined


def cached_request(url, path, timeout, description):
    if path.exists():
        saved = json.loads(path.read_text(encoding='utf-8'))
        if saved['url'] != url:
            raise ValueError('Cache URL differs from requested source')
        return saved['payload']
    print(f'Fetching Ensembl orthologs for {description}...', flush=True)
    time.sleep(0.15)
    with urlopen(Request(url, headers={'Accept': 'application/json'}), timeout=timeout) as response:
        raw = response.read()
    payload = json.loads(raw)
    if not isinstance(payload, dict) or not isinstance(payload.get('data'), list):
        raise ValueError('Unexpected Ensembl response')
    saved = dict(url=url, fetched_utc=datetime.now(timezone.utc).isoformat(),
                 response_sha256=hashlib.sha256(raw).hexdigest(), payload=payload)
    path.write_text(json.dumps(saved), encoding='utf-8')
    return payload


def pilot(source, output, sequence_cache, ortholog_cache, limit=10, timeout=30, server='current'):
    output = Path(output)
    if output.exists():
        raise ValueError('Output exists; choose a new filename')
    if limit < 1 or not 0 < timeout <= 300:
        raise ValueError('Use a positive limit and timeout up to 300 seconds')
    counts = Counter()
    with Path(source).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle)
        required = {'gene', 'transcript', 'position', 'ref', 'alt', 'validation_status'}
        if not required <= set(reader.fieldnames or []) or set(OUTPUT_FIELDS) & set(reader.fieldnames):
            raise ValueError('Use a validated training CSV without ortholog columns')
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x', encoding='utf-8', newline='') as destination:
            writer = csv.DictWriter(destination, fieldnames=reader.fieldnames + OUTPUT_FIELDS)
            writer.writeheader()
            for row in islice(reader, limit):
                print(f'Checking orthologs: {row["gene"]} {row["ref"]}{row["position"]}{row["alt"]}', flush=True)
                status, detail = 'ortholog_error', ''
                try:
                    if row['validation_status'] != 'valid' or not re.fullmatch(r'NM_\d+\.\d+', row['transcript']):
                        raise ValueError('Expected validated MANE candidate')
                    record = SeqIO.read(Path(sequence_cache) / f'{row["transcript"]}.gb', 'genbank')
                    protein = checked_protein(record, row['gene'], row['transcript'])
                    position = int(row['position'])
                    if not 1 <= position <= len(protein) or protein[position - 1] != row['ref'] or row['alt'] not in CANONICAL:
                        raise ValueError('Candidate reference or alternate is invalid')
                    observations, matched = residue_observations(fetch_homology(row['gene'], ortholog_cache, timeout, server), protein, position)
                    row['ortholog_observations'] = len(observations)
                    row['ortholog_species'] = ';'.join(sorted(observations))
                    if len(observations) >= 3:
                        status = 'available'
                        row['ortholog_ref_agreement'] = sum(v == row['ref'] for v in observations.values()) / len(observations)
                        row['ortholog_alt_agreement'] = sum(v == row['alt'] for v in observations.values()) / len(observations)
                    else:
                        status = 'insufficient_orthologs' if matched else 'no_exact_mane_alignment'
                except (ValueError, OSError, RuntimeError) as error:
                    detail = str(error)
                row.update(ortholog_status=status, ortholog_detail=detail)
                writer.writerow(row)
                destination.flush()
                counts[status] += 1
                counts['processed'] += 1
                print(f'  {status}' + (f': {detail}' if detail else ''), flush=True)
    return dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('training_csv', type=Path)
    parser.add_argument('output_csv', type=Path)
    parser.add_argument('--limit', type=int, default=10)
    parser.add_argument('--timeout', type=float, default=30)
    parser.add_argument('--server', choices=SERVERS, default='current')
    parser.add_argument('--sequence-cache', type=Path, default=Path('data/sequence_cache'))
    parser.add_argument('--ortholog-cache', type=Path, default=Path('data/ortholog_cache'))
    args = parser.parse_args()
    print(json.dumps(pilot(args.training_csv, args.output_csv, args.sequence_cache,
                          args.ortholog_cache, args.limit, args.timeout, args.server), indent=2))


if __name__ == '__main__':
    main()
