"""Research-only prediction using the frozen conservation pipeline."""
import argparse
import json
import re
import socket
from pathlib import Path

import joblib
import sklearn
from Bio import Entrez
from analysis import parse_variant, validate_variant, extract_variant_features
from conservation import fetch_homology, residue_observations
from extended_features import extended_values, EXTENDED_FEATURES
from ncbi import find_mane_transcript
from validate_dataset import SequenceCache

FEATURES = EXTENDED_FEATURES + ['ortholog_ref_agreement', 'ortholog_alt_agreement']


def parse_missense(text):
    text = text.strip().upper()
    if not re.fullmatch(r'[ARNDCQEGHILKMFPSTWYV][1-9]\d*[ARNDCQEGHILKMFPSTWYV]', text):
        raise ValueError('Use a single-letter substitution such as R175H')
    info = parse_variant(text)
    if info['position'] == 1 or info['original'] == info['new']:
        raise ValueError('The model excludes start-codon changes and unchanged residues')
    return info


def feature_vector(protein, variant, payload):
    if not validate_variant(protein, variant):
        raise ValueError('Reference residue does not match the MANE protein')
    observations, matched = residue_observations(payload, str(protein), variant['position'])
    if len(observations) < 3:
        reason = 'Too few qualifying species' if matched else 'No exact MANE ortholog alignment'
        raise ValueError(f'{reason}; no prediction available')
    raw = extract_variant_features(protein, variant)
    raw.update(ref=variant['original'], alt=variant['new'])
    values = extended_values(raw) + [
        sum(v == variant['original'] for v in observations.values()) / len(observations),
        sum(v == variant['new'] for v in observations.values()) / len(observations),
    ]
    return values, observations


def predict(gene, text, model_path, sequence_cache, ortholog_cache, timeout=60):
    if not re.fullmatch(r'[A-Za-z0-9_-]+', gene):
        raise ValueError('Enter a gene symbol')
    variant = parse_missense(text)
    # Only load a model created locally by this project.
    artifact = joblib.load(model_path)
    config = artifact['configuration']
    if (artifact['features'] != FEATURES or config['features'] != FEATURES or
            artifact['sklearn_version'] != sklearn.__version__):
        raise ValueError('Model features or scikit-learn version differ; use the matching training environment')
    if config['threshold_label_1'] != 0.5 or config['conservation']['server'] != 'oct2024':
        raise ValueError('Unsupported model configuration')
    previous_timeout, previous_tries = socket.getdefaulttimeout(), Entrez.max_tries
    socket.setdefaulttimeout(timeout)
    Entrez.max_tries = 1
    try:
        accession = find_mane_transcript(gene)
        if not accession:
            raise ValueError('No current MANE transcript found')
        cache = SequenceCache(sequence_cache)
        cache.mane[gene] = accession
        protein = cache.get(gene, accession)
        if not validate_variant(protein, variant):
            raise ValueError('Reference residue does not match the MANE protein')
        payload = fetch_homology(gene, ortholog_cache, timeout, server='oct2024')
        values, observations = feature_vector(protein, variant, payload)
    finally:
        socket.setdefaulttimeout(previous_timeout)
        Entrez.max_tries = previous_tries
    model = artifact['pipeline']
    if list(model.classes_) != [0, 1] or model.n_features_in_ != len(FEATURES):
        raise ValueError('Unexpected model classes or feature count')
    score = float(model.predict_proba([values])[0, 1])
    return dict(gene=gene, transcript=accession, variant=text.strip().upper(),
                model_label=int(score >= 0.5), model_class='pathogenic' if score >= 0.5 else 'benign',
                score_label_1=score, threshold=0.5, ortholog_species=sorted(observations),
                features=dict(zip(FEATURES, values)),
                interpretation='Research model output, not a clinical classification. Score is uncalibrated and is not a probability of disease. This variant may overlap the development dataset.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gene')
    parser.add_argument('variant')
    parser.add_argument('--model', type=Path, default=Path('data/conservation_model/model.joblib'))
    parser.add_argument('--sequence-cache', type=Path, default=Path('data/sequence_cache'))
    parser.add_argument('--ortholog-cache', type=Path, default=Path('data/ortholog_cache'))
    parser.add_argument('--timeout', type=float, default=60)
    args = parser.parse_args()
    if not 0 < args.timeout <= 300:
        parser.error('Timeout must be between zero and 300 seconds')
    try:
        result = predict(args.gene.strip(), args.variant, args.model, args.sequence_cache, args.ortholog_cache, args.timeout)
        print(json.dumps(result, indent=2))
    except (ValueError, OSError, RuntimeError) as error:
        parser.exit(1, f'Prediction unavailable: {error}\n')


if __name__ == '__main__':
    main()
