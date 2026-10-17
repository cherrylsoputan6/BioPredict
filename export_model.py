"""Verify a data-only model export against frozen training predictions."""
import argparse
import hashlib
import json
from pathlib import Path
import joblib
import numpy as np
from train_model import read_split
from evaluate_conservation import add_conservation
from portable_model import PortablePipeline


def export(source, training, output):
    output = Path(output)
    if output.exists():
        raise ValueError('Export exists; preserve existing work')
    artifact = joblib.load(source)
    pipeline = artifact['pipeline']
    scaler, classifier = pipeline.named_steps['scaler'], pipeline.named_steps['classifier']
    record = dict(format_version=1, model_id='biopredict-eight-feature-baseline-v1',
                  features=artifact['features'], configuration=artifact['configuration'],
                  sklearn_version=artifact['sklearn_version'], scaler_mean=scaler.mean_.tolist(),
                  scaler_scale=scaler.scale_.tolist(), coefficients=classifier.coef_[0].tolist(),
                  intercept=float(classifier.intercept_[0]), source_artifact_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest())
    rows, six, _ = read_split(training, extended=True)
    values = add_conservation(rows, six)
    difference = float(np.max(np.abs(pipeline.predict_proba(values) - PortablePipeline(record).predict_proba(values))))
    if difference > 1e-12:
        raise ValueError('Export differs from frozen pipeline')
    record['export_verification'] = dict(training_rows=len(rows), maximum_absolute_score_difference=difference,
                                       note='Training-only numerical equivalence check; no test evaluation or retraining')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(record, indent=2), encoding='utf-8')
    return record['export_verification']


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('training')
    parser.add_argument('output')
    args = parser.parse_args()
    print(json.dumps(export(args.source, args.training, args.output), indent=2))
