"""Portable research reports; a score is never described as disease probability."""
import json


def report_bytes(result, kind):
    report = dict(result, independent_validation='not_completed',
                  report_scope='Research demonstration; not a clinical classification or probability of disease')
    if kind == 'json':
        return json.dumps(report, indent=2).encode(), 'application/json'
    if kind != 'txt':
        raise ValueError('Unsupported report format')
    text = '\n'.join([
        'BioPredict research report', f'{result["gene"]} {result["variant"]}',
        f'MANE transcript: {result["transcript"]}', f'Model class: {result["model_class"]}',
        f'Uncalibrated label-1 score: {result["score_label_1"]}', f'Threshold: {result["threshold"]}',
        'Species: ' + ', '.join(result['ortholog_species']),
        'Model: ' + result.get('model_id', 'baseline'),
        'Model artifact SHA256: ' + result.get('model_artifact_sha256', 'not_recorded'),
        'Created UTC: ' + result.get('created_utc', 'not_recorded'),
        'Archive: ' + result.get('source_archive', 'Ensembl October 2024'),
        '', report['report_scope'], 'Independent validation: not completed.',
        result['interpretation'], '', 'Feature values:', json.dumps(result['features'], indent=2),
    ])
    return text.encode(), 'text/plain; charset=utf-8'
