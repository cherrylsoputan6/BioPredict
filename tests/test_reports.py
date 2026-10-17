import json
import unittest
from reporting import report_bytes


class ReportTests(unittest.TestCase):
    def test_reports_preserve_score_provenance_and_limits(self):
        result = dict(gene='EXAMPLE',variant='R2H',transcript='NM_123.1',model_class='benign',
                      score_label_1=.2,threshold=.5,ortholog_species=['mouse'],
                      interpretation='Uncalibrated research score',features={'x':1},
                      model_artifact_sha256='abc',created_utc='2026-10-13T00:00:00Z')
        data, mime = report_bytes(result,'json')
        self.assertEqual(json.loads(data)['score_label_1'], .2)
        self.assertEqual(json.loads(data)['independent_validation'],'not_completed')
        text, mime = report_bytes(result,'txt')
        self.assertIn(b'Model artifact SHA256: abc',text)
        self.assertIn(b'not a clinical classification',text)
        with self.assertRaises(ValueError):
            report_bytes(result,'exe')
