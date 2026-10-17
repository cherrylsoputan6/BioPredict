import unittest
from cohort_registry import novel_candidate


class CohortGateTests(unittest.TestCase):
    def row(self):
        return dict(VariationID='99', GeneID='7', ReviewStatus='reviewed by expert panel',
                    Assembly='GRCh38', ClinicalSignificance='Benign',
                    Name='NM_123.1(EXAMPLE):c.5G>A (p.Arg2His)', GeneSymbol='EXAMPLE',
                    **{'#AlleleID':'1', 'RCVaccession':'RCV1', 'LastEvaluated':'2026-10-13'})
    def test_id_and_stable_gene_exclusions(self):
        row = self.row()
        self.assertEqual(novel_candidate(row, {'99'}, set())[1], 'previously_present_variant')
        self.assertEqual(novel_candidate(row, set(), {'7'})[1], 'development_gene')
        candidate, reason = novel_candidate(row, set(), set())
        self.assertIsNone(reason)
        self.assertEqual(candidate['gene_id'], '7')
        self.assertEqual(novel_candidate(dict(row, ReviewStatus='no assertion criteria provided'),set(),set())[1], 'insufficient_label_evidence')
