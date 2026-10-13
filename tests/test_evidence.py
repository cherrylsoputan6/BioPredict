import unittest
from evidence_sensitivity import filter_rows


class EvidenceTests(unittest.TestCase):
    def test_exact_status_filter_preserves_rows(self):
        rows = [dict(review_status='no assertion criteria provided', label='1'),
                dict(review_status='criteria provided, single submitter', label='0'),
                dict(review_status='reviewed by expert panel', label='1')]
        self.assertEqual(filter_rows(rows), rows[1:])
        self.assertEqual(len(rows), 3)

    def test_missing_status_is_not_silently_kept(self):
        with self.assertRaises(ValueError):
            filter_rows([dict(review_status='', label='0')])
