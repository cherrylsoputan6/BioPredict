import unittest
from training_diagnostics import score_bins


class DiagnosticTests(unittest.TestCase):
    def test_bin_boundaries_and_empty_bins(self):
        bins = score_bins([0, 1, 1], [0.0, 0.2, 1.0])
        self.assertEqual(sum(b['rows'] for b in bins), 3)
        self.assertEqual(bins[0]['rows'], 1)
        self.assertEqual(bins[1]['rows'], 1)
        self.assertEqual(bins[4]['rows'], 1)
        self.assertEqual(bins[1]['observed_label_1_fraction'], 1)
        self.assertIsNone(bins[2]['mean_score'])
