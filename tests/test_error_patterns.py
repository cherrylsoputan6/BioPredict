import unittest
from error_patterns import summarize


class ErrorPatternTests(unittest.TestCase):
    def test_rates_use_class_denominators(self):
        base = dict(ref='A', alt='D', relative_position='0.5', hydrophobicity_change='-5.3',
                    weight_change='44.01', ortholog_ref_agreement='1', ortholog_alt_agreement='0')
        rows = [dict(base, label='0', model_label='1'), dict(base, label='0', model_label='0'),
                dict(base, label='1', model_label='0')]
        for group in summarize(rows):
            self.assertEqual(group['benign'], 2)
            self.assertEqual(group['pathogenic'], 1)
            self.assertEqual(group['false_positive_rate'], 0.5)
            self.assertEqual(group['false_negative_rate'], 1)
