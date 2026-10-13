import unittest
from extended_features import extended_values


class ExtendedFeaturesTests(unittest.TestCase):
    def test_absolute_deltas_and_matrix(self):
        row = dict(ref='A', alt='D', relative_position=0.5,
                   hydrophobicity_change=-5.3, weight_change=44.01)
        values = extended_values(row)
        self.assertAlmostEqual(values[3], 5.3)
        self.assertAlmostEqual(values[4], 44.01)
        self.assertEqual(values[5], -2)
        reverse = dict(row, ref='D', alt='A', hydrophobicity_change=5.3, weight_change=-44.01)
        reverse_values = extended_values(reverse)
        self.assertEqual(values[3:], reverse_values[3:])

    def test_invalid_identity_and_inconsistent_features(self):
        row = dict(ref='A', alt='D', relative_position=0.5,
                   hydrophobicity_change=-5.3, weight_change=44.01)
        for update in (dict(alt='*'), dict(alt='A'), dict(weight_change=0)):
            with self.assertRaises(ValueError):
                extended_values(dict(row, **update))
