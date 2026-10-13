import unittest
from compare_conservation import matched_cohort


class ConservationComparisonTests(unittest.TestCase):
    def test_reorders_and_filters_without_changing_training_values(self):
        training = [dict(variation_id=str(i), gene=f'G{i}', transcript='NM_123.1', variant='R2H', label=str(i % 2)) for i in range(3)]
        observations = [dict(r, ortholog_status='available' if i < 2 else 'ortholog_error') for i, r in enumerate(training)]
        cohort = matched_cohort(training, list(reversed(observations)))
        self.assertEqual([r['variation_id'] for r in cohort], ['0', '1'])
        with self.assertRaises(ValueError):
            matched_cohort(training, observations[:-1])
        observations[0]['label'] = '1'
        with self.assertRaises(ValueError):
            matched_cohort(training, observations)
