import unittest
from evaluate_conservation import add_conservation


class EvaluationFeatureTests(unittest.TestCase):
    def test_features_do_not_mutate_input(self):
        values = [[0.5, 1.3]]
        row = dict(ortholog_status='available', ortholog_observations='3',
                   ortholog_ref_agreement='1', ortholog_alt_agreement='0')
        self.assertEqual(add_conservation([row], values), [[0.5, 1.3, 1, 0]])
        self.assertEqual(values, [[0.5, 1.3]])
        for change in (dict(ortholog_ref_agreement='nan'), dict(ortholog_observations='2'),
                       dict(ortholog_alt_agreement='1')):
            with self.assertRaises(ValueError):
                add_conservation([dict(row, **change)], values)
