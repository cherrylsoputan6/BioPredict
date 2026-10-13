import unittest
import numpy as np
from scientific_validation import paired_gene_bootstrap
from validation_groups import connected_groups


class ScientificValidationTests(unittest.TestCase):
    def test_aliases_and_identical_proteins_stay_together(self):
        groups = connected_groups(['OLD', 'NEW', 'OTHER', 'DISTINCT'],
                                  [['GeneID:1', 'seq:a'], ['GeneID:1', 'seq:b'],
                                   ['GeneID:2', 'seq:b'], ['GeneID:3', 'seq:c']])
        self.assertEqual(len(set(groups[:3])), 1)
        self.assertNotEqual(groups[0], groups[3])

    def test_identical_models_have_zero_paired_difference(self):
        y = np.array([0, 1, 0, 1, 0, 1])
        scores = np.array([.2, .8, .6, .4, .1, .9])
        report = paired_gene_bootstrap(y, scores, scores, np.array(['A','A','B','B','C','C']), repetitions=30)
        for measure in report.values():
            self.assertEqual(measure['context_minus_baseline'], 0)
            self.assertEqual(measure['conditional_95_percentile_interval'], [0, 0])
