import unittest
from predict_variant import parse_missense, feature_vector


class PredictionTests(unittest.TestCase):
    def test_reject_unsupported_variants(self):
        self.assertEqual(parse_missense(' r2h ')['position'], 2)
        for text in ('M1I', 'R2R', 'R0H', 'R2*', 'R2Hfs', ''):
            with self.assertRaises(ValueError):
                parse_missense(text)

    def test_inference_features_and_coverage(self):
        homologies = [dict(type='ortholog_one2one', source=dict(species='homo_sapiens', align_seq='M-R'),
                          target=dict(species=species, align_seq=sequence))
                      for species, sequence in [('mus_musculus', 'MAR'), ('bos_taurus', 'MAH'), ('gallus_gallus', 'M-R')]]
        payload = dict(data=[dict(homologies=homologies)])
        vector, observations = feature_vector('MR', parse_missense('R2H'), payload)
        self.assertEqual(len(vector), 8)
        self.assertAlmostEqual(vector[1], 1.3)
        self.assertAlmostEqual(vector[2], -19.04)
        self.assertAlmostEqual(vector[6], 2 / 3)
        self.assertAlmostEqual(vector[7], 1 / 3)
        with self.assertRaises(ValueError):
            feature_vector('MA', parse_missense('R2H'), payload)
        payload['data'][0]['homologies'].pop()
        with self.assertRaises(ValueError):
            feature_vector('MR', parse_missense('R2H'), payload)
