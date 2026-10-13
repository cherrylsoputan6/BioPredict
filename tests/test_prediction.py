import unittest
from unittest.mock import Mock, patch
import sklearn
from predict_variant import parse_missense, feature_vector, predict, FEATURES


class PredictionTests(unittest.TestCase):
    def test_prediction_refuses_before_model_scoring(self):
        pipeline = Mock()
        artifact = dict(features=FEATURES, sklearn_version=sklearn.__version__,
                        configuration=dict(features=FEATURES, threshold_label_1=0.5,
                                           conservation=dict(server='oct2024')), pipeline=pipeline)
        cache = Mock()
        cache.mane = {}
        cache.get.return_value = 'MA'
        with patch('predict_variant.joblib.load', return_value=artifact), patch('predict_variant.find_mane_transcript', return_value='NM_123.1'), patch('predict_variant.SequenceCache', return_value=cache), patch('predict_variant.fetch_homology') as fetch:
            with self.assertRaisesRegex(ValueError, 'Reference residue'):
                predict('EXAMPLE', 'R2H', 'unused.joblib', 'cache', 'orthologs')
            fetch.assert_not_called()
            pipeline.predict_proba.assert_not_called()

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
