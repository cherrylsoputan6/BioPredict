import unittest
import numpy as np
from portable_model import PortablePipeline


class PortableModelTests(unittest.TestCase):
    def test_stable_logits_and_input_validation(self):
        model = PortablePipeline(dict(scaler_mean=[0], scaler_scale=[1], coefficients=[1], intercept=0))
        scores = model.predict_proba([[0], [1000], [-1000]])
        np.testing.assert_allclose(scores, [[.5,.5],[0,1],[1,0]])
        with self.assertRaises(ValueError):
            model.predict_proba([[float('nan')]])
        with self.assertRaises(ValueError):
            PortablePipeline(dict(scaler_mean=[0], scaler_scale=[0], coefficients=[1], intercept=0))
