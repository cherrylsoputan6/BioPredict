"""Data-only export of the frozen scaler/logistic pipeline; no retraining."""
import json
from pathlib import Path
import numpy as np


class PortablePipeline:
    classes_ = np.array([0, 1])

    def __init__(self, record):
        self.mean = np.asarray(record['scaler_mean'], dtype=float)
        self.scale = np.asarray(record['scaler_scale'], dtype=float)
        self.coef = np.asarray(record['coefficients'], dtype=float)
        self.intercept = float(record['intercept'])
        self.n_features_in_ = len(self.mean)
        if (self.mean.shape != self.scale.shape or self.mean.shape != self.coef.shape or
                not np.all(self.scale > 0) or not all(np.all(np.isfinite(v)) for v in (self.mean, self.scale, self.coef)) or
                not np.isfinite(self.intercept)):
            raise ValueError('Invalid frozen model parameters')

    def predict_proba(self, values):
        values = np.asarray(values, dtype=float)
        if values.ndim != 2 or values.shape[1] != self.n_features_in_ or not np.all(np.isfinite(values)):
            raise ValueError('Invalid feature vector')
        logits = ((values - self.mean) / self.scale) @ self.coef + self.intercept
        score = np.exp(-np.logaddexp(0, -logits))
        return np.column_stack([1 - score, score])


def load_portable(path):
    record = json.loads(Path(path).read_text(encoding='utf-8'))
    return dict(record, pipeline=PortablePipeline(record))
