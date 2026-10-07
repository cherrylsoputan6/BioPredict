import csv
import tempfile
import unittest
from pathlib import Path
from train_model import fit_baseline, train
from prepare_dataset import MODEL_FEATURES


class TrainingTests(unittest.TestCase):
    def test_scaler_fits_training_only(self):
        x = [[0.2, -2, -10], [0.4, 2, 10], [0.6, -1, -5], [0.8, 1, 5]]
        model = fit_baseline(x, [0, 1, 0, 1])
        before = model.named_steps['scaler'].mean_.copy()
        model.predict([[0.9, 1000, 10000]])
        self.assertEqual(model.named_steps['scaler'].mean_.tolist(), before.tolist())
        self.assertAlmostEqual(before[0], 0.5)
        self.assertEqual(before[1], 0)

    def test_artifacts_and_gene_overlap_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fields = ['gene', 'transcript', 'variant', 'label', 'validation_status'] + MODEL_FEATURES
            for name, prefix in [('train', 'A'), ('test', 'B')]:
                with (root / f'{name}.csv').open('w', newline='') as handle:
                    writer = csv.writer(handle)
                    writer.writerow(fields)
                    for i in range(4):
                        writer.writerow([f'{prefix}{i}', 'NM_123.1', 'R2H', i % 2, 'valid',
                                         0.2, (i % 2) * 2 - 1, (i % 2) * 10 - 5])
            report = train(root, root / 'out')
            self.assertEqual(report['gene_overlap'], 0)
            self.assertEqual(report['features'], MODEL_FEATURES)
            self.assertTrue((root / 'out/model.joblib').exists())
            self.assertEqual(report['logistic_regression']['confusion_matrix'], [[2, 0], [0, 2]])
            with self.assertRaises(ValueError):
                train(root, root / 'out')
            (root / 'test.csv').write_bytes((root / 'train.csv').read_bytes())
            with self.assertRaisesRegex(ValueError, 'Gene overlap'):
                train(root, root / 'bad')
            self.assertFalse((root / 'bad').exists())
