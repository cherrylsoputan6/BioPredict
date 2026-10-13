import csv
import tempfile
import unittest
from pathlib import Path
from cross_validate_model import cross_validate
from prepare_dataset import MODEL_FEATURES


class CrossValidationTests(unittest.TestCase):
    def test_grouped_folds_cover_all_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / 'train.csv', Path(directory) / 'report.json'
            with source.open('w', newline='') as handle:
                writer = csv.writer(handle)
                writer.writerow(['gene', 'transcript', 'variant', 'label', 'validation_status'] + MODEL_FEATURES)
                for gene in range(10):
                    for label in (0, 1):
                        writer.writerow([f'G{gene}', 'NM_123.1', f'R{label+2}H', label, 'valid',
                                         0.2, label * 2 - 1, label * 10 - 5])
            report = cross_validate(source, output, folds=5)
            self.assertEqual(report['rows'], 20)
            self.assertEqual(sum(f['validation_rows'] for f in report['fold_results']), 20)
            self.assertTrue(all(f['gene_overlap'] == 0 for f in report['fold_results']))
            self.assertEqual(report['out_of_fold']['confusion_matrix'], [[10, 0], [0, 10]])
            with self.assertRaises(ValueError):
                cross_validate(source, output)
            with self.assertRaises(ValueError):
                cross_validate(source, Path(directory) / 'bad.json', folds=11)
            self.assertFalse((Path(directory) / 'bad.json').exists())
            subset = cross_validate(source, Path(directory) / 'subset.json', feature_subset=['weight_change'])
            self.assertEqual(subset['features'], ['weight_change'])
            with self.assertRaises(ValueError):
                cross_validate(source, Path(directory) / 'leak.json', feature_subset=['label'])
