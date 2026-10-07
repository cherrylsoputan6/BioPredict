import csv
import tempfile
import unittest
from pathlib import Path
from dataset import FIELDS
from sample_dataset import select_sample


def row(gene, variant, label, identifier):
    value = dict.fromkeys(FIELDS, '-')
    value.update(gene=gene, transcript='NM_123.1', variant=variant,
                 ref=variant[0], alt=variant[-1], position=variant[1:-1],
                 label=str(label), variation_id=str(identifier), allele_id=str(identifier))
    return value


class SamplingTests(unittest.TestCase):
    def test_balance_diversity_conflicts_and_reproducibility(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'input.csv'
            rows = [row('A', 'R2H', 0, 1), row('A', 'R2H', 0, 2),
                    row('A', 'R3H', 0, 3), row('B', 'R2H', 0, 4),
                    row('C', 'R2H', 1, 5), row('D', 'R2H', 1, 6),
                    row('E', 'R2H', 0, 7), row('E', 'R2H', 1, 8),
                    row('F', 'M1H', 1, 9), row('E', 'R2H', 0, 10)]
            with source.open('w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=FIELDS)
                writer.writeheader()
                writer.writerows(rows)
            first, second = Path(directory) / 'one.csv', Path(directory) / 'two.csv'
            summary = select_sample(source, first, 2, 42)
            select_sample(source, second, 2, 42)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(summary['selected_label_0'], 2)
            self.assertEqual(summary['selected_label_1'], 2)
            self.assertEqual(summary['selected_genes'], 4)
            self.assertEqual(summary['duplicate_rows_removed'], 1)
            self.assertEqual(summary['conflicting_rows_excluded'], 3)
            self.assertEqual(summary['position_one_excluded'], 1)
            with self.assertRaises(ValueError):
                select_sample(source, first, 2)
            with self.assertRaises(ValueError):
                select_sample(source, Path(directory) / 'too_many.csv', 100)
            self.assertFalse((Path(directory) / 'too_many.csv').exists())
