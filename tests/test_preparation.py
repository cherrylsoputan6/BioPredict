import csv
import tempfile
import unittest
from pathlib import Path
from prepare_dataset import prepare


class PreparationTests(unittest.TestCase):
    def write_input(self, path, rows):
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def rows(self):
        return [dict(gene=f'G{i}', transcript='NM_123.1', variant=f'R{p}H',
                     label=str(i % 2), validation_status='valid', position=p,
                     protein_length=10, relative_position=p / 10,
                     hydrophobicity_change=1.3, weight_change=-19.04)
                for i in range(10) for p in (2, 3)]

    def test_group_split_and_reproducibility(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'input.csv'
            rows = self.rows()
            rows.append(dict(rows[0], validation_status='sequence_error'))
            self.write_input(source, rows)
            summary = prepare(source, root / 'one')
            prepare(source, root / 'two')
            self.assertEqual(summary['gene_overlap'], 0)
            self.assertEqual(summary['retained'], 20)
            self.assertEqual(summary['excluded_nonvalid'], 1)
            for name in ('train', 'test'):
                self.assertGreater(summary[name]['label_0'], 0)
                self.assertGreater(summary[name]['label_1'], 0)
                self.assertEqual((root / 'one' / f'{name}.csv').read_bytes(),
                                 (root / 'two' / f'{name}.csv').read_bytes())
            self.assertEqual(summary['test']['rows'], 4)
            with self.assertRaises(ValueError):
                prepare(source, root / 'one')

    def test_bad_features_and_duplicates_fail_before_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for case in ('nan', 'duplicate', 'one_class'):
                rows = self.rows()
                if case == 'nan':
                    rows[0]['weight_change'] = 'nan'
                elif case == 'duplicate':
                    rows.append(rows[0].copy())
                else:
                    for row in rows:
                        row['label'] = '0'
                source = root / f'{case}.csv'
                self.write_input(source, rows)
                with self.assertRaises(ValueError):
                    prepare(source, root / case)
                self.assertFalse((root / case).exists())
