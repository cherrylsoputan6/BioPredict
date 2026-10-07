import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from dataset import REQUIRED, build_dataset, candidate_from_row


def example(**updates):
    row = dict.fromkeys(REQUIRED, '-')
    row.update({
        '#AlleleID': '1', 'VariationID': '2', 'GeneSymbol': 'EXAMPLE',
        'Name': 'NM_123456.1(EXAMPLE):c.5G>A (p.Arg2His)',
        'Assembly': 'GRCh38', 'ClinicalSignificance': 'Pathogenic',
        'ReviewStatus': 'criteria provided, single submitter',
    })
    row.update(updates)
    return row


class DatasetTests(unittest.TestCase):
    def test_substitution_and_provenance(self):
        candidate, reason = candidate_from_row(example())
        self.assertIsNone(reason)
        self.assertEqual(candidate['variant'], 'R2H')
        self.assertEqual(candidate['transcript'], 'NM_123456.1')
        self.assertEqual(candidate['label'], 1)
        self.assertEqual(candidate['review_status'], 'criteria provided, single submitter')

    def test_labels_and_exclusions(self):
        for value in ('Benign', 'Likely benign', 'Benign/Likely benign'):
            self.assertEqual(candidate_from_row(example(ClinicalSignificance=value))[0]['label'], 0)
        for value in ('Uncertain significance', 'Conflicting classifications of pathogenicity', 'Pathogenic; risk factor'):
            self.assertEqual(candidate_from_row(example(ClinicalSignificance=value))[1], 'excluded_classification')

    def test_reject_unsupported_changes(self):
        for change in ('Arg2Ter', 'Arg2Arg', 'Arg2del', 'Arg2Hisfs', 'Arg0His', 'Xaa2His'):
            self.assertIsNone(candidate_from_row(example(Name=f'NM_123456.1(EXAMPLE):c.5G>A (p.{change})'))[0])
        self.assertEqual(candidate_from_row(example(GeneSymbol='OTHER'))[1], 'ambiguous_gene')
        self.assertEqual(candidate_from_row(example(Assembly='GRCh37'))[1], 'other_assembly')

    def test_plain_and_gzip_deduplication(self):
        for compressed in (False, True):
            with tempfile.TemporaryDirectory() as directory:
                source = Path(directory) / ('input.txt.gz' if compressed else 'input.txt')
                output = Path(directory) / 'output.csv'
                opener = gzip.open if compressed else open
                with opener(source, 'wt', newline='', encoding='utf-8') as handle:
                    writer = csv.DictWriter(handle, fieldnames=sorted(REQUIRED), delimiter='\t')
                    writer.writeheader()
                    writer.writerows([example(), example(), example(Assembly='GRCh37'), example(ClinicalSignificance='Uncertain significance')])
                counts = build_dataset(source, output)
                self.assertEqual(counts['rows_read'], 4)
                self.assertEqual(counts['retained'], 1)
                self.assertEqual(counts['duplicate_variation_id'], 1)
                with output.open(newline='') as handle:
                    self.assertEqual(len(list(csv.DictReader(handle))), 1)

    def test_invalid_header_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / 'bad.txt', Path(directory) / 'output.csv'
            source.write_text('wrong\nheader\n')
            with self.assertRaises(ValueError):
                build_dataset(source, output)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
