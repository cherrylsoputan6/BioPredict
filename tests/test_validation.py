import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from validate_dataset import checked_protein, validate_batch, SequenceCache


def record():
    result = SeqRecord(Seq('ATGCGTTAA'), id='NM_123.1')
    result.annotations.update(organism='Homo sapiens', molecule_type='DNA')
    result.features = [SeqFeature(FeatureLocation(0, 9), type='CDS',
                       qualifiers={'gene': ['EXAMPLE'], 'translation': ['MR']})]
    return result


class ValidationTests(unittest.TestCase):
    def test_explicit_gene_synonym_only(self):
        renamed = record()
        renamed.features[0].qualifiers.update(gene=['NEWNAME'], gene_synonym=['OTHER; EXAMPLE'])
        self.assertEqual(checked_protein(renamed, 'EXAMPLE', 'NM_123.1'), 'MR')
        with self.assertRaises(ValueError):
            checked_protein(renamed, 'EXAM', 'NM_123.1')
        renamed.features.append(renamed.features[0])
        with self.assertRaises(ValueError):
            checked_protein(renamed, 'EXAMPLE', 'NM_123.1')

    def test_translation_exception_excluded(self):
        unusual = record()
        unusual.features[0].qualifiers['transl_except'] = ['(pos:4..6,aa:Sec)']
        with self.assertRaisesRegex(ValueError, 'Unsupported CDS translation'):
            checked_protein(unusual, 'EXAMPLE', 'NM_123.1')

    def test_annotation_and_accession(self):
        self.assertEqual(checked_protein(record(), 'EXAMPLE', 'NM_123.1'), 'MR')
        with self.assertRaises(ValueError):
            checked_protein(record(), 'EXAMPLE', 'NM_123.2')
        bad = record()
        bad.features[0].qualifiers['translation'] = ['MH']
        with self.assertRaises(ValueError):
            checked_protein(bad, 'EXAMPLE', 'NM_123.1')

    def test_disk_cache_and_mane(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch('validate_dataset.find_mane_transcript', return_value='NM_123.1'), patch('validate_dataset.download_genbank_record', return_value=record()) as download:
                cache = SequenceCache(directory)
                self.assertIsNone(cache.get('EXAMPLE', 'NM_123.2'))
                self.assertEqual(cache.get('EXAMPLE', 'NM_123.1'), 'MR')
                self.assertEqual(cache.get('EXAMPLE', 'NM_123.1'), 'MR')
                self.assertEqual(SequenceCache(directory).get('EXAMPLE', 'NM_123.1'), 'MR')
                download.assert_called_once()

    def test_batch_features_failures_and_limit(self):
        fields = ['gene', 'transcript', 'variant', 'ref', 'alt', 'position', 'label']
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / 'in.csv', Path(directory) / 'out.csv'
            with source.open('w', newline='') as handle:
                writer = csv.writer(handle)
                writer.writerow(fields)
                writer.writerows([
                    ['EXAMPLE', 'NM_123.1', 'R2H', 'R', 'H', 2, 1],
                    ['EXAMPLE', 'NM_123.1', 'R1H', 'R', 'H', 1, 0],
                    ['OTHER', 'NM_124.1', 'R2H', 'R', 'H', 2, 1],
                    ['ERROR', 'NM_125.1', 'R2H', 'R', 'H', 2, 1],
                    ['IGNORED', 'NM_126.1', 'R2H', 'R', 'H', 2, 1],
                ])
            cache = Mock()
            cache.get.side_effect = ['MR', 'MR', None, OSError('offline')]
            counts = validate_batch(source, output, cache, limit=4)
            self.assertEqual(counts, dict(processed=4, valid=1, reference_mismatch=1,
                                          not_current_mane=1, sequence_error=1))
            with output.open(newline='') as handle:
                rows = list(csv.DictReader(handle))
            self.assertAlmostEqual(float(rows[0]['weight_change']), -19.04)
            self.assertAlmostEqual(float(rows[0]['relative_position']), 1)
            self.assertEqual(rows[1]['weight_change'], '')
            with self.assertRaises(ValueError):
                validate_batch(source, output, cache)
