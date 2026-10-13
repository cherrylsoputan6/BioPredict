import csv
import tempfile
import unittest
from pathlib import Path
from dataset import FIELDS
from retry_validation import extract, merge, retryable


class RetryTests(unittest.TestCase):
    def test_network_failures_only(self):
        for detail in ('getaddrinfo failed', 'The handshake operation timed out', 'HTTP Error 503'):
            self.assertTrue(retryable(dict(validation_status='sequence_error', validation_detail=detail)))
        self.assertFalse(retryable(dict(validation_status='sequence_error', validation_detail='Unsupported CDS translation')))

    def test_merge_preserves_success_and_requires_complete_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fields = FIELDS + ['validation_status', 'validation_detail', 'weight_change']
            rows = []
            for identifier, status, detail in [('1', 'valid', ''), ('2', 'sequence_error', 'getaddrinfo failed'), ('3', 'sequence_error', 'Unsupported CDS translation')]:
                row = dict.fromkeys(fields, '-')
                row.update(variation_id=identifier, gene='EXAMPLE', transcript='NM_123.1',
                           variant='R2H', label='1', validation_status=status, validation_detail=detail)
                rows.append(row)
            def write(path, data):
                with path.open('w', newline='') as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(data)
            source, retry = root / 'original.csv', root / 'retried.csv'
            write(source, rows)
            self.assertEqual(extract(source, root / 'candidates.csv')['retry_candidates'], 1)
            write(retry, [])
            with self.assertRaises(ValueError):
                merge(source, retry, root / 'incomplete.csv')
            recovered = dict(rows[1], validation_status='valid', validation_detail='', weight_change='-19.04')
            write(retry, [recovered])
            summary = merge(source, retry, root / 'merged.csv')
            self.assertEqual(summary['statuses'], {'valid': 2, 'sequence_error': 1})
            with (root / 'merged.csv').open(newline='') as handle:
                result = list(csv.DictReader(handle))
            self.assertEqual(result[0], rows[0])
            self.assertEqual(result[2], rows[2])
            write(retry, [dict(recovered, label='0')])
            with self.assertRaises(ValueError):
                merge(source, retry, root / 'changed.csv')
