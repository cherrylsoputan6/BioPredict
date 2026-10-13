import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch
from conservation import residue_observations, fetch_homology, SPECIES


def payload(human='M-AR', target='MVDR', species='mus_musculus'):
    return {'data': [{'homologies': [dict(type='ortholog_one2one',
        source=dict(species='homo_sapiens', align_seq=human),
        target=dict(species=species, align_seq=target))]}]}


class ConservationTests(unittest.TestCase):
    def test_species_requests_and_partial_cache_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            def fake_request(url, path, timeout, description):
                self.assertIn(';target_species=', url)
                return payload(species=url.split('target_species=')[-1])
            with patch('conservation.cached_request', side_effect=fake_request) as fetch:
                result = fetch_homology('EXAMPLE', directory, 30)
                self.assertEqual(fetch.call_count, len(SPECIES))
                self.assertEqual(len(result['data']), len(SPECIES))
            with patch('conservation.cached_request', side_effect=fake_request) as fetch:
                fetch_homology('EXAMPLE', directory, 30, server='oct2024')
                self.assertTrue(all(call.args[0].startswith('https://oct2024.rest.ensembl.org/') for call in fetch.call_args_list))
                self.assertTrue(all(call.args[1].parent == Path(directory) / 'oct2024' for call in fetch.call_args_list))

    def test_gap_aware_mapping(self):
        observed, matched = residue_observations(payload(), 'MAR', 2)
        self.assertEqual(observed, {'mus_musculus': 'D'})
        self.assertEqual(matched, 1)

    def test_isoform_mismatch_and_gap_are_not_scores(self):
        observed, matched = residue_observations(payload(human='M-QR'), 'MAR', 2)
        self.assertEqual((observed, matched), ({}, 0))
        observed, matched = residue_observations(payload(target='MV-R'), 'MAR', 2)
        self.assertEqual((observed, matched), ({}, 1))

    def test_malformed_alignment_and_conflicting_species(self):
        with self.assertRaises(ValueError):
            residue_observations(payload(target='MAR'), 'MAR', 2)
        data = payload()
        data['data'][0]['homologies'] += payload(target='MVAR')['data'][0]['homologies']
        with self.assertRaises(ValueError):
            residue_observations(data, 'MAR', 2)
