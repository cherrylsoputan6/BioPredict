import unittest
from sequence_context import context_values


class ContextTests(unittest.TestCase):
    def test_excludes_variant_site_and_handles_edges(self):
        values = context_values('GAG', 2)
        self.assertEqual(values, [-0.4, 1.0, 0.0])
        self.assertEqual(values, context_values('GRG', 2))
        self.assertEqual(context_values('AGG', 1), values)
        with self.assertRaises(ValueError):
            context_values('AGG', 4)
        with self.assertRaises(ValueError):
            context_values('AUA', 1)
