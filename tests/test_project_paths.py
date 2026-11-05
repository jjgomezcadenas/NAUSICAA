"""Historical cards open canonical files without directory compatibility links."""
import json
from pathlib import Path
import tempfile
import unittest
from ion_signal.project_paths import resolve_recorded_path, resolve_config_paths
from ion_signal.spice_storage import source_manifest_matches


class RecordedPaths(unittest.TestCase):
    def test_longest_prefix_and_chained_report_move(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'docs').mkdir()
            (root / 'docs/ltspice_paths.json').write_text(json.dumps({'moves': [
                {'old': 'LTSpice/old', 'new': 'LTSpice/reference'},
                {'old': 'LTSpice/old/report.pdf', 'new': 'latex_generated/example/report.pdf'}]}))
            (root / 'docs/generated_reports.json').write_text(json.dumps({'reports': {
                'example': {'group': 'ion_readout'}}}))
            expected = root / 'latex_generated/ion_readout/example/report.pdf'
            self.assertEqual(resolve_recorded_path('LTSpice/old/report.pdf', root=root), expected)
            self.assertEqual(resolve_recorded_path(root / 'LTSpice/old/report.pdf', root=root), expected)
            self.assertEqual(resolve_recorded_path('/external/source', root=root), Path('/external/source'))
            card = {'input': 'LTSpice/old/circuit.cir', 'nested': ['latex_generated/example'], 'label': 'example'}
            resolved = resolve_config_paths(card, root=root)
            self.assertEqual(resolved['input'], 'LTSpice/reference/circuit.cir')
            self.assertEqual(resolved['nested'], ['latex_generated/ion_readout/example'])
            self.assertEqual(card['input'], 'LTSpice/old/circuit.cir')
            self.assertTrue(source_manifest_matches({'LTSpice/old/circuit.cir': 'hash'},
                {'LTSpice/reference/circuit.cir': 'hash'}, root))
            self.assertFalse(source_manifest_matches({'LTSpice/old/circuit.cir': 'changed'},
                {'LTSpice/reference/circuit.cir': 'hash'}, root))


if __name__ == '__main__':
    unittest.main()
