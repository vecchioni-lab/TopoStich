"""Renaming must preserve both public APIs and executable module aliases."""
import importlib
import subprocess
import sys
import unittest


class RenameTests(unittest.TestCase):
    def test_new_and_legacy_modules_report_new_brand(self):
        for module in ('topotrace', 'topotrace.cgd', 'topostich', 'topostich.cgd'):
            with self.subTest(module=module):
                result = subprocess.run([sys.executable, '-m', module, '--version'],
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.strip(), 'TopoTrace 0.3.1')

    def test_public_api_identity_and_packaged_catalog(self):
        for name, symbol in (('converter', 'convert_file'), ('converter', 'ConversionError'),
                             ('cgd', 'export_cgd'), ('nucleotides', 'recognize_residue'),
                             ('graph_viewer', 'marker_xml'), ('mmcif_input', 'read_mmcif')):
            with self.subTest(module=name, symbol=symbol):
                old = importlib.import_module('topostich.' + name)
                new = importlib.import_module('topotrace.' + name)
                self.assertIs(getattr(old, symbol), getattr(new, symbol))
        from topotrace.nucleotides import classify_atom
        self.assertEqual(classify_atom('DZ', 'C1').role, 'base')
