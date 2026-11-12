"""Backbone nick stitching regressions on small synthetic crystals.

The screen is ported from NAtopo's find_backbone_join_candidates: loose P
(with own O5') to loose O3' (with own C3') across different chain copies.
"""
from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import gemmi

from test_converter import crystal
from topostich.cgd import export_cgd
from topostich.converter import ConversionError

SCHEMAS = Path(__file__).resolve().parents[1] / 'schemas'


def nucleotide(serial, residue, chain, positions, residue_name='DG'):
    """PDB records for one residue; positions maps atom name -> xyz."""
    lines = []
    for offset, (name, xyz) in enumerate(positions.items()):
        lines.append(f"ATOM  {serial+offset:5d} {name:>4} {residue_name:>3} {chain:1}{residue:4d}    "
                     f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.:6.2f}{20.:6.2f}          "
                     f" {name[0]}\n")
    return ''.join(lines)


def line_residue(x0, y=5., z=5., *, phosphate=True, o3=True, step=1.):
    """Atoms laid along x: P, O5', C5', C3', O3' (spacing step Angstrom)."""
    names = (['P'] if phosphate else []) + ["O5'", "C5'"] + (["C3'", "O3'"] if o3 else [])
    return {name: (x0 + i*step, y, z) for i, name in enumerate(names)}


def cgd_rows(text):
    atoms, edges = {}, {}
    for line in text.splitlines():
        parts = line.split()
        if parts and parts[0] == 'atom':
            atoms[parts[1]] = (int(parts[2]), [float(v) for v in parts[3:]])
            edges[parts[1]] = []
        elif parts and parts[0] == 'edge':
            edges[parts[1]].append([float(v) for v in parts[2:]])
    return atoms, edges


def assert_reciprocal(test, text, group):
    """Each edge row must have its inverse-mapped reciprocal row (any group image)."""
    atoms, edges = cgd_rows(text)
    ops = list(gemmi.find_spacegroup_by_name(group).operations())

    def same(a, b):
        return all(abs(x-y) < 1e-6 for x, y in zip(a, b))
    for source, rows in edges.items():
        test.assertEqual(atoms[source][0], len(rows), f'declared degree of {source}')
        for row in rows:
            matches = []
            for target, (_, frac) in atoms.items():
                for op in ops:
                    image = op.apply_to_xyz(frac)
                    shift = [round(r-i) for r, i in zip(row, image)]
                    if same([i+s for i, s in zip(image, shift)], row):
                        matches.append((target, op, shift))
            test.assertTrue(matches, f'edge row {source} {row} names no node image')
            target, op, shift = matches[0]
            g = gemmi.Op(op.triplet())
            g.tran = [t + 24*s for t, s in zip(g.tran, shift)]
            back = g.inverse().apply_to_xyz(atoms[source][1])
            test.assertTrue(any(same(back, r) for r in edges[target]),
                            f'missing reciprocal of {source}->{target}: expected {back}')


class StitchingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'input.pdb'
        self.out = self.root / 'graph.cgd'

    def export(self, text, **kw):
        self.source.write_text(text)
        return export_cgd(self.source, self.out, **kw)

    def translation_chain(self):
        # Residue 2 O3' (x=18) is 2.5 A from residue 1 P translated by +a (x=20.5);
        # residue 1 O3' (x=4.5) is bonded to residue 2 P (x=6, 1.5 A).
        second = {'P': (6., 5., 5.), "O5'": (7., 5., 5.), "C5'": (8., 5., 5.),
                  "C3'": (17., 5., 5.), "O3'": (18., 5., 5.)}
        return (crystal() + nucleotide(1, 1, 'A', line_residue(.5))
                + nucleotide(10, 2, 'A', second))

    def test_default_export_remains_unstitched_schema_1(self):
        r = self.export(self.translation_chain())
        self.assertEqual((r['schema_version'], r['stitching'], r['edge_count']), (1, False, 1))
        self.assertNotIn('stitched_edges', r)

    def test_lattice_translation_nick_closes_periodic_strand(self):
        r = self.export(self.translation_chain(), stitch_nicks=True)
        self.assertEqual((r['schema_version'], r['stitching']), (2, True))
        self.assertEqual((r['edge_count'], r['trace_edge_count'], r['stitched_edge_count']), (2, 1, 1))
        self.assertEqual([n['degree'] for n in r['nodes']], [2, 2])
        edge, = r['stitched_edges']
        self.assertEqual((edge['source'], edge['target'], edge['kind']), ('S1_N2', 'S1_N1', 'backbone_join'))
        self.assertEqual((edge['operation'], edge['lattice_translation']), ('x,y,z', [1, 0, 0]))
        self.assertEqual(edge['image_operation'], 'x+1,y,z')
        self.assertAlmostEqual(edge['distance_A'], 2.5, places=6)
        self.assertEqual([round(v, 9) for v in edge['target_image_fractional']], [1.225, round(1/6, 9), .125])
        self.assertEqual([round(v, 9) for v in edge['reverse_image_fractional']], [-.1, round(1/6, 9), .125])
        text = self.out.read_text()
        self.assertIn('edge S1_N2 1.225 ', text)
        self.assertIn('edge S1_N1 -0.1 ', text)
        assert_reciprocal(self, text, 'P1')
        candidate, = r['join_candidates']
        self.assertEqual(candidate['status'], 'stitched')
        self.assertEqual((candidate['o3_atom']['residue'], candidate['p_atom']['seq']), ('DG', '1'))

    def test_screw_operation_negative_translation_and_inverse(self):
        # P 1 21 1: chain B P image under -x+1,y-1/2,-z+1 lies 2.5 A from chain A O3'.
        a = {"C5'": (4.5, 12., 12.), "C3'": (5., 12., 12.5), "O3'": (6., 12., 12.)}
        b = {'P': (11.5, 27., 28.), "O5'": (11.5, 27., 26.6), "C5'": (11.5, 25.6, 26.6),
             "C3'": (12.5, 24.6, 26.6), "O3'": (13., 23.6, 26.6)}
        r = self.export(crystal(group='P 1 21 1') + nucleotide(1, 1, 'A', a)
                        + 'TER\n' + nucleotide(10, 1, 'B', b), stitch_nicks=True)
        edge, = r['stitched_edges']
        self.assertEqual((edge['source'], edge['target']), ('S1_N1', 'S2_N1'))
        self.assertEqual(edge['operation'], '-x,y+1/2,-z')
        self.assertEqual(edge['lattice_translation'], [1, -1, 1])
        self.assertEqual(edge['image_operation'], '-x+1,y-1/2,-z+1')
        self.assertEqual(edge['reverse_image_operation'], '-x+1,y+1/2,-z+1')
        src = r['nodes'][0]['fractional']
        expected = gemmi.Op('-x+1,y+1/2,-z+1').apply_to_xyz(src)
        for got, want in zip(edge['reverse_image_fractional'], expected):
            self.assertAlmostEqual(got, want, places=9)
        self.assertEqual([n['degree'] for n in r['nodes']], [1, 1])
        assert_reciprocal(self, self.out.read_text(), 'P 1 21 1')

    def test_ambiguous_partners_are_reported_not_stitched(self):
        # Two different 5' ends lie within 3.5 A of one O3' end.
        text = (crystal() + nucleotide(1, 1, 'A', line_residue(2.))
                + 'TER\n' + nucleotide(10, 1, 'B', line_residue(8.5))
                + 'TER\n' + nucleotide(20, 1, 'C', line_residue(5., 7.5, 5.)))
        r = self.export(text, stitch_nicks=True)
        self.assertEqual(r['stitched_edge_count'], 0)
        statuses = sorted(c['status'] for c in r['join_candidates'])
        self.assertEqual(statuses, ['ambiguous']*3)
        self.assertTrue(any('ambiguous' in w for w in r['warnings']))

    def test_same_chain_copy_and_linked_ends_are_not_candidates(self):
        # Last O3' 2.5 A from first P in the same image: same chain copy (NAtopo rule).
        last = {'P': (21., 5., 5.), "O5'": (22., 6., 5.), "C5'": (22., 7., 5.),
                "C3'": (12., 7.5, 5.), "O3'": (11.5, 7.5, 5.)}
        circular = (crystal(60, 60, 60) + nucleotide(1, 1, 'A', line_residue(10.))
                    + nucleotide(10, 2, 'A', line_residue(15.5)) + nucleotide(20, 3, 'A', last))
        r = self.export(circular, stitch_nicks=True)
        self.assertEqual(r['join_candidates'], [])
        self.assertEqual(r['stitched_edge_count'], 0)

    def test_segment_internal_break_is_reported(self):
        # Residues 1 and 2 share a segment but are not covalently linked
        # (O3'-P 3.5 A > 2.2 A), so both ends are loose yet internal to the trace.
        text = (crystal(60, 60, 60) + nucleotide(1, 1, 'A', line_residue(10.))
                + nucleotide(10, 2, 'A', line_residue(17.5))
                + 'TER\n' + nucleotide(20, 1, 'B', line_residue(14., 7.5, 5.)))
        r = self.export(text, stitch_nicks=True)
        self.assertEqual(r['stitched_edge_count'], 0)
        self.assertIn('segment_internal', {c['status'] for c in r['join_candidates']})

    def test_o5_fallback_is_opt_in_and_maps_missing_o3_residue(self):
        # Chain B residue 1 has no P and no O3'; residue 2 holds the first B node.
        b1 = line_residue(9.5, phosphate=False, o3=False)
        text = (crystal(60, 60, 60) + nucleotide(1, 1, 'A', line_residue(2.))
                + 'TER\n' + nucleotide(10, 1, 'B', b1)
                + nucleotide(20, 2, 'B', line_residue(12., 9., 5.)))
        r = self.export(text, stitch_nicks=True)
        candidate, = r['join_candidates']
        self.assertEqual((candidate['kind'], candidate['status']), ('backbone_join_o5', 'o5_fallback_disabled'))
        self.assertEqual(r['stitched_edge_count'], 0)
        self.out.unlink()
        r = self.export(text, stitch_nicks=True, include_o5_fallback=True)
        edge, = r['stitched_edges']
        self.assertEqual((edge['source'], edge['target'], edge['kind']), ('S1_N1', 'S2_N1', 'backbone_join_o5'))
        self.assertAlmostEqual(edge['distance_A'], 3.5)

    def test_cutoff_and_explicit_links(self):
        r = self.export(self.translation_chain(), stitch_nicks=True, join_cutoff=2.4)
        self.assertEqual((r['join_candidates'], r['stitched_edge_count']), ([], 0))
        self.out.unlink()
        # A deposited O3'-P link between two segments makes both ends non-loose
        # and is stitched as a same-image deposited link.
        text = (crystal(60, 60, 60) + nucleotide(1, 1, 'A', line_residue(10.))
                + 'TER\n' + nucleotide(10, 1, 'B', line_residue(15.6))
                + 'CONECT    5   10\n')
        r = self.export(text, stitch_nicks=True)
        edge, = r['stitched_edges']
        self.assertEqual((edge['kind'], edge['image_operation']), ('deposited_link', 'x,y,z'))
        self.assertEqual([n['degree'] for n in r['nodes']], [1, 1])

    def test_r3_centering_operations(self):
        # O3' of chain A is 2.5 A from the P of chain B after x+2/3,y+1/3,z+1/3.
        cell = gemmi.UnitCell(100, 100, 90, 90, 90, 120)
        o3 = gemmi.Fractional(.10, .20, .30)
        p = gemmi.Fractional(.10-2/3+.025, .20-1/3, .30-1/3)
        def orth(f, dx=0., dy=0., dz=0.):
            v = cell.orthogonalize(f)
            return (v.x+dx, v.y+dy, v.z+dz)
        a = {"C3'": orth(o3, -1.2), "O3'": orth(o3)}
        pp = orth(p)
        b = {'P': pp, "O5'": (pp[0], pp[1]+1.5, pp[2]), "C5'": (pp[0], pp[1]+3., pp[2]),
             "C3'": (pp[0]-1, pp[1]+4, pp[2]), "O3'": (pp[0]-2, pp[1]+5, pp[2])}
        cmm = self.root/'copies.cmm'
        r = self.export(crystal(100, 100, 90, 90, 90, 120, 'H 3') + nucleotide(1, 1, 'A', a)
                        + 'TER\n' + nucleotide(10, 1, 'B', b), stitch_nicks=True, crystal_cmm=cmm)
        edge, = r['stitched_edges']
        self.assertEqual(edge['operation'], 'x+2/3,y+1/3,z+1/3')
        self.assertEqual(edge['lattice_translation'], [0, 0, 0])
        assert_reciprocal(self, self.out.read_text(), 'R3:H')
        joins = [m for m in ET.parse(cmm).getroot() if 'nick joins' in m.get('name')]
        self.assertEqual(len(joins), 1)
        self.assertEqual(len(joins[0].findall('link')), 9)

    def test_asu_cmm_draws_join_to_image_marker(self):
        cmm = self.root/'graph.cmm'
        r = self.export(self.translation_chain(), stitch_nicks=True, cmm=cmm)
        root = ET.parse(cmm).getroot()
        joins = [m for m in root if 'nick joins' in m.get('name')]
        self.assertEqual(len(joins), 1)
        markers = joins[0].findall('marker')
        self.assertEqual(len(joins[0].findall('link')), 1)
        self.assertEqual(sorted(round(float(m.get('x')), 6) for m in markers), [18., 24.5])
        self.assertEqual([round(v, 6) for v in r['stitched_edges'][0]['target_image_cartesian']], [24.5, 5., 5.])

    def test_cli_and_receipt_schema_2(self):
        self.source.write_text(self.translation_chain())
        receipt = self.root/'receipt.json'
        p = subprocess.run([sys.executable, '-m', 'topostich.cgd', str(self.source), '--out', str(self.out),
                            '--stitch-nicks', '--join-cutoff', '3.5', '--receipt', str(receipt)],
                           text=True, capture_output=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn('1 stitched nick', p.stdout)
        r = json.loads(receipt.read_text())
        self.assertEqual(r['schema_version'], 2)
        try:
            import jsonschema
        except ImportError:
            self.skipTest('jsonschema not installed')
        jsonschema.validate(r, json.loads((SCHEMAS/'cgd-receipt-2.schema.json').read_text()))
        self.out.unlink(); receipt.unlink()
        export_cgd(self.source, self.out, receipt=receipt)
        jsonschema.validate(json.loads(receipt.read_text()),
                            json.loads((SCHEMAS/'cgd-receipt.schema.json').read_text()))

    def test_invalid_parameters(self):
        for kw in ({'join_cutoff': 0}, {'o5_cutoff': -1}, {'join_cutoff': float('nan')}):
            with self.subTest(kw=kw), self.assertRaises(ConversionError):
                self.export(self.translation_chain(), stitch_nicks=True, **kw)
            self.assertFalse(self.out.exists())


if __name__ == '__main__':
    unittest.main()
