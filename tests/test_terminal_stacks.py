"""Chain-end pi-stack stitching regressions.

Every pi stack between two chain ends (first or last nucleotide of a segment)
is stitched, in any symmetry image and for any polarity, unless the two ends
are already backbone joined. Same-polarity (5'-5', 3'-3') stacks are flagged.
"""
import importlib.util
import math
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

from test_converter import crystal
from test_stitching import assert_reciprocal, nucleotide
from topostich.cgd import export_cgd

TOOLS = Path(__file__).resolve().parents[1] / 'tools'


def ring(cx, cy, z):
    """Cytosine ring atoms centered at (cx, cy, z) in the z plane."""
    return {name: (cx + 1.39*math.cos(math.radians(60*k)), cy + 1.39*math.sin(math.radians(60*k)), z)
            for k, name in enumerate(('N1', 'C2', 'N3', 'C4', 'C5', 'C6'))}


def strand(serial, chain, zs, *, cx=10., cy=10., dx=5., phosphate_5=False, drop_p=()):
    """A strand stacked along z at heights zs (5'->3'); backbone beside it at x = cx + dx.

    Consecutive residues are bonded O3'(i)-P(i+1) at 1.8 A. The 5' residue has
    no phosphate unless phosphate_5; residues in drop_p (1-based) lose their P.
    """
    text = ''
    for i, z in enumerate(zs):
        d = 1 if zs[-1] >= zs[0] else -1
        atoms = ring(cx, cy, z)
        backbone = ([('P', -.8)] if (i or phosphate_5) and i+1 not in drop_p else []) + \
            [("O5'", -.2), ("C5'", .4), ("C3'", .2), ("O3'", .8)]
        for name, dz in backbone:
            atoms[name] = (cx + dx, cy, z + d*dz)
        text += nucleotide(serial + 20*i, i+1, chain, atoms, 'DC')
    return text + 'TER\n'


class TerminalStackTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'input.pdb'
        self.out = self.root / 'graph.cgd'

    def export(self, text, **kw):
        self.source.write_text(text)
        if self.out.exists():
            self.out.unlink()
        kw.setdefault('stitch_terminal_stacks', True)
        return export_cgd(self.source, self.out, **kw)

    def test_three_five_stack_across_lattice_translation(self):
        # A1 (5') -> A2 (3'), c = 6.8 A: A2's top face stacks on A1 one cell up.
        r = self.export(crystal(20, 20, 6.8) + strand(1, 'A', (1.7, 5.1)))
        self.assertEqual((r['schema_version'], r['stitching']), (2, True))
        edge, = r['stitched_edges']
        self.assertEqual((edge['kind'], edge['source'], edge['target']), ('terminal_stack', 'S1_N2', 'S1_N1'))
        self.assertEqual(edge['image_operation'], 'x,y,z+1')
        self.assertAlmostEqual(edge['distance_A'], 3.4, places=6)
        self.assertEqual({k: edge['evidence'][k] for k in ('polarity', 'same_polarity', 'source_end', 'target_end')},
                         {'polarity': "5'-3'", 'same_polarity': False, 'source_end': "3'", 'target_end': "5'"})
        self.assertEqual([n['degree'] for n in r['nodes']], [2, 2])
        statuses = sorted((c['status'], c['image_operation']) for c in r['terminal_stack_candidates'])
        self.assertEqual(statuses, [('backbone_joined', 'x,y,z'), ('stitched', 'x,y,z+1')])
        self.assertEqual(r['join_candidates'], [])
        assert_reciprocal(self, self.out.read_text(), 'P1')

    def test_same_polarity_stacks_are_flagged(self):
        # 3'-3': A runs up to z=13.4, B runs down to z=16.8; their 3' ends stack.
        three = (crystal(60, 60, 60) + strand(1, 'A', (10, 13.4)) + strand(100, 'B', (20.2, 16.8), dx=-5))
        # 5'-5': A runs down from z=13.4, B runs up from z=16.8; their 5' ends stack.
        five = (crystal(60, 60, 60) + strand(1, 'A', (13.4, 10)) + strand(100, 'B', (16.8, 20.2), dx=-5))
        for text, polarity, end in ((three, "3'-3'", "3'"), (five, "5'-5'", "5'")):
            with self.subTest(polarity=polarity):
                r = self.export(text)
                edge, = r['stitched_edges']
                self.assertEqual((edge['evidence']['polarity'], edge['evidence']['same_polarity']), (polarity, True))
                self.assertEqual((edge['evidence']['source_end'], edge['evidence']['target_end']), (end, end))
                self.assertEqual(edge['image_operation'], 'x,y,z')
                summary = r['stitch_parameters']['terminal_stacks']
                self.assertEqual(summary['stitched_polarity_counts'], {polarity: 1})
                self.assertEqual(summary['candidate_status_counts'], {'backbone_joined': 2, 'stitched': 1})
                self.assertTrue(any("same-polarity" in w for w in r['warnings']))

    def test_internal_residues_are_never_candidates(self):
        # A2 (internal) stacks on A3 (3' end) with no O3'-P bond: not a chain-end pair.
        r = self.export(crystal(60, 60, 60) + strand(1, 'A', (10, 13.4, 16.8), drop_p=(3,)))
        self.assertEqual(r['terminal_stack_candidates'], [])
        self.assertEqual(r['stitch_parameters']['terminal_stacks']['chain_ends'], 2)

    def test_single_residue_segments_are_counted_not_stitched(self):
        # A2 (3' end) stacks on lone residue B1, which is not a chain end.
        r = self.export(crystal(60, 60, 60) + strand(1, 'A', (6.6, 10)) + strand(100, 'B', (13.4,), dx=-5))
        self.assertEqual([c['status'] for c in r['terminal_stack_candidates']], ['backbone_joined'])
        self.assertEqual(r['stitched_edges'], [])
        self.assertEqual(r['stitch_parameters']['terminal_stacks']['single_residue_segments'], 1)

    def test_nick_joined_ends_are_not_stacked_again(self):
        # A2 (3') stacks on B1 (5'); B1's P is 1.8 A from A2's O3' across segments.
        text = crystal(60, 60, 60) + strand(1, 'A', (6.6, 10)) + strand(100, 'B', (13.4, 16.8), phosphate_5=True)
        r = self.export(text, stitch_nicks=True)
        self.assertEqual([e['kind'] for e in r['stitched_edges']], ['backbone_join'])
        across = [c for c in r['terminal_stack_candidates'] if c['residue1']['chain'] != c['residue2']['chain']]
        self.assertEqual([c['status'] for c in across], ['backbone_joined'])
        r = self.export(text)
        across = [e for e in r['stitched_edges'] if e['source'][:2] != e['target'][:2]]
        self.assertEqual([(e['kind'], e['evidence']['polarity']) for e in across], [('terminal_stack', "5'-3'")])

    def test_stack_geometry_limits(self):
        for dz, dx, ok in ((3.4, 0, True), (4.6, 0, False), (3.4, 3.1, False), (1.9, 0, False)):
            with self.subTest(dz=dz, dx=dx):
                r = self.export(crystal(60, 60, 60) + strand(1, 'A', (6.6, 10))
                                + strand(100, 'B', (10+dz, 13.4+dz), cx=10+dx, dx=-10))
                across = [c for c in r['terminal_stack_candidates'] if c['residue1']['chain'] != c['residue2']['chain']]
                self.assertEqual(len(across), int(ok))

    def test_unmapped_end_has_no_plane_and_is_reported(self):
        text = (crystal(60, 60, 60) + strand(1, 'A', (6.6, 10))
                + strand(100, 'B', (13.4, 16.8), dx=-5).replace(' DC B   1', 'XQZ B   1'))
        r = self.export(text)
        self.assertFalse([c for c in r['terminal_stack_candidates'] if c['residue1']['chain'] != c['residue2']['chain']])
        issues = r['stitch_parameters']['terminal_stacks']['plane_diagnostics']
        self.assertEqual([(d['residue'], d['issue']) for d in issues], [('XQZ', 'unmapped_component')])

    def test_cmm_draws_terminal_stack_model(self):
        cmm, copies = self.root/'graph.cmm', self.root/'copies.cmm'
        self.export(crystal(20, 20, 6.8) + strand(1, 'A', (1.7, 5.1)), cmm=cmm, crystal_cmm=copies)
        sets = [m.get('name') for m in ET.parse(cmm).getroot()]
        self.assertTrue(any(name.endswith('terminal stacks') for name in sets))
        self.assertFalse(any('nick joins' in name for name in sets))
        stacks = [m for m in ET.parse(copies).getroot() if 'terminal stacks' in m.get('name')][0]
        self.assertEqual(len(stacks.findall('link')), 1)

    def twofold(self):
        # 5' end A1 (ring plane x=1.7) stacks on its own image under -x,y,-z.
        a1 = {name: (1.7, 10 + 1.39*math.cos(math.radians(60*k)), 1.39*math.sin(math.radians(60*k)))
              for k, name in enumerate(('N1', 'C2', 'N3', 'C4', 'C5', 'C6'))}
        a2 = {name: (x+3.4, y, z) for name, (x, y, z) in a1.items()}
        a1.update({"O5'": (1.5, 16, 3), "C5'": (1.6, 16, 3), "C3'": (1.9, 16, 3), "O3'": (2.5, 16, 3)})
        a2.update({'P': (4.1, 16, 3), "O5'": (4.5, 16, 3), "C5'": (4.7, 16, 3), "C3'": (5.3, 16, 3), "O3'": (5.5, 16, 3)})
        return crystal(30, 30, 30, group='P 1 2 1') + nucleotide(1, 1, 'A', a1, 'DC') + nucleotide(20, 2, 'A', a2, 'DC')

    def test_stack_on_own_twofold_mate_is_one_neighbor(self):
        r = self.export(self.twofold())
        edge, = r['stitched_edges']
        self.assertEqual((edge['source'], edge['target'], edge['image_operation']), ('S1_N1', 'S1_N1', '-x,y,-z'))
        self.assertEqual(edge['evidence']['polarity'], "5'-5'")
        self.assertEqual(edge['reverse_image_operation'], edge['image_operation'])
        self.assertEqual([n['degree'] for n in r['nodes']], [2, 1])
        rows = [line for line in self.out.read_text().splitlines() if line.startswith('edge S1_N1 ')]
        self.assertEqual(len(rows), 2)
        assert_reciprocal(self, self.out.read_text(), 'P 1 2 1')
        spec = importlib.util.spec_from_file_location('check_stitched_cgd', TOOLS/'check_stitched_cgd.py')
        audit = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(audit)
        records = [['NODES', '4'], ['EDGE', 'S1_N1', 'S1_N2', 'Vector(0,0,0)'],
                   ['EDGE', 'S1_N1', 'S1_N2', 'Vector(0,0,0)'], ['EDGE', 'S1_N1', 'S1_N1', 'Vector(0,0,0)']]
        self.assertEqual(audit.compare(records, r)[0], [])

    def test_systre_audit_comparison_logic(self):
        spec = importlib.util.spec_from_file_location('check_stitched_cgd', TOOLS/'check_stitched_cgd.py')
        audit = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(audit)
        r = self.export(crystal(20, 20, 6.8) + strand(1, 'A', (1.7, 5.1)))
        records = [['NODES', '2'], ['EDGE', 'S1_N1', 'S1_N2', 'Vector(0,0,0)'],
                   ['EDGE', 'S1_N2', 'S1_N1', 'Vector(0,0,1)'], ['COMPONENT', '1', '2', '1']]
        problems, summary = audit.compare(records, r)
        self.assertEqual(problems, [])
        self.assertEqual(summary['components'], {'dimension 1, 2 nodes': 1})
        problems, _ = audit.compare(records[:2] + records[3:], r)
        self.assertEqual(problems, ['S1_N1-S1_N2: Systre 1 edges, receipt implies 2'])


if __name__ == '__main__':
    unittest.main()
