"""Junction (non-backbone pi-stack) stitching regressions.

The rule is NAtopo's nonbackbone_pi_stack_v2 over its find_pi_interactions
screen: an accepted base stack not already joined by a verified O3'-P bond.
'face' scope stitches only stacks on a base face without backbone continuation.
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


def residue(cx, cy, z, *, backbone=(), dx=5.):
    """A cytosine ring centered at (cx, cy, z) plus backbone atoms beside it.

    backbone maps atom name -> z offset; those atoms sit at x = cx + dx.
    """
    atoms = {}
    for k, name in enumerate(('N1', 'C2', 'N3', 'C4', 'C5', 'C6')):
        angle = math.radians(60*k)
        atoms[name] = (cx + 1.39*math.cos(angle), cy + 1.39*math.sin(angle), z)
    for name, dz in backbone:
        atoms[name] = (cx + dx, cy, z + dz)
    return atoms


# 3' side of a residue (C3', O3') and 5' side (P, O5', C5') laid along z.
THREE = (("C3'", .2), ("O3'", .8))
FIVE = (('P', -.8), ("O5'", -.2), ("C5'", .4))


class JunctionTests(unittest.TestCase):
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
        return export_cgd(self.source, self.out, **kw)

    def column(self):
        # A1 (z=1.7) -> A2 (z=5.1) bonded O3'-P 1.6 A; c = 6.8 A, so A2's free top
        # face stacks on A1 one cell up: a continuous stacked column.
        a1 = residue(10, 10, 1.7, backbone=(("O5'", -.4), ("C5'", -.2), ("C3'", .2), ("O3'", .8)))
        a2 = residue(10, 10, 5.1, backbone=(('P', -1.), ("O5'", -.6), ("C5'", -.4), ("C3'", .2), ("O3'", .4)))
        return crystal(20, 20, 6.8) + nucleotide(1, 1, 'A', a1, 'DC') + nucleotide(20, 2, 'A', a2, 'DC')

    def test_face_junction_across_lattice_translation(self):
        r = self.export(self.column(), stitch_junctions=True)
        self.assertEqual((r['schema_version'], r['stitching']), (2, True))
        edge, = r['stitched_edges']
        self.assertEqual((edge['kind'], edge['source'], edge['target']), ('junction_stack', 'S1_N1', 'S1_N2'))
        self.assertEqual(edge['image_operation'], 'x,y,z-1')
        self.assertAlmostEqual(edge['distance_A'], 3.4, places=6)
        self.assertEqual(edge['evidence']['reason'], 'stack_on_face_without_backbone_continuation')
        self.assertEqual([n['degree'] for n in r['nodes']], [2, 2])
        candidate, = r['junction_candidates']
        self.assertEqual(candidate['status'], 'stitched')
        self.assertEqual(r['stitch_parameters']['junctions']['verified_backbone_stacks_excluded'], 1)
        self.assertEqual(r['join_candidates'], [])
        assert_reciprocal(self, self.out.read_text(), 'P1')

    def test_face_scope_skips_stacks_on_backbone_faces(self):
        # A1-A2-A3 stack along z; D1 stacks on A2's top face, which A3 occupies.
        a1 = residue(10, 10, 10, backbone=THREE)
        a2 = residue(10, 10, 13.4, backbone=FIVE + THREE)
        a3 = residue(10, 10, 16.8, backbone=FIVE)
        d1 = residue(12.5, 10, 16.8, backbone=FIVE + THREE, dx=-12.)
        text = (crystal(60, 60, 60) + nucleotide(1, 1, 'A', a1, 'DC') + nucleotide(20, 2, 'A', a2, 'DC')
                + nucleotide(40, 3, 'A', a3, 'DC') + 'TER\n' + nucleotide(60, 1, 'D', d1, 'DC'))
        r = self.export(text, stitch_junctions=True)
        candidate, = r['junction_candidates']
        self.assertEqual((candidate['reason'], candidate['status']), ('non_backbone_stack', 'not_face_junction'))
        self.assertEqual(r['stitched_edge_count'], 0)
        self.assertEqual(r['stitch_parameters']['junctions']['verified_backbone_stacks_excluded'], 2)
        r = self.export(text, stitch_junctions=True, junction_scope='all')
        edge, = r['stitched_edges']
        self.assertEqual((edge['source'], edge['target'], edge['image_operation']), ('S1_N2', 'S2_N1', 'x,y,z'))
        self.assertAlmostEqual(edge['evidence']['lateral_shift_A'][0], 2.5, places=6)

    def test_stack_geometry_limits(self):
        # Base B1 above A1: 3.4 A passes; tilted, too far, or too shifted fails.
        a1 = residue(10, 10, 10, backbone=THREE)
        for dz, dx, ok in ((3.4, 0, True), (4.6, 0, False), (3.4, 3.1, False), (1.9, 0, False)):
            b1 = residue(10+dx, 10, 10+dz, backbone=FIVE, dx=-10.)
            with self.subTest(dz=dz, dx=dx):
                r = self.export(crystal(60, 60, 60) + nucleotide(1, 1, 'A', a1, 'DC') + 'TER\n'
                                + nucleotide(20, 1, 'B', b1, 'DC'), stitch_junctions=True, junction_scope='all')
                self.assertEqual(len(r['junction_candidates']), int(ok))

    def test_nick_and_stack_on_same_pair_are_not_duplicated(self):
        # A1 3' end stacks on B1 5' end and its O3' is 2.5 A from B1's P.
        a1 = residue(10, 10, 10, backbone=(("C3'", .2), ("O3'", .9)))
        b1 = residue(10, 10, 13.4, backbone=(('P', 0.), ("O5'", -.2), ("C5'", .4), ("C3'", 1.), ("O3'", 1.6)))
        text = crystal(60, 60, 60) + nucleotide(1, 1, 'A', a1, 'DC') + 'TER\n' + nucleotide(20, 1, 'B', b1, 'DC')
        r = self.export(text, stitch_nicks=True, stitch_junctions=True, junction_scope='all')
        self.assertEqual([e['kind'] for e in r['stitched_edges']], ['backbone_join'])
        self.assertEqual(r['junction_candidates'][0]['status'], 'duplicate_edge')
        self.assertEqual([n['degree'] for n in r['nodes']], [1, 1])

    def test_unmapped_residue_has_no_plane_and_is_reported(self):
        a1 = residue(10, 10, 10, backbone=THREE)
        b1 = residue(10, 10, 13.4, backbone=FIVE, dx=-10.)
        r = self.export(crystal(60, 60, 60) + nucleotide(1, 1, 'A', a1, 'DC') + 'TER\n'
                        + nucleotide(20, 1, 'B', b1, 'XQZ'), stitch_junctions=True, junction_scope='all')
        self.assertEqual(r['junction_candidates'], [])
        issues = r['stitch_parameters']['junctions']['plane_diagnostics']
        self.assertEqual([(d['residue'], d['issue']) for d in issues], [('XQZ', 'unmapped_component')])

    def test_cmm_draws_junction_model(self):
        cmm, copies = self.root/'graph.cmm', self.root/'copies.cmm'
        self.export(self.column(), stitch_junctions=True, cmm=cmm, crystal_cmm=copies)
        sets = [m.get('name') for m in ET.parse(cmm).getroot()]
        self.assertTrue(any(name.endswith('junction stacks') for name in sets))
        self.assertFalse(any('nick joins' in name for name in sets))
        junction = [m for m in ET.parse(copies).getroot() if 'junction stacks' in m.get('name')][0]
        self.assertEqual(len(junction.findall('link')), 1)

    def test_stack_on_own_twofold_mate_is_one_neighbor(self):
        # Terminal A1 ring (plane x=1.7) stacks on its own image under -x,y,-z.
        ring = {name: (1.7, 10 + 1.39*math.cos(math.radians(60*k)), 1.39*math.sin(math.radians(60*k)))
                for k, name in enumerate(('N1', 'C2', 'N3', 'C4', 'C5', 'C6'))}
        a1 = dict(ring, **{"C3'": (1.9, 16, 3), "O3'": (2.5, 16, 3)})
        a2 = {name: (x+3.4, y, z) for name, (x, y, z) in ring.items()}
        a2.update({'P': (4.1, 16, 3), "O5'": (4.5, 16, 3), "C5'": (4.7, 16, 3), "C3'": (5.3, 16, 3), "O3'": (5.5, 16, 3)})
        r = self.export(crystal(30, 30, 30, group='P 1 2 1') + nucleotide(1, 1, 'A', a1, 'DC')
                        + nucleotide(20, 2, 'A', a2, 'DC'), stitch_junctions=True)
        edge, = r['stitched_edges']
        self.assertEqual((edge['source'], edge['target'], edge['image_operation']), ('S1_N1', 'S1_N1', '-x,y,-z'))
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
        r = self.export(self.column(), stitch_junctions=True)
        records = [['NODES', '2'], ['EDGE', 'S1_N1', 'S1_N2', 'Vector(0,0,0)'],
                   ['EDGE', 'S1_N2', 'S1_N1', 'Vector(0,0,1)'], ['COMPONENT', '1', '2', '1']]
        problems, summary = audit.compare(records, r)
        self.assertEqual(problems, [])
        self.assertEqual(summary['components'], {'dimension 1, 2 nodes': 1})
        problems, _ = audit.compare(records[:2] + records[3:], r)
        self.assertEqual(problems, ['S1_N1-S1_N2: Systre 1 edges, receipt implies 2'])


if __name__ == '__main__':
    unittest.main()
