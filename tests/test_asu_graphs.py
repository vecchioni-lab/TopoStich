"""ASU quotient graphs (chain and residue resolution) of covalent bonds and stitched joins.

Ported from NAtopo's asu_quotient/residue_quotient; the symmetry algebra and
periodicity are checked directly, then through export_cgd.
"""
import json
from pathlib import Path
import tempfile
import unittest

import gemmi

from test_converter import crystal
from test_stitching import SCHEMAS, nucleotide
from test_terminal_stacks import strand
from topostich import asu_graphs as q
from topostich.cgd import export_cgd
from topostich.converter import ConversionError

I = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


def op(triplet):
    return q.from_gemmi(gemmi.Op(triplet))


class AlgebraTests(unittest.TestCase):
    def test_compose_inverse_and_codes(self):
        screw = op('-x,y+1/2,-z')
        self.assertEqual(q.compose(screw, screw), (I, (0, 24, 0)))
        self.assertEqual(q.compose(screw, q.inverse(screw)), q.IDENTITY)
        ops = [(1, op('x,y,z')), (2, op('-x,y+1/2,-z'))]
        self.assertEqual(q.symmetry_code(op('x,y,z+1'), ops), ('1_556', '1', [0, 0, 1]))
        self.assertEqual(q.symmetry_code(op('-x+1,y-1/2,-z'), ops), ('2_645', '2', [1, -1, 0]))
        self.assertEqual(q.symmetry_code(op('x+5,y,z'), ops)[0], '1(5, 0, 0)')
        with self.assertRaises(ValueError):
            q.symmetry_code(op('-x,-y,z'), ops)

    def test_generated_translations(self):
        basis, order = q.generated_translations([op('-x,y+1/2,-z')])
        self.assertEqual((basis, order), ([[0, 24, 0]], 2))
        basis, order = q.generated_translations([op('-x,-y,z')])  # finite 2-fold cluster
        self.assertEqual((basis, order), ([], 2))
        self.assertEqual(q._hermite_basis([(48, 0, 0), (24, 24, 0), (0, 0, 24)]), [[24, 24, 0], [0, 48, 0], [0, 0, 24]])

    def test_periodicity_and_interpenetration(self):
        ops = [(1, q.IDENTITY)]
        loops = [dict(source='a', target='a', operation=op(t)) for t in ('x+2,y,z', 'x,y+1,z', 'x,y,z+1')]
        part, = q.periodicity(['a'], loops, ops)
        self.assertEqual((part['dimension'], part['interpenetrating_copies']), (3, 2))
        self.assertEqual(part['translation_basis'], ['[2, 0, 0]', '[0, 1, 0]', '[0, 0, 1]'])
        # Two nodes joined twice: the cycle carries the translation.
        edges = [dict(source='a', target='b', operation=q.IDENTITY),
                 dict(source='b', target='a', operation=op('x,y,z-1'))]
        part, single = q.periodicity(['a', 'b', 'c'], edges, ops)
        self.assertEqual((part['nodes'], part['dimension'], part['translation_basis']), (['a', 'b'], 1, ['[0, 0, 1]']))
        self.assertEqual((single['nodes'], single['dimension']), (['c'], 0))

    def test_numbering_follows_input_operators(self):
        class Fake:
            symops = {2: '-x,y+1/2,-z', 1: 'x,y,z'}
        self.assertEqual([i for i, _ in q.numbered_operations(Fake, [])], [1, 2])
        Fake.symops = {}
        ops = list(gemmi.find_spacegroup_by_name('P 1 21 1').operations())
        self.assertEqual([i for i, _ in q.numbered_operations(Fake, ops)], [1, 2])


class ExportTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'input.pdb'
        self.schema = json.loads((SCHEMAS / 'asu-graphs.schema.json').read_text())

    def export(self, text, **kw):
        if text is not None:
            self.source.write_text(text)
        out = self.root / 'out'
        if out.exists():
            for f in out.iterdir():
                f.unlink()
        else:
            out.mkdir()
        r = export_cgd(self.source, out/'g.cgd', asu_graphs=out/'asu.json', asu_chain_html=out/'chain.html',
                       asu_residue_html=out/'residue.html', **kw)
        document = json.loads((out/'asu.json').read_text())
        try:
            import jsonschema
        except ImportError:
            pass
        else:
            jsonschema.validate(document, self.schema)
        return r, document, out

    def nick_chain(self):
        # Residue 2 O3' meets residue 1 P one cell along a (2.5 A); 1 -> 2 is a polymer link.
        atoms = [('P', 1, .5), ("O5'", 1, 1.5), ("C5'", 1, 2.5), ("C3'", 1, 3.5), ("O3'", 1, 4.5),
                 ('P', 2, 6.), ("O5'", 2, 7.), ("C5'", 2, 8.), ("C3'", 2, 17.), ("O3'", 2, 18.)]
        return crystal(20, 30, 40) + ''.join(
            nucleotide(i, res, 'A', {name: (x, 5., 5.)}) for i, (name, res, x) in enumerate(atoms, 1))

    def test_nick_self_join_makes_a_rod(self):
        r, doc, out = self.export(self.nick_chain(), stitch_nicks=True)
        chain = doc['chain_graph']
        self.assertEqual([n['id'] for n in chain['nodes']], ['S1'])
        self.assertEqual(chain['nodes'][0]['internal_covalent_edges'], 1)
        edge, = chain['edges']
        self.assertEqual((edge['source'], edge['target'], edge['type']), ('S1', 'S1', 'backbone_join'))
        self.assertIn(edge['symmetry_code'], ('1_655', '1_455'))
        part, = chain['periodicity']
        self.assertEqual((part['dimension'], part['translation_basis']), (1, ['[1, 0, 0]']))
        residues = doc['residue_graph']
        self.assertEqual([(e['source'], e['target'], e['type'], e['symmetry_code']) for e in residues['edges']],
                         [('S1:DG1', 'S1:DG2', 'covalent', '1_555'), ('S1:DG1', 'S1:DG2', 'backbone_join', '1_455')])
        self.assertTrue(residues['edges'][0]['intra_chain'])
        self.assertEqual(residues['summary'], ['A (2) residues: 1-periodic rod/chain along [1, 0, 0]'])
        self.assertEqual(r['asu_graphs']['chain']['periodicity'], chain['summary'])
        self.assertEqual([a['name'] for a in r['artifacts']], ['g.cgd', 'asu.json', 'chain.html', 'residue.html'])

    def test_unstitched_graph_is_covalent_only(self):
        r, doc, _ = self.export(self.nick_chain())
        self.assertEqual(r['schema_version'], 1)
        self.assertEqual(doc['stitching'], {'nicks': False, 'terminal_stacks': False})
        self.assertEqual(doc['chain_graph']['edges'], [])
        self.assertEqual(doc['chain_graph']['summary'], ['1 single node(s) with no join, finite: A'])
        self.assertEqual(doc['residue_graph']['edge_type_counts'], {'covalent': 1})
        self.assertEqual(doc['residue_graph']['periodicity'][0]['dimension'], 0)

    def test_gap_splits_residue_pieces_but_not_the_chain(self):
        text = crystal(60, 60, 60) + strand(1, 'A', (10, 13.4, 16.8), drop_p=(3,))
        _, doc, _ = self.export(text)
        self.assertEqual([len(p['nodes']) for p in doc['residue_graph']['periodicity']], [2, 1])
        self.assertEqual(doc['chain_graph']['nodes'][0]['internal_covalent_edges'], 1)

    def test_terminal_stacks_record_polarity(self):
        r, doc, out = self.export(crystal(60, 60, 60) + strand(1, 'A', (10, 13.4))
                                  + strand(100, 'B', (20.2, 16.8), dx=-5), stitch_terminal_stacks=True)
        stack, = [e for e in doc['chain_graph']['edges'] if e['type'] == 'terminal_stack']
        self.assertEqual((stack['source'], stack['target'], stack['symmetry_code']), ('S1', 'S2', '1_555'))
        self.assertEqual((stack['polarities'], stack['same_polarity']), (["3'-3'"], True))
        self.assertEqual(stack['residue_pairs'], ["A:DC2 <-> B:DC2 (3'-3')"])
        residue_stack, = [e for e in doc['residue_graph']['edges'] if e['type'] == 'terminal_stack']
        self.assertIn("polarity 3'-3'", residue_stack['contacts'][0])
        page = (out/'chain.html').read_text()
        self.assertIn('stroke-dasharray', page)
        self.assertIn('Stitched chain-end pi stack', page)

    def test_deposited_links_typed_and_untyped(self):
        dna = strand(1, 'A', (10, 13.4))
        o3 = (15., 10., 14.2)  # A2 O3'
        lig = (f"HETATM  200  C1  LIG C   1    {o3[0]+1.4:8.3f}{o3[1]:8.3f}{o3[2]:8.3f}  1.00 20.00           C\n"
               f"HETATM  201 MG    MG D   1    {o3[0]:8.3f}{o3[1]+2.1:8.3f}{o3[2]:8.3f}  1.00 20.00          MG\n"
               f"HETATM  202  O   HOH E   1    {o3[0]:8.3f}{o3[1]-1.0:8.3f}{o3[2]:8.3f}  1.00 20.00           O\n")
        serial = [line[6:11].strip() for line in dna.splitlines() if line[12:16].strip() == "O3'"][-1]
        conect = (f'CONECT{int(serial):5d}  200\nCONECT{int(serial):5d}  201\nCONECT{int(serial):5d}  202\n')
        _, doc, _ = self.export(crystal(60, 60, 60) + dna + lig + conect)
        self.assertEqual([n['label'] for n in doc['chain_graph']['nodes']], ['A', 'C'])
        edge, = doc['chain_graph']['edges']
        self.assertEqual((edge['type'], edge['residue_pairs']), ('covalent', ['A:DC2 <-> C:LIG1']))
        self.assertEqual(doc['chain_graph']['nodes'][1]['kind'], 'other')

    def test_hydrogen_bond_annotation_is_not_covalent(self):
        # 5KEK/9ZXK carry _struct_conn hydrog records at 2.27 A; only covalent types count.
        from test_mmcif_input import HEADER, TAGS, site
        rows = (site(1, "O3'", label='A', chain='A', xyz=(1., 1., 1.)) + site(2, "C3'", label='A', chain='A', xyz=(1., 1., 2.4))
                + site(3, "O3'", label='B', chain='B', xyz=(3.27, 1., 1.)) + site(4, "C3'", label='B', chain='B', xyz=(3.27, 1., 2.4)))
        conn = ('loop_\n_struct_conn.id\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n'
                '_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n'
                '_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\n')
        for kind, expected in (('hydrog', []), ('covale', ['covalent'])):
            with self.subTest(kind=kind):
                self.source = self.root / f'{kind}.cif'
                self.source.write_text(HEADER + TAGS + rows + conn + f'c1 {kind} A 1 "O3\'" B 1 "O3\'"\n')
                _, doc, _ = self.export(None)
                self.assertEqual([e['type'] for e in doc['chain_graph']['edges']], expected)

    def test_html_is_self_contained(self):
        _, doc, out = self.export(self.nick_chain(), stitch_nicks=True)
        for name in ('chain.html', 'residue.html'):
            page = (out/name).read_text()
            self.assertIn('<svg', page)
            self.assertNotIn('<script src', page)
            self.assertNotIn('<link', page)
            self.assertNotIn('http://', page.replace('http://www.w3.org/2000/svg', ''))
            self.assertNotIn('https://', page)
        self.assertIn('1_455', (out/'residue.html').read_text())

    def test_outputs_are_never_overwritten(self):
        self.source.write_text(self.nick_chain())
        taken = self.root/'taken.html'
        taken.write_text('keep')
        with self.assertRaises(ConversionError):
            export_cgd(self.source, self.root/'g.cgd', asu_chain_html=taken)
        self.assertEqual(taken.read_text(), 'keep')
        self.assertFalse((self.root/'g.cgd').exists())
        with self.assertRaises(ConversionError):
            export_cgd(self.source, self.root/'g.cgd', asu_graphs=self.root/'same', asu_chain_html=self.root/'same')


if __name__ == '__main__':
    unittest.main()
