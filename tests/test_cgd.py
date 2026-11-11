"""Explicit graph and viewing contracts; no scientific inputs are distributed."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

from test_converter import atom, crystal
from topostich.converter import ConversionError, convert_file
from topostich.cgd import export_cgd


class CgdTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'input.pdb'
        self.out = self.root / 'graph.cgd'

    def export(self, text, **kw):
        self.source.write_text(text)
        return export_cgd(self.source, self.out, **kw)

    def test_chain_ter_missing_and_exact_dialect(self):
        r = self.export(crystal() + atom(1,"O3'",1,xyz=(1,2,3))
                        + atom(2,"C4'",2,xyz=(2,3,4))
                        + atom(3,"O3'",3,xyz=(5,6,7))
                        + atom(4,"O3'",4,xyz=(8,9,10)) + 'TER\n'
                        + atom(5,"O3'",5,xyz=(9,10,11))
                        + atom(6,"O3'",1,'B',xyz=(10,11,12)))
        self.assertEqual((r['node_count'],r['edge_count'],r['segment_count']),(5,2,3))
        self.assertEqual([n['degree'] for n in r['nodes']],[1,2,1,0,0])
        self.assertEqual(len(r['missing_atoms']),1)
        self.assertEqual(r['edges'],[['S1_N1','S1_N2'],['S1_N2','S1_N3']])
        lines = self.out.read_text().splitlines()
        self.assertEqual(lines[0],'crystal')
        self.assertEqual(lines[-1],'end')
        self.assertEqual([x.split()[0] for x in lines[:5]],['crystal','name','cell','group','atom'])
        edges = [x.split() for x in lines if x.startswith('edge ')]
        self.assertEqual(len(edges),4)
        self.assertTrue(all(len(x)==5 for x in edges))
        self.assertEqual(edges[0][1], 'S1_N1')
        self.assertEqual([float(x) for x in edges[0][2:]],[.25,.2,.175])

    def test_unwrapped_coordinates_and_cmm_identity(self):
        cmm = self.root/'graph.cmm'
        r = self.export(crystal()+atom(1,"O3'",xyz=(-1,2,3))
                        +atom(2,"O3'",2,xyz=(21,4,5)),cmm=cmm)
        self.assertAlmostEqual(r['nodes'][0]['fractional'][0],-.05)
        self.assertAlmostEqual(r['nodes'][1]['fractional'][0],1.05)
        xml=ET.parse(cmm).getroot()
        markers=xml.findall('.//marker'); links=xml.findall('.//link')
        self.assertEqual((len(markers),len(links)),(2,1))
        self.assertEqual(float(markers[0].get('x')),-1)
        self.assertIn('chain=A',markers[0].get('note'))
        self.assertEqual((links[0].get('id1'),links[0].get('id2')),('1','2'))

    def test_mmcif_uses_same_graph_without_pdb_intermediate(self):
        self.source.write_text(crystal()+atom(1,"O3'",xyz=(1,2,3))+atom(2,"O3'",2,xyz=(4,5,6)))
        native=self.root/'native.cif'
        convert_file(self.source,native,output_format='mmcif')
        text=native.read_text().replace('1 2 3','1.123456789 2 3')
        native.write_text(text)
        r=export_cgd(native,self.out)
        self.assertEqual(r['source_format'],'mmcif')
        self.assertAlmostEqual(r['nodes'][0]['cartesian'][0],1.123456789,places=9)
        self.assertEqual(r['edge_count'],1)

    def test_r3_cell_and_symmetry_copies(self):
        cmm=self.root/'crystal.cmm'
        r=self.export(crystal(100,100,90,90,90,120,'H 3')
                      +atom(1,"O3'",xyz=(10,17,28))+atom(2,"O3'",2,xyz=(15,20,30)),crystal_cmm=cmm)
        self.assertEqual(r['space_group'],'R3:H')
        self.assertIn('group R3:H\n',self.out.read_text())
        xml=ET.parse(cmm).getroot()
        self.assertEqual(len(xml.findall('.//marker')),18)
        self.assertEqual(len(xml.findall('.//link')),9)
        self.assertEqual(len(r['symmetry_operations']),9)

    def test_reject_symmetry_coincidences_and_special_positions(self):
        for text in [crystal()+atom(1,"O3'",xyz=(1,2,3))+atom(2,"O3'",2,'B',xyz=(21,2,3)),
                     crystal(100,100,90,90,90,120,'H 3')+atom(1,"O3'",xyz=(0,0,20))]:
            with self.subTest(text=text), self.assertRaisesRegex(ConversionError,'symmetry|special'):
                self.export(text)
            self.assertFalse(self.out.exists())

    def test_altloc_model_and_input_unchanged(self):
        text=crystal()+'MODEL        1\n'+atom(1,"O3'",xyz=(1,2,3))+'ENDMDL\nMODEL        2\n'+atom(2,"O3'",xyz=(4,5,6),alt='A')+'ENDMDL\n'
        self.source.write_text(text)
        with self.assertRaises(ConversionError):
            export_cgd(self.source,self.out,model=2)
        r=export_cgd(self.source,self.out,model=2,altloc='A')
        self.assertEqual(r['nodes'][0]['cartesian'],[4.,5.,6.])
        self.assertEqual(self.source.read_text(),text)

    def test_no_partial_outputs_or_overwrites(self):
        text=crystal()+atom(1,"O3'",xyz=(1,2,3))
        occupied=self.root/'occupied.cmm'; occupied.write_text('keep')
        for kw in ({'cmm':occupied},{'cmm':self.out},{'name':'bad\nend'}):
            with self.subTest(kw=kw),self.assertRaises(ConversionError):
                self.export(text,**kw)
            self.assertFalse(self.out.exists())
        self.assertEqual(occupied.read_text(),'keep')
        with self.assertRaises(ConversionError):
            self.export(atom(1,"O3'"),cmm=self.root/'absent.cmm')
        self.assertFalse((self.root/'absent.cmm').exists())

    def test_installed_compatible_cli_and_receipt(self):
        self.source.write_text(crystal()+atom(1,"O3'",xyz=(1,2,3)))
        p=subprocess.run([sys.executable,'-m','topostich.cgd',str(self.source),'--out',str(self.out),
                          '--cmm',str(self.root/'graph.cmm'),'--receipt',str(self.root/'receipt.json'),'--json'],
                         text=True,capture_output=True)
        self.assertEqual(p.returncode,0,p.stderr)
        r=json.loads(p.stdout)
        self.assertEqual(r,json.loads((self.root/'receipt.json').read_text()))
        self.assertEqual(r['node_count'],1)
        self.assertEqual(r['edge_count'],0)


if __name__=='__main__':
    unittest.main()
