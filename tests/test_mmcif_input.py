import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from topotrace.converter import ConversionError, convert_file, read_pdb
from test_converter import atom, crystal, loop_rows


HEADER = '''data_example
_entry.id example
_cell.length_a 20
_cell.length_b 30
_cell.length_c 40
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
_symmetry.space_group_name_H-M 'P 1'
loop_
_space_group_symop.id
_space_group_symop.operation_xyz
1 x,y,z
'''
TAGS = '''loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
'''

def site(serial, name, *, label='X', chain='A', seq='1', alt='.', model='1', xyz=(1.1234567, 2., 3.), b=20.):
    return f'ATOM {serial} {name[0]} "{name}" {alt} DG {label} {seq} {chain} {seq} ? {xyz[0]} {xyz[1]} {xyz[2]} 1 {b} {model}\n'


class MmcifInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.input = self.root/'input.cif'
        self.output = self.root/'output.cif'

    def convert(self, rows, **kwargs):
        self.input.write_text(HEADER + TAGS + rows)
        return convert_file(self.input, self.output, **kwargs)

    def test_native_precision_and_content_detection(self):
        self.input = self.root/'misleading.pdb'
        r = self.convert(site(1,"C4'")+site(2,"O5'")+site(3,"O3'"), atoms=["O3'","C4'","O5'"])
        self.assertEqual((r.atom_count,r.bond_count),(3,2))
        row=loop_rows(self.output.read_text(),'_atom_site_')[0]
        self.assertAlmostEqual(float(row['_atom_site_fract_x'])*20,1.1234567,places=8)

    def test_label_asym_boundaries_preserve_repeated_author_chain(self):
        r=self.convert(site(1,"O5'",label='X')+site(2,"O3'",label='X')
                       +site(3,"O5'",label='Y')+site(4,"O3'",label='Y'),atoms=["O5'","O3'"])
        self.assertEqual((r.atom_count,r.bond_count,r.segment_count),(4,2,2))

    def test_model_ordinals_and_altloc_selection(self):
        self.input.write_text(HEADER+TAGS+site(1,'P',model='5')+site(1,'P',model='9',alt='A')+site(2,'P',model='9',alt='B'))
        with self.assertRaisesRegex(ConversionError,'altloc'):
            convert_file(self.input,self.output,model=2,atoms=['P'])
        r=convert_file(self.input,self.output,model=2,atoms=['P'],altloc='B')
        self.assertEqual((r.model,r.atom_count,r.bond_count),(2,1,0))
        self.assertTrue(r.warnings)

    def test_label_identity_fallback_and_insertion_codes(self):
        rows=site(1,"O5'",chain='?',seq='2')+site(2,"O3'",chain='?',seq='2')
        r=self.convert(rows.replace('2 ? 1.1234567','2 A 1.1234567'),atoms=["O5'","O3'"])
        self.assertEqual(r.bond_count,1)
        self.assertIn('chain=X',self.output.read_text())
        self.assertIn('DG:2A',self.output.read_text())

    def test_mmcif_connections_roundtrip_without_bond_perception(self):
        original=self.root/'original.pdb'
        original.write_text(crystal()+atom(1,"O5'")+atom(2,"O3'")+ 'CONECT    1    2\n')
        convert_file(original,self.input,output_format='mmcif')
        r=convert_file(self.input,self.output)
        self.assertEqual(r.bond_count,1)

    def test_mmcif_to_pdb_is_checked_and_preserves_coordinates(self):
        self.output=self.root/'out.pdb'
        r=self.convert(site(1,"O5'",xyz=(1.25,2,3))+site(2,"O3'",xyz=(4,5,6)),atoms=["O5'","O3'"],output_format='pdb')
        readback=read_pdb(self.output.read_text())
        self.assertEqual(readback.cell,(20.,30.,40.,90.,90.,90.))
        self.assertEqual(readback.atoms[0].xyz,(1.25,2.,3.))
        self.assertEqual(r.bond_count,1)

    def test_pdb_field_overflow_fails_before_output(self):
        for rows in (site(1,'P',chain='LONG'),site(1,'P',b=1200.1234)):
            self.input.write_text(HEADER+TAGS+rows)
            with self.assertRaisesRegex(ConversionError,'PDB'):
                convert_file(self.input,self.output,output_format='pdb')
            self.assertFalse(self.output.exists())

    def test_invalid_coordinate_partial_cell_and_core_input_are_rejected(self):
        for text in (HEADER+TAGS+site(1,'P',xyz=('nan',0,0)),
                     HEADER.replace('_cell.length_b 30\n','')+TAGS+site(1,'P'),
                     "data_core\n_atom_site_fract_x 0.25\n"):
            self.input.write_text(text)
            with self.assertRaises(ConversionError): convert_file(self.input,self.output)
            self.assertFalse(self.output.exists())

    def test_explicit_input_format_and_cli_base_rejection(self):
        self.input.write_text(HEADER+TAGS+site(1,"C4'")+site(2,'N9'))
        r=subprocess.run([sys.executable,'-m','topotrace',str(self.input),'--input-format','mmcif',
            '-atom',"C4'",'-out',str(self.output),'--json'],capture_output=True,text=True)
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(json.loads(r.stdout)['atom_count'],1)
        rejected=subprocess.run([sys.executable,'-m','topotrace',str(self.input),'-atom','N9',
            '-out',str(self.root/'reject.cif')],capture_output=True,text=True)
        self.assertEqual(rejected.returncode,2)
        self.assertIn('nucleobase',rejected.stderr.lower())
        self.assertFalse((self.root/'reject.cif').exists())

    def test_extra_sugar_is_recognized_but_not_given_a_fake_linear_order(self):
        self.input.write_text(crystal()+atom(1,"C1'"))
        with self.assertRaisesRegex(ConversionError,'sugar.*linear'):
            convert_file(self.input,self.output,atoms=["C1'"])

    def test_native_small_nonstandard_transform_is_not_treated_as_pdb_rounding(self):
        transform=''.join(f'_atom_sites.fract_transf_matrix[{i}][{j}] {v}\n'
            for i,row in enumerate(((.05,.0000001,0),(0,1/30,0),(0,0,.025)),1)
            for j,v in enumerate(row,1))
        transform+=''.join(f'_atom_sites.fract_transf_vector[{i}] {v}\n' for i,v in enumerate((.000001,0,0),1))
        self.input.write_text(HEADER+transform+TAGS+site(1,'P',xyz=(10,10,10)))
        convert_file(self.input,self.output)
        row=loop_rows(self.output.read_text(),'_atom_site_')[0]
        self.assertAlmostEqual(float(row['_atom_site_fract_x']),.500002,places=10)

    def test_symmetry_and_unknown_image_connections_are_reported_not_fabricated(self):
        conn='''loop_
_struct_conn.id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_symmetry
one X 1 "O5'" X 1 "O3'" 1_555 1_655
two X 1 "O5'" X 1 "O3'" 1_555 ?
'''
        self.input.write_text(HEADER+TAGS+site(1,"O5'")+site(2,"O3'")+conn)
        r=convert_file(self.input,self.output)
        self.assertEqual(r.bond_count,0)
        self.assertEqual(len(r.warnings),2)

    def test_modified_base_name_rejection_uses_actual_residue_mapping(self):
        self.input.write_text((HEADER+TAGS+site(1,'N24')).replace(' DG ', ' 4E9 '))
        with self.assertRaisesRegex(ConversionError,'nucleobase.*N24'):
            convert_file(self.input,self.output,atoms=['N24'])
        # Full conversion still retains bases and unknown atoms.
        r=convert_file(self.input,self.output)
        self.assertEqual(r.atom_count,1)

    def test_cif_tag_names_are_case_insensitive(self):
        text=HEADER+TAGS+site(1,"O5'")+site(2,"O3'")
        text='\n'.join((line.split(' ',1)[0].upper()+(' '+line.split(' ',1)[1] if ' ' in line else ''))
                       if line.startswith('_') else line for line in text.splitlines())+'\n'
        self.input.write_text(text)
        r=convert_file(self.input,self.output,atoms=["O3'","O5'"])
        self.assertEqual((r.atom_count,r.bond_count),(2,1))


if __name__ == '__main__':
    unittest.main()
