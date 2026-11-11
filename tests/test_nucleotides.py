import unittest

from topostich.nucleotides import classify_atom, recognize_residue, lookup_residue


class NucleotideRecognitionTests(unittest.TestCase):
    def test_standard_sugar_and_phosphate_aliases(self):
        for name in ("C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'", "O4'", "O5'", 'C4*', 'O3′'):
            self.assertEqual(classify_atom('DG', name).role, 'sugar', name)
        for name in ('P', 'OP1', 'OP2', 'OP3', 'O1P', 'O2P', 'O3P'):
            self.assertEqual(classify_atom('DG', name).role, 'phosphate', name)
        self.assertEqual(classify_atom('DG', 'O1P').normalized_name, 'OP1')

    def test_actual_modified_names_are_looked_up_by_residue(self):
        self.assertEqual(classify_atom('DZ', 'C1').canonical_positions, ('N1',))
        self.assertEqual(classify_atom('DZ', 'C1').role, 'base')
        self.assertEqual(classify_atom('DG', 'C1').role, 'unknown')
        self.assertEqual(classify_atom('S6G', 'C5 ').canonical_positions, ('C5',))
        self.assertEqual(classify_atom('4E9', 'N24').canonical_positions, ('N2',))
        self.assertEqual(classify_atom('DG', 'N24').role, 'unknown')

    def test_presence_is_reported_without_inventing_missing_atoms(self):
        result = recognize_residue('DZ', ["C1'", "C4'", 'C1', 'N3', 'UNLISTED'])
        self.assertTrue(result.nucleotide_like)
        self.assertEqual(result.present_base_atoms, ('C1', 'N3'))
        self.assertEqual(result.present_sugar_atoms, ("C1'", "C4'"))
        self.assertEqual(result.unknown_atoms, ('UNLISTED',))
        self.assertIn('O2', result.missing_mapped_base_atoms)

    def test_unknown_residue_is_not_assigned_a_canonical_base(self):
        result = recognize_residue('ZZZ', ["C1'", "C4'", 'N9'])
        self.assertTrue(result.nucleotide_like)
        self.assertEqual(result.mapping_status, 'unmapped')
        self.assertEqual(classify_atom('ZZZ', 'N9').role, 'unknown')
        self.assertFalse(recognize_residue('HOH', ['O']).nucleotide_like)
        self.assertFalse(recognize_residue('LIG', ['N1']).nucleotide_like)

    def test_duplicate_rows_and_role_collisions_remain_explicit(self):
        self.assertGreater(len(lookup_residue('A40')), 1)
        self.assertEqual(recognize_residue('A40', ['N9']).mapping_status, 'ambiguous')
        matched = recognize_residue('A40', ['N9', 'N2'])
        self.assertEqual(matched.mapping_status, 'matched')
        self.assertEqual(matched.definitions[0].source_sheet, 'D')
        self.assertEqual(classify_atom('IOO', 'OP1').role, 'ambiguous')
        self.assertEqual(classify_atom('DG', 'OP1').role, 'phosphate')
        self.assertEqual(classify_atom('PSU', "C1'").role, 'ambiguous')
        self.assertEqual(classify_atom('DG', "C1'").role, 'sugar')

    def test_metadata_is_not_an_atom_and_workbook_stays_a_lookup(self):
        self.assertEqual(classify_atom('C', "CYTIDINE-5'-MONOPHOSPHATE").role, 'unknown')
        self.assertEqual(classify_atom('C', 'C5-mod').role, 'unknown')
        self.assertEqual(classify_atom('DF', 'S2').role, 'base')
        # NASolve's component-specific S1 correction is not silently globalized.
        self.assertEqual(classify_atom('DF', 'S1').role, 'unknown')


if __name__ == '__main__':
    unittest.main()
