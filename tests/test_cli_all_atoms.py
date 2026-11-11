"""All-atom CLI conversion keeps unfiltered coordinates and explicit links."""

import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from test_converter import atom, crystal, loop_rows
from test_mmcif_input import HEADER, TAGS, site


class AllAtomCliTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "input.pdb"

    def run_cli(self, output, *options):
        return subprocess.run(
            [sys.executable, "-m", "topostich", str(self.source),
             *options, "-out", str(output), "--json"],
            text=True, capture_output=True,
        )

    def test_all_preserves_every_pdb_site_and_does_not_infer_links(self):
        water = atom(5, "O", residue=2, record="HETATM")
        water = water[:17] + "HOH" + water[20:]
        ion = atom(6, "MG", residue=3, element="MG", record="HETATM")
        ion = ion[:17] + " MG" + ion[20:]
        original = (crystal() + atom(1, "C4'", xyz=(-1.25, 32.5, 43.75))
                    + atom(2, "N9", alt="A") + atom(3, "N9", alt="B")
                    + atom(4, "O3'", xyz=(8, 9, 10)) + water + ion)
        self.source.write_text(original)
        explicit = self.root / "all.cif"
        completed = self.run_cli(explicit, "-atom", "-all", "-cif")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        receipt = json.loads(completed.stdout)
        self.assertEqual((receipt["atom_count"], receipt["bond_count"]), (6, 0))
        self.assertEqual(receipt["atom_selection"], [])
        self.assertEqual(receipt["schema_version"], 1)
        rows = loop_rows(explicit.read_text(), "_atom_site_")
        self.assertEqual(len(rows), 6)
        for axis, length, coordinate in zip("xyz", (20, 30, 40), (-1.25, 32.5, 43.75)):
            self.assertAlmostEqual(float(rows[0][f"_atom_site_fract_{axis}"]) * length,
                                   coordinate, places=7)
        implicit = self.root / "implicit.cif"
        completed = self.run_cli(implicit)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(explicit.read_bytes(), implicit.read_bytes())
        self.assertEqual(self.source.read_text(), original)

    def test_all_keeps_native_mmcif_precision_and_declared_connections(self):
        self.source = self.root / "input.mmcif"
        connection = '''loop_
_struct_conn.id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_symmetry
one X 1 "C1'" X 1 N9 1_555 1_555
'''
        self.source.write_text(HEADER + TAGS + site(1, "C1'") + site(2, "N9") + connection)
        output = self.root / "all.cif"
        completed = self.run_cli(output, "-atom", "-all")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        receipt = json.loads(completed.stdout)
        self.assertEqual((receipt["atom_count"], receipt["bond_count"]), (2, 1))
        row = loop_rows(output.read_text(), "_atom_site_")[0]
        self.assertAlmostEqual(float(row["_atom_site_fract_x"]) * 20, 1.1234567, places=8)

    def test_all_aliases_keep_explicit_model_and_altloc_selection(self):
        self.source.write_text(crystal() + "MODEL        1\n" + atom(1, "P")
                               + "ENDMDL\nMODEL        2\n"
                               + atom(1, "N9", alt="A", xyz=(4, 5, 6))
                               + atom(2, "N9", alt="B", xyz=(7, 8, 9)) + "ENDMDL\n")
        for index, flags in enumerate((("-atom", "-all"), ("--atom", "-all"),
                                       ("-atom", "all"), ("--atom=all",), ("--atom=-all",))):
            with self.subTest(flags=flags):
                output = self.root / f"alias-{index}.cif"
                completed = self.run_cli(output, *flags, "--model", "2", "--altloc", "B")
                self.assertEqual(completed.returncode, 0, completed.stderr)
                receipt = json.loads(completed.stdout)
                self.assertEqual((receipt["model"], receipt["atom_count"], receipt["bond_count"]),
                                 (2, 1, 0))
                row = loop_rows(output.read_text(), "_atom_site_")[0]
                self.assertAlmostEqual(float(row["_atom_site_fract_x"]) * 20, 7)

    def test_all_cannot_be_mixed_with_named_atoms(self):
        self.source.write_text(crystal() + atom(1, "O3'"))
        output = self.root / "rejected.cif"
        for flags in (("-atom", "-all", "-atom", "O3'"),
                      ("-atom", "O3'", "--atom", "all")):
            with self.subTest(flags=flags):
                completed = self.run_cli(output, *flags)
                self.assertEqual(completed.returncode, 2)
                self.assertIn("cannot be combined", completed.stderr)
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
