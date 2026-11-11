import hashlib
import json
import math
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile
import unittest

from topostich.converter import ConversionError, convert_file


def atom(serial, name, residue=1, chain="A", xyz=(1., 2., 3.), *,
         alt="", insertion="", element=None, record="ATOM"):
    element = element or name[0]
    return (f"{record:<6}{serial:5d} {name:>4}{alt:1} DG {chain:1}"
            f"{residue:4d}{insertion:1}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"{1.:6.2f}{20.:6.2f}          {element:>2}\n")


def crystal(a=20., b=30., c=40., alpha=90., beta=90., gamma=90., group="P 1"):
    return (f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{alpha:7.2f}{beta:7.2f}"
            f"{gamma:7.2f} {group:<11}   1\n")


def loop_rows(text, prefix):
    """Read our simple scalar CIF loops independently of the writer."""
    lines = iter(text.splitlines())
    for line in lines:
        if line != "loop_":
            continue
        headers = []
        for line in lines:
            if line.startswith("_"):
                headers.append(line)
            else:
                break
        values = []
        while line.strip() and not line.startswith("#"):
            values.extend(shlex.split(line))
            line = next(lines, "")
        if headers and headers[0].startswith(prefix):
            assert len(values) % len(headers) == 0
            return [dict(zip(headers, values[n:n+len(headers)]))
                    for n in range(0, len(values), len(headers))]
    return []


class ConversionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "input.pdb"
        self.output = self.root / "output.cif"

    def convert(self, text, **kwargs):
        self.source.write_text(text)
        result = convert_file(self.source, self.output, **kwargs)
        self.assertEqual(self.source.read_text(), text)
        return result, self.output.read_text()

    def test_core_cif_preserves_all_atoms_without_filter(self):
        result, text = self.convert(crystal() + atom(1, "O5'") + atom(2, "N9"))
        rows = loop_rows(text, "_atom_site_")
        self.assertEqual(result.atom_count, 2)
        self.assertEqual(result.bond_count, 0)
        self.assertEqual(rows[0]["_atom_site_type_symbol"], "O")
        self.assertAlmostEqual(float(rows[0]["_atom_site_fract_x"]), 0.05)
        self.assertAlmostEqual(float(rows[0]["_atom_site_fract_y"]), 2/30)
        self.assertAlmostEqual(float(rows[0]["_atom_site_U_iso_or_equiv"]), 20/(8*math.pi**2))
        self.assertNotIn("_atom_site.Cartn_x", text)

    def test_backbone_order_missing_atoms_and_no_stitching(self):
        records = [atom(1, "O3'", xyz=(9, 0, 0)), atom(2, "C4'", xyz=(6, 0, 0)),
                   atom(3, "O5'", xyz=(3, 0, 0)), atom(4, "O3'", 2, xyz=(15, 0, 0)),
                   atom(5, "O5'", 2, xyz=(12, 0, 0)), "TER\n",
                   atom(6, "O5'", 3), atom(7, "O3'", 3),
                   atom(8, "O5'", 1, "B"), atom(9, "O3'", 1, "B")]
        result, text = self.convert(crystal() + "".join(records), atoms=["O3'", "O5'", "C4'"])
        self.assertEqual((result.atom_count, result.bond_count, result.segment_count), (9, 6, 3))
        self.assertEqual(result.missing_count, 3)
        rows = loop_rows(text, "_atom_site_")
        self.assertEqual([float(row["_atom_site_fract_x"]) for row in rows[:5]], [.15, .3, .45, .6, .75])
        bonds = loop_rows(text, "_geom_bond_")
        self.assertEqual(len(bonds), 6)
        self.assertEqual(float(bonds[0]["_geom_bond_distance"]), 3.)
        self.assertEqual(len(loop_rows(text, "_topol_link.")), 6)

    def test_entire_missing_residue_is_bridged_and_numbering_not_sorted(self):
        result, _ = self.convert(crystal() + atom(1, "O3'", 10)
                                 + atom(2, "N9", 11) + atom(3, "O5'", 9),
                                 atoms=["O5'", "O3'"])
        self.assertEqual((result.atom_count, result.bond_count, result.missing_count), (2, 1, 4))

    def test_insertion_codes_and_legacy_star_names(self):
        result, text = self.convert(crystal() + atom(1, "O5*", 1) + atom(2, "O3*", 1)
                                   + atom(3, "O5'", 1, insertion="A"), atoms=["O3'", "O5*"])
        self.assertEqual((result.atom_count, result.bond_count), (3, 2))
        self.assertIn("1A", text)

    def test_nonorthogonal_cell_and_fractional_symmetry(self):
        header = crystal(105.547, 105.547, 96.268, gamma=120., group="H 3")
        header += "REMARK 290       1555   X,Y,Z\nREMARK 290       2555   -Y,X-Y,Z\n"
        _, text = self.convert(header + atom(1, "O5'", xyz=(29.654, 6.408, 3.337)))
        row = loop_rows(text, "_atom_site_")[0]
        fy = 6.408 / (105.547 * math.sin(math.radians(120)))
        fx = 29.654 / 105.547 - fy * math.cos(math.radians(120))
        self.assertAlmostEqual(float(row["_atom_site_fract_x"]), fx)
        self.assertAlmostEqual(float(row["_atom_site_fract_y"]), fy)
        self.assertIn("'R 3'", text)
        self.assertIn("-y,x-y,z", text)

    def test_nonstandard_scale_and_translation_are_used(self):
        scale = "".join(f"SCALE{i+1}    {row[0]:10.6f}{row[1]:10.6f}{row[2]:10.6f}     {row[3]:10.5f}\n"
                        for i, row in enumerate(((0, .05, 0, .25), (-1/30, 0, 0, 0), (0, 0, .025, 0))))
        _, text = self.convert(crystal() + scale + atom(1, "P", xyz=(3, 4, 8)))
        row = loop_rows(text, "_atom_site_")[0]
        self.assertAlmostEqual(float(row["_atom_site_fract_x"]), .45)
        self.assertAlmostEqual(float(row["_atom_site_fract_y"]), -.099999)
        self.assertAlmostEqual(float(row["_atom_site_fract_z"]), .2)

    def test_missing_cell_rejected_before_output(self):
        self.source.write_text(atom(1, "P"))
        with self.assertRaisesRegex(ConversionError, "CRYST1"):
            convert_file(self.source, self.output)
        self.assertFalse(self.output.exists())

    def test_no_selection_is_an_error_not_empty_success(self):
        self.source.write_text(crystal() + atom(1, "N9"))
        with self.assertRaisesRegex(ConversionError, "No atoms"):
            convert_file(self.source, self.output, atoms=["O3'"])
        self.assertFalse(self.output.exists())

    def test_models_never_join(self):
        result, _ = self.convert(crystal() + "MODEL        1\n" + atom(1, "P")
                                 + "ENDMDL\nMODEL        2\n" + atom(1, "P", xyz=(9, 8, 7))
                                 + "ENDMDL\n", atoms=["P"], model=2)
        self.assertEqual((result.atom_count, result.bond_count), (1, 0))

    def test_altloc_requires_choice_for_trace_but_plain_conversion_retains_both(self):
        pdb = crystal() + atom(1, "P", alt="A") + atom(2, "P", alt="B")
        self.source.write_text(pdb)
        with self.assertRaisesRegex(ConversionError, "altloc"):
            convert_file(self.source, self.output, atoms=["P"])
        result = convert_file(self.source, self.output, atoms=["P"], altloc="B")
        self.assertEqual(result.atom_count, 1)
        result = convert_file(self.source, self.root / "all.cif")
        self.assertEqual(result.atom_count, 2)

    def test_existing_output_and_input_never_overwritten(self):
        self.source.write_text(crystal() + atom(1, "P"))
        self.output.write_text("keep me")
        with self.assertRaises(FileExistsError):
            convert_file(self.source, self.output)
        self.assertEqual(self.output.read_text(), "keep me")
        with self.assertRaises(FileExistsError):
            convert_file(self.source, self.source)

    def test_mmcif_and_pdb_export_connections(self):
        self.source.write_text(crystal() + atom(12, "O5'") + atom(42, "O3'", xyz=(5, 6, 7)))
        for fmt in ("pdb", "mmcif"):
            target = self.root / ("out." + fmt)
            result = convert_file(self.source, target, atoms=["O3'", "O5'"], output_format=fmt)
            self.assertEqual(result.bond_count, 1)
            text = target.read_text()
            if fmt == "pdb":
                self.assertIn("CONECT    1    2", text)
                self.assertEqual(text.count("\nTER"), 1)
            else:
                self.assertEqual(len(loop_rows(text, "_atom_site.")), 2)
                self.assertEqual(len(loop_rows(text, "_struct_conn.")), 1)
                self.assertIn('"O5\'"', text)

    def test_conect_preserved_only_when_present_in_unfiltered_input(self):
        result, text = self.convert(crystal() + atom(10, "C1'") + atom(20, "N9")
                                   + "CONECT   10   20\nCONECT   20   10\n")
        self.assertEqual(result.bond_count, 1)
        self.assertEqual(len(loop_rows(text, "_geom_bond_")), 1)

    def test_duplicate_identity_rejected_for_trace(self):
        self.source.write_text(crystal() + atom(1, "P") + atom(2, "P"))
        with self.assertRaisesRegex(ConversionError, "Duplicate"):
            convert_file(self.source, self.output, atoms=["P"])

    def test_invalid_coordinate_and_partial_scale_fail_cleanly(self):
        for pdb in (crystal() + atom(1, "P", xyz=(float("nan"), 0, 0)),
                    crystal() + "SCALE1        0.050000  0.000000  0.000000        0.00000\n" + atom(1, "P")):
            self.source.write_text(pdb)
            with self.assertRaises(ConversionError):
                convert_file(self.source, self.output)
            self.assertFalse(self.output.exists())

    def test_cli_quoted_primes_format_aliases_and_json_receipt(self):
        self.source.write_text(crystal() + atom(1, "O5'") + atom(2, "O3'"))
        command = [sys.executable, "-m", "topostich", str(self.source), "-atom", "O3'",
                   "-atom", "O5'", "-out", str(self.output), "-cif", "--json"]
        completed = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        receipt = json.loads(completed.stdout)
        self.assertEqual(receipt["atom_count"], 2)
        self.assertEqual(receipt["schema_version"], 1)
        self.assertEqual(receipt["source_sha256"], hashlib.sha256(self.source.read_bytes()).hexdigest())
        self.assertEqual(receipt["output_sha256"], hashlib.sha256(self.output.read_bytes()).hexdigest())


if __name__ == "__main__":
    unittest.main()
