# Nucleotide atom recognition

The 2026-11-09 request adds native PDBx/mmCIF input and an internal recognition
layer for future sugar ribbons. The user explicitly retained the existing
linear backbone CLI; sugar-ring connectivity and ribbons are future work.

## Lookup semantics

NARestraints `residue_library.py`, `builder.record_for_residue`, and
`builder.pair_residue_from_pdb` establish the flow: use the structure's residue
code, resolve canonical-role columns to the actual names in that residue's row,
and include only names present in the coordinates. `dna_geometry.py` likewise
distinguishes a mapped heterocycle from a possible sugar attachment. NASolve
`docs/architecture.md` separates atom-role mappings from component bond and
stereochemical dictionaries. TopoStich follows that distinction and does not
construct chemical bonds or sugar stereochemistry from this spreadsheet.

The packaged `src/topostich/data/nucleotide_atoms.json` is a snapshot of
`NARestraints/restraints/data/Ligands.xlsx`, SHA-256
`57a566c710f45df06dc8b5174f135dfea6f75db781755905a8570523fb57d6e5`.
It was read from the user-local working workbook, not asserted to equal any
published NARestraints release. The original workbook and its lockfile are
untouched. Only atom mappings, required residue metadata and source-row
provenance are extracted; the workbook itself is not distributed.

There are 600 residue codes and 619 rows from the 13 named chemistry sheets.
The `sort` workspace sheet is excluded, matching NARestraints' authoritative
sheet list. Column headers such as `N9` are canonical positions. Cell values
are actual atom names: DZ's canonical N1 maps to C1; 4E9's N2 maps to N24.
That does not make C1 or N24 a base atom in every other residue.

Header roles, substituent mappings and annotation cells stay distinct.
`C5-mod` and `C7-mod` cells are included only when their values are atom-like
identifiers, not descriptive phrases. Excluded cells are recorded by location.
Names are trimmed; star and Unicode prime aliases normalize to ASCII prime.
Recognition also maps O1P/O2P/O3P to OP1/OP2/OP3. Conversion preserves the source
atom identity and does not rename a structure from the workbook.

## Internal API

- `lookup_residue(code)` returns every workbook definition for that code.
- `classify_atom(code, name, observed_atoms=...)` reports a normalized name,
  role, mapped canonical positions, candidate roles and mapping status.
- `recognize_residue(code, observed_names)` groups observed sugar, phosphate,
  base, unknown and ambiguous atoms and reports missing mapped base names.
  Call it for each model/segment/chain/residue/insertion identity separately.

Recognition covers conventional sugar and phosphate heavy atoms and the base
atoms mapped in the workbook. Hydrogens, arbitrary new sugar analog atom names
and unmapped modifications can remain unknown. This is not complete CCD
perception. Known library membership or conventional sugar names provide
`nucleotide_like` evidence; this is not a claim of validated chemical topology.

For duplicate definitions, score distinct mapped names actually observed in the
residue, retaining all rows tied for the highest score. Equivalent mappings are
identified; conflicting ties are `ambiguous`. No first-row tie-break or union
of incompatible chemistry is presented as one resolved map. Only unanimously
mapped absent names are reported as missing when a tie remains unresolved.

A sugar/phosphate name that occurs as a mapped base value remains ambiguous
for that residue (e.g. the workbook's IOO OP1 or PSU C1′ entries). Such a cell
does not poison a global name blacklist. Component-specific NASolve adaptations
such as DF O2→S1 are not copied into every DF record: the snapshot retains its
source mapping, and unmapped observed names remain visible as unknown.

CLI rejection is a small consumer of this API. Canonical base headers and
observed residue-specific base names produce a nucleobase diagnostic. Extra
recognized sugar atoms receive an explicit outside-linear-path diagnostic.
No-filter conversion still includes all atoms, including bases and unknowns.
The receipt remains schema 1; recognition data are available through the Python
API without changing existing machine consumers.

## Refreshing the snapshot

Install the development extra `.[catalog]`, then run:

```sh
python tools/build_atom_catalog.py /explicit/path/to/Ligands.xlsx /new/path/nucleotide_atoms.json
```

The tool is read-only on the workbook and refuses an existing destination.
Review the resulting mapping diff, source hash, excluded cells and regression
tests before replacing the packaged snapshot. Run the complete build/install
gate afterward. No runtime absolute path or workbook installation is required.
