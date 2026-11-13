# TopoStich

Convert a DNA/RNA **PDB or PDBx/mmCIF** into a **core CIF** for topology software.
Optionally select backbone atoms and connect them in strand order. Python
3.10+; Gemmi is installed automatically. Phenix, Coot, Excel and NASolve are
not required.

Conversion and within-strand tracing never stitch chains. The separate
`topostich-cgd` graph exporter can optionally stitch **backbone nicks**
(`--stitch-nicks`, NAtopo's P…O3′ join rule) and **chain ends whose bases
pi-stack** (`--stitch-terminal-stacks`, any polarity: 5′-3′, 5′-5′ or 3′-3′)
across chains and symmetry images. It can also write NAtopo's **ASU quotient
graphs** at chain and residue resolution for the covalent bonds and stitched
joins (`--asu-graphs`, `--asu-chain-html`, `--asu-residue-html`). Hydrogen
bonds and other contacts are not stitched.

For agents and incoming collaborators, start with [the component map and next
symmetry-connection task](docs/agent-handoff.json).

## Install

From a checkout:

```sh
python3 -m pip install .
topostich --version
```

Once this local version is published, it can also be installed from GitHub:

```sh
python3 -m pip install git+https://github.com/vecchioni-lab/TopoStich.git
```

For development use `python3 -m pip install -e .`. Prefer a virtual environment.
`python3 -m topostich` is equivalent to the installed `topostich` command.
GitHub publication is distinct from a PyPI release; `pip install TopoStich`
from PyPI is not promised.

## Convert

For the separate **one O3′ node per nucleotide** Systre CGD exporter and
interactive ChimeraX graph views, see [CGD usage](docs/cgd.md):

```sh
topostich-cgd 8D93.pdb --out 8D93_O3.cgd --cmm 8D93_chains.cmm \
  --crystal-cmm 8D93_symmetry.cmm --receipt 8D93_graph.json
```

This keeps chain/TER segments separate. Add `--stitch-nicks` to join
unambiguous backbone nicks and `--stitch-terminal-stacks` to join chain ends
whose bases stack (unless those ends are already backbone joined; 5′-5′ and
3′-3′ stacks are flagged). Joins to symmetry-related or translated copies
become explicit periodic CGD edges (see [nick](docs/cgd.md#backbone-nick-stitching)
and [chain-end stack](docs/cgd.md#chain-end-stack-stitching) stitching). The
[ASU graphs](docs/cgd.md#asu-quotient-graphs) fold every covalent bond and
stitched join back onto one asymmetric unit, labelled with PDB symmetry codes,
and report each part's periodicity:

```sh
topostich-cgd 7SE0.cif --stitch-nicks --stitch-terminal-stacks --out 7SE0_stitched.cgd \
  --cmm 7SE0_stitched.cmm --receipt 7SE0_stitched.json \
  --asu-graphs 7SE0_asu.json --asu-chain-html 7SE0_chains.html --asu-residue-html 7SE0_residues.html
```

`tools/check_stitched_cgd.py` audits a stitched file against its receipt with
Systre's own parser ([native audit](docs/cgd.md#native-systre-audit)). The
following existing conversion commands continue to work unchanged.

Keep every atom in the first model, with its original position:

```sh
topostich 8D93.pdb -atom -all -out 8D93.cif -cif
topostich input.mmcif -atom -all -out all_atoms.cif
```

`-atom -all` performs ordinary file conversion. It includes sugar, base,
phosphate, water, ion and other atoms, and retains alternate locations. Core CIF
expresses the same positions as fractional coordinates; no coordinates are
wrapped or adjusted and no trace links are added. Existing explicit connections
are handled as described below. Omitting `-atom` or using `-atom all` gives the
same result. Do not combine `all` with named atom selections. As before, one
model is written (the first by default); `--model` and `--altloc` remain available.

Make a reduced backbone trace:

```sh
topostich 8D93.pdb -atom "O3'" -atom "O5'" -atom "C4'" -out 8D93_backbone.cif -cif
```

The same command accepts native mmCIF, retaining its coordinate precision:

```sh
topostich input.mmcif -atom "O5'" -atom "C4'" -atom "O3'" -out trace.cif
```

Input syntax is detected from content, including `.cif` files. Use
`--input-format pdb` or `--input-format mmcif` to specify it explicitly.
mmCIF input must contain one macromolecular coordinate data block. Core CIF
and CCD component dictionaries are not inputs for this reader.

Quote atom names containing primes. Atom values do **not** need a leading dash.
`--atom` and `--out` are also accepted. Choose `-pdb`, `-mmCif`, or
`--format cif|pdb|mmcif` for other outputs; the default is core CIF regardless
of the output suffix. Every output must be a new path. Inputs are never edited.

Selected atoms follow **P → O5′ → C5′ → C4′ → C3′ → O3′**, then the next
residue in input-file order. CLI flag order and atom-record order within a
residue do not change this path. Numbering gaps and insertion codes are retained;
residues are not numerically re-sorted. Legacy `*` atom names match primes.
Selections outside that linear path are rejected instead of inventing an order.
Use `-atom -all` or omit `-atom` for ordinary conversion of all atoms, including
non-backbone atoms.

Missing selected atoms are skipped and the remaining selected atoms are joined
within the same chain/TER segment. A wholly missing selected residue is bridged
if the surrounding coordinates belong to that segment. Missing counts are
reported for nucleotide-like residues represented in the input; absent residues
cannot be counted from coordinates alone. Chain changes, `TER`, and model
boundaries are never bridged. Link lengths are the actual endpoint distances;
no chemical distance cutoff, minimization, or ideal geometry is applied.

Without an atom filter, TopoStich preserves explicit input `CONECT` edges and
unambiguous same-image mmCIF `_struct_conn` edges. Unsupported symmetry-image
or unresolved mmCIF links are omitted with receipt warnings. No additional
bonds are guessed. Preserved edges can include cross-chain connections
already declared in the input; TopoStich does not create new ones.

Use `--model 2` for the second model in file order. Ordinary conversion retains
alternate locations. A trace containing selected alternate atoms requires an
explicit `--altloc A` (blank/common atoms plus alternate A) to avoid connecting
mutually exclusive coordinates.

mmCIF author atom/residue/chain identifiers are retained, with label fields as
fallbacks. `label_asym_id` changes split traces even when author chain IDs
repeat. Insertion codes and distinct label residue identities are retained.
PDB output checks field widths rather than truncating long chain IDs or
overflowing values. It rounds coordinates to 0.001 Å and occupancy/B values to
0.01 as required by PDB; core CIF and mmCIF retain native precision.

## Nucleotide recognition underneath the trace

`topostich.nucleotides` recognizes conventional sugar and phosphate heavy-atom
names and looks up base atoms by **residue code → canonical role → actual atom
name present in the structure**. Its packaged snapshot comes from the
NARestraints workbook: 600 residue codes and 619 source rows. No live Excel
or NARestraints installation is needed.

```python
from topostich.nucleotides import recognize_residue

info = recognize_residue("DZ", ["C1'", "C4'", "C1", "N3"])
print(info.present_sugar_atoms)  # ("C1'", "C4'")
print(info.present_base_atoms)   # ("C1", "N3")
```

`classify_atom` reports sugar, phosphate, base, unknown or ambiguous roles.
`recognize_residue` reports observed/missing mapped atoms and candidate workbook
definitions. Observed atom names select the best matching row, as in NARestraints.
Tied conflicting mappings remain explicit. Unknown residues are not assigned a
guessed canonical base. The workbook is an atom-role lookup, not a bond dictionary.

Selecting a known nucleobase atom gives a specific rejection. `-atom -all` or
omitting `-atom` still converts all atoms. Extra sugar atoms C1′, C2′, O2′ and O4′
are recognized internally but remain outside the CLI's linear trace; no ring or ribbon
connections are added. See [recognition and provenance](docs/nucleotide-recognition.md).

## What “CIF” means here

* **Core CIF** uses crystal cell/symmetry, fractional coordinates and tags such
  as `_atom_site_fract_x`. This is the intended default.
* **PDBx/mmCIF** uses macromolecular categories such as `_atom_site.Cartn_x`.
  Changing a `.cif` filename does not convert between dictionaries.
* **CCD** describes chemical components and is itself distributed in mmCIF.
  It is not a whole-crystal core CIF or a refinement restraint dictionary.

Core output requires real cell dimensions and a space-group symbol in PDB
`CRYST1` or the corresponding mmCIF categories.
Missing cell information or the PDB's unit-cell placeholder is an error; TopoStich
does not invent a P1 box. Coordinates are fractionalized without wrapping into
the unit cell. Conventional, rounded `SCALE` records use the more precise
cell-derived transform; nonstandard rotations/translations use the supplied
complete `SCALE`. Full-precision mmCIF matrices and nonzero shifts are retained.
PDB's `H` lattice notation becomes `R` with the same hexagonal
cell. Fractional symmetry operators in `REMARK 290` are retained. If those
operators are absent, the space-group symbol is retained and this is reported;
the downstream reader must resolve its operators.

Core CIF includes `_geom_bond_*` and IUCr `_topol_*` records for the same explicit
graph. Virtual trace edges are **generic links**, with no chemical bond order.
mmCIF uses `_struct_conn` with unspecified connection type and explanatory
details; PDB uses `CONECT`. These derivatives are for coordinates/graphs, not
chemical restraint generation or deposition. Arbitrary PDB metadata, anisotropic
displacement tensors, assembly/MTRIX expansion, and original atom serials are
not carried through; isotropic B values, occupancies and atom identities are.

ToposPro 5.5.2.0 release notes describe CIF connectivity import. **A live
TopoStich → ToposPro import has not yet been validated.** Inspect the imported
adjacency before running automatic bond detection, which may replace the
deliberately long links. Core-CIF readability and graph retention are separate
checks. Keep your full model when studying spatial entanglement: a reduced
straight-line trace need not preserve every geometric feature.

## Machine-readable output and project intent

```sh
topostich input.pdb -atom "O3'" -atom "O5'" -out trace.cif --json > receipt.json
```

The schema-1 receipt contains version, selected model, canonical atom selection,
counts, missing atom identities, warnings and input/output SHA-256. It records
filenames rather than machine-specific operational paths. The CIF carries the
input hash and a comment mapping each short site label to its source identity.
Exit status is 0 on success and 2 on input/output errors.

Read [project intent](docs/intent.json), [architecture](docs/architecture.md),
[contributor contract](AGENTS.md), and [validation handoff](docs/development-handoff.md).
**Every substantial patch must remain installable:** test, build wheel and source
distribution, and install/test the built wheel in a clean environment.

## Tests

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m unittest discover -s tests -v
python3 tools/check_install.py
```

Tests generate synthetic coordinates in temporary directories and require no
downloads. The installation check needs setuptools and wheel in its build
environment and installs the result and Gemmi into a separate temporary virtual
environment. For an offline check, cache compatible dependency wheels and use
`python3 tools/check_install.py --wheelhouse /path/to/wheels`.

## Format references

* [wwPDB CRYST1 and SCALE specification](https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html)
* [IUCr core CIF fractional coordinates](https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_fract_.html)
* [wwPDB Chemical Component Dictionary](https://www.wwpdb.org/data/ccd)
* [IUCr topology dictionary](https://www.iucr.org/__data/iucr/cifdic_html/3/TOPOLOGY_CIF/CTOPOL.html)
* [ToposPro release notes](https://topospro.com/software/topospro/release/)
* [Public example 8D93](https://www.rcsb.org/structure/8D93)

No NASolve datasets or generated crystallographic runs are distributed here.
A software license has not yet been selected by the project owner.
