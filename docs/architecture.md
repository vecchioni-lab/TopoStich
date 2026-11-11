# Architecture

`cgd.py` provides the separate `topostich-cgd` CLI and `export_cgd` API. It reuses
the established readers, O3′ selection and fractional transform without altering
the original converter. Gemmi validates the named space-group setting and cell;
a conservative symmetry-position guard rejects ambiguous overlaps and special
positions. The graph is serialized as CRYSTAL with named-source neighbor-image
edges. `graph_viewer.py` writes native ChimeraX marker/link XML from that same
graph, either the observed asymmetric unit or intact conventional symmetry copies.
All payloads serialize before exclusive output creation; partial writes are
rolled back on I/O failure. A separate schema-1 CGD receipt maps node identities,
coordinates, edges, cell, operations and artifact hashes. Existing conversion
receipt and nucleotide catalog schemas are unchanged. See `cgd.md` for limits.

`src/topostich/converter.py` owns fixed-column PDB parsing, input dispatch, atom selection,
coordinate transforms, graph construction and the three output writers.
`src/topostich/cli.py` maps CLI arguments to `convert_file` and serializes the
schema-1 receipt. `mmcif_input.py` uses Gemmi for CIF syntax and normalizes native
mmCIF into the shared atom/model container without a PDB intermediate. Conversion
does not invoke external scientific programs or the network.

The CLI accepts `-atom -all` and `-atom all` as explicit aliases for the existing
unfiltered `convert_file(..., atoms=None)` path. It binds the dash-prefixed value
before argument parsing and rejects mixtures of `all` and named selections.
This does not change the converter, per-model behavior or receipt schema:
`atom_selection: []` continues to mean all atoms in the selected model, subject
to an explicit alternate-location choice. No atom-role rejection or backbone
graph construction runs in this mode.

Input sniffing uses content, with `--input-format` as an explicit override.
One mmCIF data block is required. The reader preserves row order and native
numeric precision, resolves auth fields with label fallbacks, and separates
label-asym/author-chain boundaries. Model selection is ordinal in encounter
order. The shared `Atom` retains label residue identity so repeated author
numbers cannot silently merge distinct mmCIF residues.

mmCIF same-image `_struct_conn` records are resolved using their declared
label/auth identifiers and alternate/insertion constraints. Ambiguous
endpoints stop. Unsupported image links or absent endpoints are omitted with
warnings in unfiltered conversion. Trace mode constructs its own within-segment
edges and does not add input connection records. No symmetry stitching or bond
perception is introduced.

The conversion reads a snapshot of source bytes and records its SHA-256.
Parsing retains atom identities including insertion code, alternate location,
and a segment ID that changes at every chain/TER boundary. One model is selected
explicitly or defaults to the first. Backbone filtering groups atoms by residue,
then uses a fixed ordered linear path. Adjacent retained atoms receive an edge
only when they share a segment. No nearest-neighbor calculation is involved.
The unfiltered path uses the input's unambiguous CONECT or same-image
`_struct_conn` records instead.

Core CIF needs fractional coordinates. The conventional orthogonalization
matrix is inverted analytically from the six cell parameters. A complete SCALE
matrix takes precedence for nonstandard coordinate frames. A SCALE equal to
the rounded conventional transform is replaced with the more precise
cell-derived transform. Native full-precision mmCIF transforms and nonzero
shifts are retained. Fractional coordinates are not wrapped, so link images
do not change. Fractional REMARK 290 operators and their IDs are copied, not
generated from a partial private space-group table. PDB H lattice naming is
normalized to R without changing its hexagonal cell or supplied operations.

Core graph output includes conventional geometry rows and IUCr topology rows;
the latter mark edges as generic links. PDB uses fresh serials and CONECT;
mmCIF uses unique segment-based label_asym_id plus original auth identities and
struct_conn rows with unknown chemical connection type. File format checks
do not prove a target application retains the links. Check that separately
against the target version before claiming compatibility.

All validation/serialization completes before exclusive output creation.
Existing paths, including symlinks, are refused. Output is self-contained;
source names/hashes are historical provenance, not operational dependencies.
The optional JSON receipt goes to stdout, leaving the user's choice of where
to save it explicit. No machine-local workspace registry is created.

`nucleotides.py` owns residue-aware atom roles independently of tracing. Its
packaged JSON catalog is built from the NARestraints workbook by an explicit,
read-only development tool. Runtime lookup never opens a host workbook or
imports NASolve/NARestraints. Matching uses the residue code and observed atom
names. Ties, role collisions and unknown atoms remain explicit; the library
does not infer chemical bonds. See `nucleotide-recognition.md` for provenance
and the internal API contract.

Known limits: PDB/PDBx-mmCIF input only; one model per output; linear conventional
nucleic-acid backbone selection; no reconstruction of absent residues, biological
assemblies, or MTRIX copies; no preservation of anisotropic displacement tensors
or arbitrary PDB annotation; no full crystallographic dictionary validation.
Malformed cells/coordinates, incomplete SCALE records, ambiguous selected
alternates/duplicates and empty selections stop with a useful diagnostic.

PDB output from mmCIF validates fixed-field limits before opening the output;
it reports PDB's coordinate/occupancy/B-factor rounding. Native CIF output does
not pass through these fields.

Packaging is a `src`-layout setuptools project with Gemmi as its only runtime
dependency. JSON atom data are included in both source and wheel distributions.
`tools/check_install.py` copies only distributable source into temporary staging,
builds an sdist, deletes that staging tree, builds a wheel from extracted source,
and tests the installed wheel and its dependencies in a clean virtual environment.
An optional wheelhouse enables a fully offline dependency installation. Tests and build
products stay outside scientific working directories.
