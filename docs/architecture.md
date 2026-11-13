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

`stitching.py` holds the opt-in backbone nick rule ported from NAtopo
(`nucleotide_interactions.find_backbone_join_candidates`, with the loose-end
rules of `prepare_nucleotides` and the polymer links of `build_atom_graph`).
It is called only by `export_cgd` with `stitch_nicks` or `stitch_terminal_stacks`. Per residue it finds
loose P (with own O5′), loose O3′ (with own C3′) and phosphate-less O5′ ends; a
consecutive O3′–P pair at 0.8–2.2 Å or an explicit input link at 0.8–2.3 Å
makes an end non-loose (mmCIF `_struct_conn` links count only with a covalent
`conn_type_id`; untyped PDB CONECT counts). Every space-group operation and every lattice
translation that can bring a loose O3′ within the cutoff of a 5′ end is
enumerated exactly (no finite cell window). Candidates on the same chain copy
are excluded. Joins are stitched only when one-to-one and between segment
termini; the source is the 3′ end in its stored position and the target is the
5′ end's first O3′ node in image `g = operation + n`. The reciprocal CGD row
uses `g⁻¹`, and degrees count distinct neighbor images. Stitched receipts use
schema 2; unstitched output is byte-identical to 0.3.0.

Chain-end stack stitching (`stitch_terminal_stacks=True`) uses the same
`StitchContext`, so nicks are registered first: a stack between two ends that
a stitched nick (or a verified same-image O3′–P link) already joins is
`backbone_joined`, and duplicates of existing edges are detected. Chain ends are
the first (5′) and last (3′) backbone-bearing residues of segments with at
least two of them. Base planes are SVD fits (a dependency-free Jacobi
eigen-solve) of complete ring heavy atoms, with ring names from the packaged
NARestraints catalog. Stacks are enumerated as symmetry orbits: end pairs
`a ≤ b`, every operation and lattice translation within the 5.5 Å centroid
cutoff, and one of `g`/`g⁻¹` for an end against itself. NAtopo's angle,
separation and shift limits apply. A 5′-3′ stack is oriented 3′ → 5′; each end
maps to its segment's outermost O3′ node. CGD rows are de-duplicated per node,
so a join to the node's own image under an involution is one neighbor. The
viewer draws nick joins and chain-end stacks in separate models with explicit
endpoint markers. `tools/CheckCgd.java --stitched` plus
`tools/check_stitched_cgd.py` compare Systre's parsed adjacency with the receipt.

`asu_graphs.py` ports NAtopo's `asu_quotient.py` and `residue_quotient.py`
without NetworkX/NumPy. Operations are exact integer pairs `(R, t·24)`. Raw
residue joins are the segment-internal O3′–P polymer links, deposited covalent
links (0.8–2.3 Å, covalent `conn_type_id` or untyped CONECT, no metal/H/water
endpoint) and the `StitchContext.stitched_links` recorded by `add_edge` (the
actual 3′/5′ residues, not their O3′ nodes). Each join is oriented as
`min((s, t, g), (t, s, g⁻¹))` by node order and operation key and merged per
`(source, target, g, type)`; chain-resolution edges skip same-segment identity
joins (counted as `internal_covalent_edges`). Symmetry codes are written
`op_(k,l,m)` with signed cell shifts (the PDB `op_klm` form is kept as
`pdb_symmetry_code`) and use the input's
operator numbering (`_space_group_symop`/REMARK 290) or Gemmi's order. The
periodicity is NAtopo's voltage-graph method: BFS spanning-tree potentials,
cycle elements, Schreier generation of the translation subgroup, its Hermite
basis, and the index for interpenetrating copies. `asu_view.py` renders
NAtopo's two layouts as inline SVG with a small script for filters, details and
pan/zoom; it uses no external resources. `asu_png.py` reuses those layouts and
edge geometry to draw PNG images with matplotlib (optional `png` extra, imported
only when a PNG is requested; `ARROW_SIZE` is shared by both renderers).

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
