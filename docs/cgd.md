# Chain-aware O3′ graphs for Systre and ToposPro

TopoStich 0.3.0 adds a separate `topostich-cgd` command. Each observed O3′
becomes one nucleotide node; consecutive retained nodes link within their
original chain/TER segment. This is a topological trace, not a chemical model.
Existing `topostich` conversions remain unchanged.

```sh
topostich-cgd 8D93.pdb --out 8D93_O3.cgd \
  --cmm 8D93_chains.cmm --crystal-cmm 8D93_symmetry.cmm \
  --receipt 8D93_graph.json
```

Native PDBx/mmCIF input works through the same command. `--model`, `--altloc`,
and `--input-format` retain the converter's meanings. `python -m topostich.cgd`
is equivalent to the installed command. All output paths must be new and have
existing parents. The command does not overwrite inputs or earlier outputs.

## Exact dialect

The output uses the screenshot's Systre/Gavrog **CGD CRYSTAL** grammar:

```text
crystal
name example_O3
cell 20 30 40 90 90 90
group P1
atom S1_N1 1 0.05 0.1 0.1
edge S1_N1 0.2 0.2 0.2
atom S1_N2 1 0.2 0.2 0.2
edge S1_N2 0.05 0.1 0.1
end
```

Each `atom` declares a unique label, degree, and fractional position. Each
`edge` names its source and gives its neighbor's fractional image coordinates.
Both incident directions are written, matching the supplied examples; Systre
deduplicates these into undirected edges. No alternative EDGE grammar is used.
Coordinates are never wrapped. The cell and named space-group setting are
retained; `R3:H` explicitly denotes R3 in hexagonal axes, as in PDB `H 3`.

Labels `S<segment>_N<ordinal>` identify nodes; the JSON receipt maps them to
original chain, residue, insertion, alternate and mmCIF label identities.
Degrees come from the explicit graph (usually 1 at ends, 2 internally).
Chain/TER/label-asym boundaries remain open. Missing O3′ atoms are skipped and
bridged within the existing segment, with missing selections reported.

The named group implies full symmetry expansion. Explicit input operators,
when present, must agree with that setting. No link joins copies or chain ends.
Near symmetry-equivalent nodes and special positions are rejected in this first
version so the reader cannot silently merge identities or change degrees.
The conservative collision guard uses 0.003 in each conventional fractional
coordinate, modulo lattice translations, against all group images. Systre's
own matcher uses 0.001 in its primitive setting. Future support for special
positions must account for site stabilizers and expanded edge multiplicities.

## Backbone nick stitching

`--stitch-nicks` (API: `export_cgd(..., stitch_nicks=True)`) joins strands
across backbone nicks. The rule is the NAtopo project's backbone-join screen,
ported without NAtopo, NetworkX or SciPy dependencies:

| Step | Rule (defaults from NAtopo's `nucleotide_interactions.json`) |
|---|---|
| Loose P | Residue has P and its own O5′; P has no covalent link to another residue |
| Loose O3′ | Residue has O3′ and its own C3′; O3′ has no covalent link to another residue |
| Existing link | Consecutive residues in a segment with O3′–P 0.8–2.2 Å, or an input CONECT/same-image `_struct_conn` at 0.8–2.3 Å |
| Backbone join | Loose P…loose O3′ **< 3.5 Å** (`--join-cutoff`) on a different chain copy: another segment, or any other symmetry/lattice image |
| O5′ fallback | Phosphate-less residue's O5′ (own C5′) …loose O3′ **< 6.0 Å** (`--o5-cutoff`); reported only unless `--include-o5-fallback` |

Images are searched exactly: every conventional space-group operation
(including centering) and every lattice translation that can bring the
endpoints within the cutoff. NAtopo's finite cell window is not needed. On 17
nucleic-acid crystals in the NAtopo workspace, TopoStich found exactly the same
candidate orbits, relative images and distances (to 10⁻⁴ Å) as NAtopo's own
`find_backbone_join_candidates` with its radius-2 window.

A candidate becomes a graph edge only when it is **one-to-one** (each end
has a single enabled partner across all images) and joins a segment's 3′
terminus to another copy's 5′ terminus. An input O3′–P link that already joins two
segments in the same image is stitched as `deposited_link`. Every candidate is
listed in the receipt's `join_candidates` with one status: `stitched`, `ambiguous`,
`segment_internal`, `o5_fallback_disabled`, `no_graph_node` or `duplicate_edge`.
Nothing is chosen silently.

Each stitched edge goes from the upstream 3′ end's node (stored position) to
the downstream strand's first O3′ node in the image `g = operation + n`:

```text
atom S2_N7 2 ...
edge S2_N7 <g applied to S2_N1>        # forward join row
atom S2_N1 2 ...
edge S2_N1 <g⁻¹ applied to S2_N7>      # reciprocal row: inverse affine image
```

The receipt (schema 2, `schemas/cgd-receipt-2.schema.json`) records
`operation`, integer `lattice_translation`, `image_operation`,
`reverse_image_operation`, both image coordinates (fractional and Cartesian),
distance, endpoint descriptions and the evidence for each join. Degrees count
distinct neighbor images. A join from a node to its own image under an
involution has identical forward and reciprocal rows; it is written once and
counted once. Without stitching options, CGD, CMM and schema-1 receipts are
byte-identical to 0.3.0.

Limits: joins are geometric proposals, not chemical bonds. The same chain copy
cannot join itself in its own image (NAtopo's interchain rule).

## Junction stitching

`--stitch-junctions` (API: `stitch_junctions=True`) joins strands where base
stacking continues a helix without a backbone bond, for example across a nick
whose 5′ phosphate is missing, or where helix ends stack across a crystal
contact. The rule is NAtopo's `nonbackbone_pi_stack_v2` junction classifier
over its pi-stack screen (`find_pi_interactions`, `fit_base_planes`,
`junctions.find_junction_candidates`):

| Step | Rule (NAtopo defaults) |
|---|---|
| Base plane | SVD plane of the complete ring heavy atoms: 9 for purines, 6 for pyrimidines; centroid = ring-atom mean |
| Ring atom names | TopoStich's NARestraints-derived residue catalog maps canonical ring positions to actual names, including modified residues; unmapped/ambiguous residues are reported in `plane_diagnostics` |
| Pi stack | Centroids **< 5.5 Å**, normals within **30°**, both normal separations **2.0–4.5 Å**, both lateral shifts **≤ 3.0 Å** |
| Junction candidate | A stack whose residues are not joined by a verified O3′–P bond (same image) |
| Face reason | `stack_on_face_without_backbone_continuation` if, for either endpoint, every verified backbone neighbor lies on the opposite face (0.5 Å side tolerance); otherwise `non_backbone_stack` |

`--junction-scope face` (default) stitches only the face reason: the helix
continues by stacking past a strand end. `--junction-scope all` stitches every
junction candidate, including stacks on a face that already has backbone.
Stacks are not one-to-one by nature, so every qualifying stack is stitched.
A stack between the same two nodes in the same image as an existing trace or nick
edge is reported as `duplicate_edge`. Each candidate appears in
`junction_candidates` with `stitched`, `not_face_junction`, `no_graph_node` or
`duplicate_edge`, its measurements and both endpoints' backbone context. The
edge joins the two residues' O3′ nodes; `distance_A` is the ring-centroid distance.

On the same 17 crystals, TopoStich's junction candidates, relative images,
face/non-face reasons, centroid distances (to 10⁻⁴ Å) and counts of excluded
backbone stacks are identical to NAtopo's full `build_interaction_templates`
pipeline (radius-2 window, `bases.xlsx`).

Both CMM views draw stitched joins in separate models: gold **nick joins** and
magenta **junction stacks**. The ASU view draws each join from its source node
to an explicit marker at the target image. The symmetry view draws it for every
operation, so join endpoints can lie in translated copies outside the
displayed cells.

| Structure | Nick joins | Face junctions | Systre components, nicks + junctions |
|---|---|---|---|
| 1BNA | 0 | 0 | 8 finite |
| 1UHX | 0 | 1 | 6 × 1-periodic |
| 1V3P | 0 | 2 | 4 × 1-periodic |
| 5KEK / 5KEO | 2 | 4 | 6 × 1-periodic, 2 × 3-periodic |
| 6YL5 / 6YLB | 0 | 5 | 33 finite |
| 7SE0 | 2 | 4 | 3 × 1-periodic, 1 × 3-periodic |
| 7U3R | 2 | 4 | 4 finite |
| 8EPG | 2 | 2 | 2 × 3-periodic |
| 9Q5J | 2 | 2 | 5 × 1-periodic |
| 9Y3B | 0 | 8 | 4 finite, 8 × 1-periodic |
| 9ZXK | 0 | 4 | 4 finite, 1 × 1-periodic |
| 9OSM | 1 | 3 | 1 finite, 4 × 1-periodic |
| 82B refinement | 2 | 2 | 4 × 3-periodic |
| 8A refinement | 3 | 11 | 1 × 1-periodic, 1 × 3-periodic |
| arch-solids refinement | 4 | 2 | 4 finite |

Counts are per primitive cell from Systre 19.6.0. They describe the stitched
graph only; they are not a scientific verdict on each crystal's connectivity.
O5′ fallback joins and `all`-scope junctions change several of these.

## Native Systre audit

`tools/CheckCgd.java` opens a CGD with Systre's own parser. Its default mode is
the strict **unstitched** baseline: every edge joins consecutive nodes of one
segment and every component is finite. `--stitched` keeps the Systre-side checks
(parse, declared degrees) and prints Systre's primitive-cell edges and
components. `tools/check_stitched_cgd.py` compiles and runs it, then compares
with the receipt: the same node pairs, and for each pair exactly one edge per
primitive-cell copy (half for a self-reversing join).

```sh
python tools/check_stitched_cgd.py graph.cgd graph.json --systre /path/to/Systre-19.6.0.jar
```

All 85 exports (17 crystals × unstitched, nicks, junctions, both, and
everything including O5′ fallback with `all` junctions) passed. Systre's edge
counts matched the receipts exactly. The audit found one bug during
development: duplicate rows for a junction between a base and its own 2-fold
mate. That bug is fixed and covered by a regression test.

## Explore in ChimeraX

Open the companion **`.cmm`** file in ChimeraX using File → Open. This uses its
native marker/link reader, with the same explicit edges and no bond guessing.
Drag to rotate, scroll to zoom, and toggle individual strands in the Models
panel. Each chain/segment has a color; endpoints have larger spheres.

Control-click a node to select it, then show its original identity with:

```text
label sel atoms attribute marker_note
```

Use `label delete` to clear labels. Marker numbers within each model are the
CGD `N` ordinals, not the source residue numbers; the note and JSON receipt
retain the original identities. Editing a ChimeraX view does not update CGD.

`--cmm` shows the observed asymmetric unit in its original Cartesian frame.
`--crystal-cmm` shows one intact copy per conventional-cell symmetry operation
in the standard cell frame; stitched joins appear only in their own models. This is a finite inspection
window, not every nearby translated cell. Copies are deliberately not wrapped
atom by atom; wrapping would tear strands across cell faces. Each copy/segment
is separately toggleable. In a nonstandard input Cartesian frame the two views
may differ by rigid frame/origin convention, while CGD retains fractional sites.

## Validation and limits

For 8D93, the observed graph has **42 nodes, 38 undirected edges and four
segments**: A=12, B=7, C=14, D=9. There are no missing O3′ selections.
The cell is 105.547, 105.547, 96.268 Å; 90°, 90°, 120°; R3:H.
The symmetry CMM contains 378 nodes and 342 links in 36 independently
toggleable models (nine conventional-cell operations × four segments).

Native Systre **19.6.0** parsed the graph into its primitive cell: 126 nodes,
114 edges, degrees 1/2 as declared, and twelve finite components. Each input
segment has three rotated copies in the primitive cell. The audit checks
consecutive node labels and absence of cross-segment edges. These counts are
different from the conventional-cell CMM because R centering has multiplicity 3.

Systre 19.6.0 emits `Unit cell parameters illegal for this group` for these
valid hexagonal parameters. A direct diagnostic of its own metric checker
reproduced the warning from cosine roundoff: using exact cos(90°)=0 and
cos(120°)=-1 removes it, with a matrix difference norm of only 7.95e-12 Å².
Gemmi also accepts the cell/group combination. The exported cell is unchanged;
the warning is recorded, not suppressed. No Systre topology identification is
claimed: unstitched finite strands are not yet a connected 3-periodic net.

ChimeraX **1.12 (2026-06-12)** opened both CMM files through its native reader
with the expected nodes, links and notes. Live **ToposPro import is untested**.
Nick and junction stitching are available (above). Special-position
handling, reviewed decisions for ambiguous joins, live ToposPro import and
ribbon construction remain future work.

The baseline native parser audit is reproducible with a separately obtained
Systre jar and Java SDK (neither is a runtime dependency):

```sh
javac -cp /path/to/Systre-19.6.0.jar -d /tmp/cgd-check tools/CheckCgd.java
java -cp /path/to/Systre-19.6.0.jar:/tmp/cgd-check CheckCgd /path/to/graph.cgd
```

Sources: [Systre parser](https://github.com/odf/gavrog/blob/Systre-19.6.0/src/org/gavrog/joss/pgraphs/io/NetParser.java),
[ToposPro CGD import](https://topospro.com/software/topospro/manual/creating-a-new-database-from-cif-shelx-or-systre-formats/),
and the installed ChimeraX 1.12 marker-format documentation and native CMM reader.
