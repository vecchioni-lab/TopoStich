# Chain-aware O3′ graphs for Systre and ToposPro

TopoTrace provides a separate `topotrace-cgd` command. Each observed O3′
becomes one nucleotide node; consecutive retained nodes link within their
original chain/TER segment. This is a topological trace, not a chemical model.
Existing `topotrace` conversions remain unchanged.

```sh
topotrace-cgd 8D93.pdb --out 8D93_O3.cgd \
  --cmm 8D93_chains.cmm --crystal-cmm 8D93_symmetry.cmm \
  --receipt 8D93_graph.json
```

Native PDBx/mmCIF input works through the same command. `--model`, `--altloc`,
and `--input-format` retain the converter's meanings. `python -m topotrace.cgd`
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
in the standard cell frame, with no added joins. This is a finite inspection
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
Nick/junction selection, symmetry-image stitching, special-position handling,
automatic sticky-end inference, and ribbon construction remain future work.

The optional native parser audit is reproducible with a separately obtained
Systre jar and Java SDK (neither is a runtime dependency):

```sh
javac -cp /path/to/Systre-19.6.0.jar -d /tmp/cgd-check tools/CheckCgd.java
java -cp /path/to/Systre-19.6.0.jar:/tmp/cgd-check CheckCgd /path/to/graph.cgd
```

Sources: [Systre parser](https://github.com/odf/gavrog/blob/Systre-19.6.0/src/org/gavrog/joss/pgraphs/io/NetParser.java),
[ToposPro CGD import](https://topospro.com/software/topospro/manual/creating-a-new-database-from-cif-shelx-or-systre-formats/),
and the installed ChimeraX 1.12 marker-format documentation and native CMM reader.
