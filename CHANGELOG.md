# Changelog

## [Unreleased]

- Version 0.4.0: opt-in backbone nick and chain-end stack stitching for
  `topostich-cgd`, plus ASU quotient graphs. Nicks
  (`--stitch-nicks`, `--join-cutoff`, `--o5-cutoff`, `--include-o5-fallback`;
  API `export_cgd(stitch_nicks=True, ...)`). New `stitching.py` ports NAtopo's
  loose P…O3′ backbone-join rule (3.5 Å, different chain copy, O5′ fallback
  6.0 Å) with exact search over all symmetry operations and lattice translations.
  Only one-to-one joins between segment termini are stitched; ambiguous,
  internal and fallback candidates are reported with a status.
- CGD writes stitched edges as neighbor-image rows with inverse-affine
  reciprocals and recomputed degrees. New schema-2 receipt
  (`schemas/cgd-receipt-2.schema.json`) for stitched exports. CMM views draw
  joins in a separate model. Unstitched output and schema-1 receipts are
  byte-identical to 0.3.0.
- Opt-in chain-end stack stitching (`--stitch-terminal-stacks`; API
  `stitch_terminal_stacks=True`). Every pi stack (NAtopo's
  `find_pi_interactions` geometry, ring names from the packaged residue
  catalog) between the first/last nucleotides of two segments is stitched in
  any symmetry image and for any polarity, unless the ends are already backbone
  joined (verified O3′–P link or a stitched nick join). 5′-5′ and 3′-3′ stacks
  are flagged (`polarity`, `same_polarity`). Candidates are reported in
  `terminal_stack_candidates`. This replaces the NAtopo junction rule
  (`--stitch-junctions`, `--junction-scope`) drafted earlier in this
  unreleased version; it was never released.
- ASU quotient graphs (`--asu-graphs JSON`, `--asu-chain-html`,
  `--asu-residue-html`; schema `schemas/asu-graphs.schema.json`). New
  `asu_graphs.py` ports NAtopo's `asu_quotient`/`residue_quotient` (the folded
  chain and residue views of `main_NAtopo.ipynb`) for covalent bonds and
  stitched joins only: one edge per symmetry orbit with its PDB `op_klm` code,
  and per-part periodicity (rank of the cycle-generated translations,
  interpenetrating copies). `asu_view.py` writes self-contained SVG/HTML views
  with NAtopo's layouts. The receipt gains an optional `asu_graphs` summary.
- mmCIF `_struct_conn.conn_type_id` is now kept (`Pdb.conect_types`); only
  covalent types (or untyped PDB CONECT) count as deposited covalent links for
  nick endpoints and ASU graphs. Short `hydrog` records in 5KEK/9ZXK had been
  read as covalent. Conversion and default CGD output are unchanged.
- CGD rows are de-duplicated per node: a join to a node's own image under an
  involution is one neighbor (found by the new Systre audit).
- `tools/CheckCgd.java --stitched` and `tools/check_stitched_cgd.py`: native
  Systre audit of stitched output against the receipt. The default Java mode
  remains the strict unstitched baseline.
- `examples/8D93/`: stitched (`_stiched`) and unstitched (`_unstiched`) CGD,
  receipts, CMM views and ASU graphs for PDB 8D93, added on request as
  reference outputs (input not included; not packaged).
- Stitching, chain-end stack and ASU-graph regressions and installed CLI smoke
  checks. Agreement with NAtopo on 17 crystals (nick candidates, ASU covalent
  and join orbits, symmetry codes, periodicity) and the Systre audit of 85
  exports are documented in `docs/cgd.md`.

- Add `docs/agent-handoff.json` with a component map, project intent and the
  next collaborator-engine integration task for symmetry-related backbone links.
  Distinguish implemented behavior from proposed image/identity interfaces;
  document required graph-schema, viewer and validator work. Documentation only.

- Version 0.3.0: add a separate installable `topostich-cgd` command for the
  exact CRYSTAL named-source/neighbor-image dialect. Represent each nucleotide
  by its observed O3′ and preserve chains/TER segments, model selection,
  missing-atom bridging, crystal cell and symmetry setting. No stitching yet.
- Add native ChimeraX CMM graph views, including intact symmetry copies,
  and a distinct schema-1 graph receipt with original nucleotide identities.
  Reject ambiguous symmetry overlaps/special positions before writing.
- Add regression/install checks and an optional native Systre parser audit.
  Existing CIF/PDB/mmCIF commands and receipts remain compatible. Scientific
  outputs and native-tool logs remain local and excluded from distributions.

- Record the agreed future CGD CRYSTAL dialect: one O3′ node per nucleotide,
  explicit named-source/neighbor-image edges, and pending nick/junction rules.
  Design documentation only; current conversion behavior is unchanged.

- Version 0.2.1: accept explicit `-atom -all` (also `-atom all` / `--atom all`)
  for ordinary all-atom conversion of the selected PDB/mmCIF model. Preserve
  positions and existing explicit connections without adding backbone links.
  Reject mixing `all` with named selections. Omitting `-atom`, trace commands,
  model/alternate options and schema-1 receipts retain their existing behavior.
- Add all-atom CLI regressions and an installed-command smoke check; keep the
  package installable through the source-distribution and clean-wheel checks.
- Version 0.2.0: accept native PDBx/mmCIF or PDB, detect syntax from content,
  and support explicit `--input-format` while retaining the existing CLI aliases.
- Preserve native mmCIF precision, model/alternate identities and label-chain
  boundaries. Read same-image explicit connections; report omitted image or
  unresolved links. Check PDB output field limits before creating a file.
- Add residue-aware sugar/phosphate/base recognition using a provenance-tracked
  NARestraints workbook snapshot. Match observed atom names, expose unknowns and
  ambiguous mappings, and give specific base-selection errors. Extra sugar atoms
  are recognized internally; the linear trace and no-stitching policy are unchanged.
- Package the lookup JSON and add Gemmi as a declared runtime dependency; Excel
  and NARestraints remain unnecessary at runtime. Keep schema-1 receipts unchanged.
- Extend regressions and clean installation checks to mmCIF and packaged lookup
  data, with an optional offline dependency wheelhouse.

## [0.1.0] - 2026-11-08

- Initial installable `topostich` command; PDB input and core CIF (default),
  PDBx/mmCIF or PDB output, with no runtime dependencies.
- Optional repeated atom selections follow conventional backbone order and
  bridge missing selections within a chain/TER segment. No stitching is inferred.
- Preserve cell, fractional symmetry records, occupancies and isotropic B values;
  handle nonorthogonal cells and nonstandard SCALE transforms.
- Explicit graph records and schema-1 JSON receipts with input/output hashes.
- Human and machine-readable intent, regression tests and required build/install
  checks after every substantial patch. No prior artifact schemas need migration.

Live ToposPro import remains unvalidated. No scientific datasets are distributed.
