# Changelog

## [Unreleased]

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
