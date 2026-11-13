# TopoStich contributor contract

Simon requested on 2026-11-08: **TopoStich must remain installable after every
substantial patch.** This is a continuing requirement, not a release-only task.

Before edits, inspect `git status --short`, `README.md`, `docs/intent.json`,
`docs/architecture.md`, `docs/development-handoff.md`, and `CHANGELOG.md`.
The active user request takes precedence over old handoffs. Preserve unrelated
dirty files and all original scientific inputs.

Read `docs/agent-handoff.json` for the component map and next task. The
collaborator's connection logic lives in the NAtopo project. 0.4.0 ports its
P…O3′ backbone-join rule into `stitching.py` (opt-in `topostich-cgd
--stitch-nicks`, requested 2026-11-11) and its folded-ASU chain/residue views
into `asu_graphs.py`/`asu_view.py` (covalent bonds and stitching only,
requested 2026-11-13). NAtopo's junction rule was removed on 2026-11-13 in
favour of chain-end stacking (`--stitch-terminal-stacks`): every pi stack
between two chain ends, any polarity, unless already backbone joined. Keep the
ported rules aligned with NAtopo rather than inventing replacements. Audit
stitched output with `tools/check_stitched_cgd.py`.

## Scope and scientific behavior

- The product converts PDB/PDBx-mmCIF to core CIF, optionally selecting a linear
  DNA/RNA backbone trace. PDB and PDBx/mmCIF are alternative coordinate outputs.
- Respect chain, TER and model boundaries. The coordinate converter and the
  default CGD export never stitch. The only authorized stitching is opt-in in
  `topostich-cgd`: NAtopo-rule backbone nick joins (`--stitch-nicks`) and
  pi stacks between chain ends (`--stitch-terminal-stacks`; NAtopo stack
  geometry, any end polarity, same-polarity flagged). Both report unstitched
  candidates with a status instead of choosing silently. No
  hydrogen-bond stitching, general bond guessing, or chemical repairs. Missing
  selected atoms are bridged within their chain/TER segment.
- Keep crystal coordinates, symmetry and graph connectivity distinct. Do not
  label a format check as a live ToposPro compatibility test.
- Keep the package independent of NASolve/Phenix/Coot. Do not add a dependency
  for functionality already adequately served by the standard library.
- Gemmi owns CIF syntax parsing. The packaged NARestraints-derived atom lookup
  maps residue codes and observed names; it is not a chemical bond dictionary.
  Keep unknown/ambiguous mappings explicit. Sugar-ring graphs/ribbons remain
  outside the current trace scope even though extra sugar atoms are recognized.
- Never overwrite inputs or prior outputs. Do not commit scientific coordinates,
  reflections, generated runs, secrets, or machine-local environment state.
  Exception requested 2026-11-13: the 8D93 reference outputs in
  `examples/8D93/` (generated files only, no input structure; not packaged).
- Preserve explicit CLI paths and valid shell quoting for atom-name primes.

## Required checks after each substantial patch

1. Add a regression test and observe it fail before the implementation change.
2. Run focused tests, then the complete suite:
   `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python3 -m unittest discover -s tests -v`.
3. Run `python3 tools/check_install.py`. This must build a source distribution,
   build a wheel from its extracted source after removing the original staging
   tree, install the wheel in a fresh virtual environment, then run the installed
   CLI and the test suite without source-tree PYTHONPATH or editable metadata.
4. Update the changelog, machine-readable intent/schema when affected, and the
   concise validation handoff. Record exact commands, runtime and results;
   don't claim unrun external checks passed.

Every handoff states behavior before/after, files and schemas, validation,
compatibility/migration, excluded generated artifacts and remaining limitations.
Keep new work in isolated branches named `codex/<arboreal-name>-<purpose>`.
