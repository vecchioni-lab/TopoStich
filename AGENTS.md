# TopoStich contributor contract

Simon requested on 2026-11-08: **TopoStich must remain installable after every
substantial patch.** This is a continuing requirement, not a release-only task.

Before edits, inspect `git status --short`, `README.md`, `docs/intent.json`,
`docs/architecture.md`, `docs/development-handoff.md`, and `CHANGELOG.md`.
The active user request takes precedence over old handoffs. Preserve unrelated
dirty files and all original scientific inputs.

Read `docs/agent-handoff.json` for the component map and next task. Simon clarified
on 2026-11-11 that the collaborator already has the connection-finding logic
engine. The next integration should consume its symmetry-aware backbone links
and write them in the existing CGD dialect; do not invent a replacement engine.
Current 0.3.0 behavior remains unstitched until that integration is implemented.

## Scope and scientific behavior

- The product converts PDB/PDBx-mmCIF to core CIF, optionally selecting a linear
  DNA/RNA backbone trace. PDB and PDBx/mmCIF are alternative coordinate outputs.
- Respect chain, TER and model boundaries. No stitching, symmetry seam joining,
  sticky-end inference, distance-based bond guessing, or chemical repairs are
  authorized by the initial scope. Missing selected atoms are bridged within
  their existing chain/TER segment.
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
