# Development handoff

## 0.4.0 nick and junction stitching — 2026-11-11

Branch `codex/rowan-nick-stitching`. The user asked for backbone nick stitching,
then junction stitching, using the tools already built in the NAtopo project,
plus a Systre audit that works on stitched output.

**Before/after.** 0.3.0 wrote unstitched strands only. 0.4.0 adds opt-in
`topostich-cgd --stitch-nicks` and `--stitch-junctions [--junction-scope face|all]`.
`src/topostich/stitching.py` ports NAtopo's backbone-join screen and its
pi-stack/`nonbackbone_pi_stack_v2` junction rule without NAtopo, NetworkX,
NumPy or SciPy imports. Gemmi remains the only dependency. It searches every
operation and lattice translation exactly. Nicks are stitched only when
one-to-one between segment termini; junctions are stitched when they are face
junctions (default) or any non-backbone stack (`all`). Defaults and all
existing commands are unchanged.

**Files.** New: `stitching.py`, `tests/test_stitching.py` (12),
`tests/test_junctions.py` (8), `schemas/cgd-receipt-2.schema.json`,
`tools/check_stitched_cgd.py`. Changed: `cgd.py` (API/CLI, de-duplicated
neighbor rows/degrees, schema-2 fields), `graph_viewer.py` (nick-join and
junction-stack models), `tools/CheckCgd.java` (`--stitched` mode; default mode
unchanged), `tools/check_install.py` (installed stitched smoke checks), version
files, README, AGENTS, architecture, `cgd.md`, intent, agent handoff, changelog.
Schema 1 is unchanged and still used for unstitched receipts.

**Validation** (Python 3.11.2, Gemmi 0.7.5, macOS arm64, scratch venv):

- Fail-before: the stitching and junction modules failed before their
  implementation. The 2-fold self-junction regression fails with the
  pre-fix row logic (declared degree 3) and passes after the fix.
- `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python -m unittest discover -s tests -v`:
  **68 passed** (48 previous + 20 new).
- `python tools/check_install.py`: PASS, including installed nick and junction
  smoke runs. The installed suite skips one jsonschema assertion because the
  clean venv has no jsonschema; all real receipts were validated separately.
- NAtopo cross-validation (NAtopo venv, Python 3.14, radius-2 window,
  `bases.xlsx`) on 17 nucleic-acid crystals: 1BNA, 1UHX, 1V3P, 5KEK, 5KEO,
  6YL5, 6YLB, 7SE0, 7U3R, 8EPG, 9Q5J, 9Y3B, 9ZXK and 9OSM, plus three local
  refinement PDBs. The results were **identical** to NAtopo:
  - **Nicks:** the same candidate orbits, relative images and distances (≤10⁻⁴ Å)
    as `find_backbone_join_candidates`.
  - **Junctions:** the same candidate orbits, images, face/non-face reasons and
    centroid distances as `build_interaction_templates` → `find_junction_candidates`.
  - **Excluded stacks:** the same count of backbone-adjacent stacks excluded.
- Native Systre 19.6.0 (jar SHA-256
  `0d272e98a1a21669bc67a809b95c014ba2a2fb39a6fd2147039201216ab0d48d`, size
  matches the GitHub release asset, which publishes no digest; system JDK 26.0.1).
  `tools/check_stitched_cgd.py` passed on all **85 exports** (17 crystals ×
  unstitched, nicks, junctions, both, all+O5′ fallback). Systre edge counts
  matched the receipts exactly. All receipts validated against schema 1/2.
  Baseline `CheckCgd` still reproduces 0.3.0's 7SE0 result (126 nodes, 114 edges).
- Default (unstitched) CGD, both CMM views and receipts are byte-identical to
  0.3.0 for 7SE0, 5KEK, 9Y3B and the 8A PDB.

**Found by the audit.** A base stacked on its own 2-fold mate (9Y3B, 6YL5,
6YLB) produced identical forward and reciprocal rows. TopoStich declared
degree 3 where Systre saw 2. Rows are now de-duplicated per node.

**Not run.** ToposPro import; opening the new CMM models in ChimeraX.

**Remaining.** Reviewed confirm/remove decisions for ambiguous nicks and for
junctions (NAtopo's decision log could feed this); special positions; live
ToposPro; scientific review of which junction scope fits each crystal.

## Agent onboarding and next task — 2026-11-11

`docs/agent-handoff.json` maps the components and defines the next integration:
the collaborator already owns the logic engine that finds symmetry-related
backbone connections. TopoStich will map those decisions to node/image identities
and serialize the exact CGD dialect, with matching viewing and provenance.
The engine interface remains to be agreed; the suggested contract is explicitly
a proposal. It flags current edge-schema, viewer and baseline-audit limitations.
README, AGENTS, intent and changelog now point to that task. This update changes
documentation only: runtime/schema versions stay unchanged, generated artifacts
are excluded, JSON/references/whitespace are checked, and the runtime/install
results below are previous implementation results rather than newly rerun tests.

## Behavior before and after

0.2.1 had coordinate conversion, all-atom mode and optional backbone traces.
0.3.0 adds a separate `topostich-cgd` command and module. It writes one observed
O3′ node per nucleotide in the exact CRYSTAL named-source/neighbor-image grammar.
Consecutive retained nodes connect only inside the original segment. Missing
O3′ atoms are skipped/bridged as before; chain, TER, model and mmCIF label-asym
boundaries remain respected. Crystal cell and the full named group are retained.
No nicks, junctions, sticky ends or symmetry copies are stitched.

Companion CMM views use explicit markers and links, chain colors, larger end
markers, source identity notes and separately toggleable models. One view is
the asymmetric unit; another applies each conventional-cell symmetry operation
to intact traces. Existing CIF/PDB/mmCIF behavior and lookup catalog are unchanged.
See `cgd.md` for CLI usage, exact dialect, interpretation and viewing commands.

## Changed files and schemas

- New `src/topostich/cgd.py`: API/CLI, symmetry/position checks, graph receipt,
  exact CGD writer and exclusive multi-output writes.
- New `src/topostich/graph_viewer.py`: native ChimeraX marker/link XML writer.
- New `tests/test_cgd.py`: eight regressions covering boundaries, missing nodes,
  original coordinates, mmCIF precision, alternate/model selection, symmetry,
  CMM links/identities, no-overwrite behavior, CLI and receipt output.
- New `schemas/cgd-receipt.schema.json`: separate schema **1**, self-contained
  node/edge graph, original identities, cell, operations and artifact hashes.
- New `tools/CheckCgd.java`: optional native Systre parser/degree/adjacency audit.
- `tools/check_install.py`: installed CGD command, receipt and CMM smoke checks.
- `pyproject.toml`, `src/topostich/__init__.py`: 0.3.0 and new entry point.
- `MANIFEST.in`: include the optional Java audit source in the sdist.
- README, changelog, architecture, intent, `cgd.md`, and this handoff: updated.

Existing conversion receipt schema 1 and nucleotide catalog schema 1 are
unchanged. No migration, new Python runtime dependency or NASolve change.

## Exact validation and results

Runtime: Python 3.11.9, Gemmi 0.7.5, setuptools 65.5.0, wheel 0.45.1,
macOS arm64. Native checks: Systre 19.6.0, Temurin OpenJDK 21.0.13+9,
ChimeraX 1.12 (2026-06-12). Java was downloaded and checksum-verified into
scratch storage; it is not a package dependency or system installation.

`TEST_PYTHON` is the isolated interpreter and `WHEELHOUSE` the offline Gemmi
wheel directory. Exact host paths, commands and logs are retained locally in
`local/8D93_O3_CGD_2026-11-11/validation/validation.json`; operational artifact
references there use the bundle anchor and relative paths.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src "$TEST_PYTHON" -m unittest discover -s tests -p test_cgd.py -v
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src "$TEST_PYTHON" -m unittest discover -s tests -v
"$TEST_PYTHON" tools/check_install.py --wheelhouse "$WHEELHOUSE"
javac -cp "$SYSTRE_JAR" -d "$AUDIT_DIR" tools/CheckCgd.java
java -cp "$SYSTRE_JAR:$AUDIT_DIR" CheckCgd "$BUNDLE/8D93_O3_chain_aware.cgd"
"$CHIMERAX" --nogui --exit --script "$AUDIT_DIR/check_chimerax.py"
```

- Fail-before: new regression module failed importing the nonexistent exporter.
- Focused suite: **8 passed**. Complete source suite: **48 passed**, retaining
  all 40 previous regressions.
- Install check: sdist built; staging source deleted; wheel built from relocated
  extracted source; installed in a fresh venv offline. `pip check`, old/new
  installed commands, coordinate-format smoke checks, packaged catalog and
  all **48 installed tests** passed without source PYTHONPATH.
- CGD receipt schema checked and 8D93 receipt validated with jsonschema 4.26.0.
- Native Systre parse: **126 primitive-cell nodes, 114 edges**; declared degrees
  agree (24 endpoints, 102 internal). Every edge joins consecutive labels in
  one segment. Twelve finite components: three copies of each size 7/9/12/14.
- Systre's cell warning was diagnosed using its own metric checker: floating
  cosine values give error 1; exact 90°/120° cosines give error 0. Metric
  difference norm is 7.95e-12 Å². The cell was neither changed nor rescaled.
- Native ChimeraX complete file-open checks: ASU **4 models, 42 markers,
  38 links**; symmetry view **36 models, 378 markers, 342 links**. All marker
  identity notes loaded. GUI renders completed and were visually inspected;
  desktop automation itself was unavailable. No Phenix/Coot run was required.

## Artifacts and remaining limits

8D93 comes from the existing checksum-verified full-backbone batch manifest.
Original source is unchanged. The ignored local bundle contains source copy,
CGD, two CMM views, JSON mapping, PNG previews, instructions and validation
logs. None of these scientific artifacts or downloaded validator/runtime
binaries are committed or included in the package distribution.

The current graph is a collection of finite strands, not a connected
three-periodic net. ToposPro import still needs collaborator testing. Symmetry
copies are displayed but not joined. Special-position/stabilizer handling,
ambiguous symmetry overlaps and endpoint-image/junction rules remain open.
Cells without trustworthy group metadata are rejected. Native parser checking
is optional and does not run on ordinary conversions. Continue keeping the
package installable after each substantial patch.
