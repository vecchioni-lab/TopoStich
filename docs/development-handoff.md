# Development handoff

## TopoTrace rename — 2026-11-13

TopoStich 0.3.0 is renamed to **TopoTrace 0.3.1**. Canonical package/commands are
`topotrace` and `topotrace-cgd`; old `topostich` commands, executable modules and
public imports remain aliases of the same implementation. Scientific behavior,
CLI options and conversion/CGD/catalog schema 1 remain unchanged. Branding and
version strings in newly written outputs change their bytes/hashes.

Changed files: package moved to `src/topotrace`, compatibility wrappers retained
in `src/topostich`; pyproject/MANIFEST, install checker, tests, README, AGENTS,
architecture/CGD/recognition docs, JSON intent/component map, schema titles and
changelog updated. Agent handoff document is version 2 because product-specific
responsibility keys changed. The atom catalog bytes are unchanged.

Migration: uninstall the old TopoStich distribution before installing TopoTrace
so two distributions do not own the same compatibility command/module files.
Existing structures and receipts need no migration. GitHub target:
https://github.com/vecchioni-lab/TopoTrace. PyPI publication is not claimed.

Validation: Python 3.12.14, Gemmi 0.7.5, setuptools 84.2.1, wheel 0.48.0,
macOS arm64. `PYTHON` below was the temporary isolated interpreter and
`WHEELHOUSE` the temporary downloaded dependency directory (not project inputs).

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src "$PYTHON" -m unittest discover -s tests -p test_rename.py -v
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src "$PYTHON" -m unittest discover -s tests -v
"$PYTHON" tools/check_install.py --wheelhouse "$WHEELHOUSE"
git diff --check
```

- The new two-test rename regression failed before implementation: topotrace
  was absent and legacy commands still reported TopoStich 0.3.0.
- Focused: **2 passed**; complete source suite: **50 passed**.
- Install: sdist built, original staging removed, wheel built from relocated
  extracted source and installed in a fresh environment offline. Both old/new
  command names, module aliases, public API identity, packaged catalog and
  format/CGD/CMM smoke checks passed; **50 installed tests passed**.
- JSON parses, component paths and whitespace validated. No live ToposPro,
  Systre or ChimeraX run was needed or claimed for this naming-only patch.

Scientific inputs, previous generated outputs, environments and local bundles
are excluded from the commit. Symmetry stitching remains the collaborator-engine
integration task; the rename adds no connectivity or scientific inference.

## Historical 0.3.0 implementation handoff

The following is prior evidence under the original TopoStich name. Its
publication-pending wording describes that earlier point in time.

Updated 2026-11-11. Local package **0.3.0**, implemented on
`codex/linden-cgd-export`; documentation handoff on `codex/linden-agent-handoff`.
GitHub publication is pending. No remote CI, release, ToposPro import or net
identification is claimed.

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
