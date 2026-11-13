# Development handoff

## 0.4.0 nick and chain-end stack stitching, ASU graphs — 2026-11-13

Branch `codex/rowan-nick-stitching` (PR #1). The user asked for backbone nick
stitching using the tools already built in the NAtopo project, plus a Systre
audit that works on stitched output. On 2026-11-13 they asked to **remove
junction stitching completely**, to stitch instead **any chain ends whose bases
pi-stack** (unless already backbone joined), for every polarity (5′-3′, 5′-5′,
3′-3′) with same-polarity stacks noted, and to generate NAtopo's **ASU graphs at
chain and residue resolution** (code from `main_NAtopo.ipynb`) for covalent
bonds and stitching only.

**Before/after.** 0.3.0 wrote unstitched strands only. 0.4.0 adds opt-in
`topostich-cgd --stitch-nicks`, `--stitch-terminal-stacks`, and
`--asu-graphs/--asu-chain-html/--asu-residue-html`. The junction rule and its
flags (`--stitch-junctions`, `--junction-scope`) were drafted earlier on this
branch and removed before release. `stitching.py` ports NAtopo's backbone-join
screen and its pi-stack geometry; `asu_graphs.py` ports `asu_quotient.py` /
`residue_quotient.py`; `asu_view.py` renders their layouts as self-contained
SVG/HTML. No NAtopo, NetworkX, NumPy, SciPy or Bokeh imports; Gemmi remains
the only dependency. Defaults and all existing commands are unchanged.

**Files.** New: `stitching.py`, `asu_graphs.py`, `asu_view.py`,
`tests/test_stitching.py` (12), `tests/test_terminal_stacks.py` (10),
`tests/test_asu_graphs.py` (12), `schemas/cgd-receipt-2.schema.json`,
`schemas/asu-graphs.schema.json`, `tools/check_stitched_cgd.py`. Changed:
`cgd.py` (API/CLI, de-duplicated neighbor rows/degrees, schema-2 fields, ASU
outputs), `converter.py`/`mmcif_input.py` (`Pdb.conect_types`,
`covalent_connections`), `graph_viewer.py` (nick-join and terminal-stack
models), `tools/CheckCgd.java` (`--stitched` mode), `tools/check_install.py`,
schema 1 (optional `asu_graphs` summary only), version files and docs.

**Validation** (Python 3.11.2, Gemmi 0.7.5, macOS arm64, scratch venv):

- `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python -m unittest discover -s tests -v`:
  **82 passed** (48 from 0.3.0 + 34 new).
- Fail-before: the new modules did not exist before; the 2-fold self-stack
  regression fails with the pre-fix row logic (declared degree 3); the
  `hydrog` regression fails when every `_struct_conn` counts as covalent.
- `python tools/check_install.py`: PASS (installed nick, terminal-stack and
  ASU-graph smoke runs; the installed suite skips one jsonschema assertion
  because the clean venv has no jsonschema).
- NAtopo cross-validation on 17 crystals (1BNA, 1UHX, 1V3P, 5KEK, 5KEO, 6YL5,
  6YLB, 7SE0, 7U3R, 8EPG, 9Q5J, 9Y3B, 9ZXK, 9OSM and three local refinement
  PDBs), NAtopo venv running the notebook pipeline through
  `build_residue_quotient`:
  - **Nicks:** identical candidate orbits, relative images and distances
    (≤10⁻⁴ Å) to `find_backbone_join_candidates`.
  - **ASU graphs:** identical covalent orbits, backbone-join orbits (against
    every TopoStich join candidate), symmetry codes and operator numbering to
    NAtopo's `G_residue_quotient`; NAtopo's `periodicity()` on TopoStich's
    edges gives the same parts, ranks, bases and interpenetration counts at
    chain and residue resolution. This found the `hydrog` bug below.
- Native Systre 19.6.0 (jar SHA-256
  `0d272e98a1a21669bc67a809b95c014ba2a2fb39a6fd2147039201216ab0d48d`, size
  matches the GitHub release asset, which publishes no digest; system JDK 26.0.1).
  `tools/check_stitched_cgd.py` passed on all **85 exports** (17 crystals ×
  unstitched, nicks, stacks, both, all+O5′ fallback) with exact edge counts.
  Chain-graph dimensions agree with Systre's components in every export. All
  receipts and ASU JSON validated against their schemas.
- Default (unstitched) CGD, both CMM views and receipts are byte-identical to
  0.3.0 (except `tool_version`) for 1BNA, 5KEK, 6YL5, 7SE0, 9Y3B, 9ZXK and the
  82B and arch-solids PDBs.
- Both HTML views opened in headless Chrome for 7SE0 and 9Y3B and rendered
  as intended; their scripts pass `node --check`.

**Found by the checks.** (1) A base stacked on its own 2-fold mate (9Y3B,
6YL5, 6YLB) produced identical forward and reciprocal rows; TopoStich declared
degree 3 where Systre saw 2. Rows are now de-duplicated per node. (2) mmCIF
`hydrog` `_struct_conn` records at 2.26–2.27 Å (5KEK A8–D11, 9ZXK D8–E26)
were read as covalent links because the reader dropped `conn_type_id`. It is
now kept; only covalent types count.

**Not run.** ToposPro import; opening the CMM models in ChimeraX; the HTML
views in a browser other than headless Chrome.

**Remaining.** Reviewed confirm/remove decisions for ambiguous nicks (NAtopo's
decision log could feed this); special positions; live ToposPro; scientific
review of whether every chain-end stack should be stitched (stacks are not
exclusive, so an end may stack on two others).

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
