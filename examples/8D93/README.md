# 8D93 example outputs

Generated with TopoStich 0.4.0 from PDB entry 8D93 (`8d93.cif`, SHA-256
`9f67eb6a7193d1ff06b68105e4ba8a066f4b02bef55787fc096404a99fa1d7bb`; the input
itself is not included). They were added on request as reference outputs.
They are not packaged, and they are regenerated rather than edited.

```sh
topostich-cgd 8d93.cif --stitch-nicks --stitch-terminal-stacks --name 8D93_stiched \
  --out 8D93_stiched.cgd --receipt 8D93_stiched.json \
  --cmm 8D93_stiched.cmm --crystal-cmm 8D93_stiched_symmetry.cmm \
  --asu-graphs 8D93_stiched_asu.json \
  --asu-chain-html 8D93_stiched_chains.html --asu-residue-html 8D93_stiched_residues.html
topostich-cgd 8d93.cif --name 8D93_unstiched \
  --out 8D93_unstiched.cgd --receipt 8D93_unstiched.json \
  --cmm 8D93_unstiched.cmm --crystal-cmm 8D93_unstiched_symmetry.cmm \
  --asu-graphs 8D93_unstiched_asu.json \
  --asu-chain-html 8D93_unstiched_chains.html --asu-residue-html 8D93_unstiched_residues.html
```

| File | Content |
|---|---|
| `*.cgd` | Systre CRYSTAL graph, one O3′ node per nucleotide |
| `*.json` | Receipt: node identities, edges, stitched joins and every candidate with its status |
| `*.cmm`, `*_symmetry.cmm` | ChimeraX marker views: asymmetric unit and intact symmetry copies |
| `*_asu.json` | ASU quotient graphs at chain and residue resolution, with periodicity |
| `*_chains.html`, `*_residues.html` | Interactive views of those graphs, open in a browser |

**Unstitched:** 42 nodes and 38 edges in four separate finite strands (A=12,
B=7, C=14, D=9).

**Stitched:** 42 nodes and 42 edges. The four added edges are:

- two backbone nick joins, both with code `3_555`: A:DA112 O3′ to D:DC113 P at
  2.93 Å, and B:DT114 to B:DA105 at 3.12 Å, which joins B to its 3-fold copies in a closed finite ring;
- two 5′-3′ chain-end stacks, both with code `9_555`: C:DC114 on C:DT101 at
  3.38 Å, and D:DA121 on A:DG101 at 3.91 Å.

Two further end stacks were not stitched because those ends are already
nick-joined. The ASU graph gives A+D a 1-periodic rod along [2/3, 1/3, 1/3],
C a rod along [0, 0, 1], and B finite.

Both files pass `tools/check_stitched_cgd.py` with Systre 19.6.0 (the edge
counts match the receipts).
