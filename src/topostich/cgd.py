"""Chain-aware O3-prime graphs in the named-source CRYSTAL CGD dialect.

This is a separate exporter: it does not change the existing CIF/PDB/mmCIF
converter. Explicit edges represent traces, not chemical bond assertions.
Optional nick (``stitch_nicks``) and junction (``stitch_junctions``) stitching
add symmetry-aware joins found by the NAtopo-derived rules in ``stitching.py``;
the default is unchanged.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import sys

import gemmi

from . import __version__
from .converter import ConversionError, _fractional_transform, _select, read_structure
from .graph_viewer import marker_xml
from .stitching import JUNCTION_SCOPES, StitchContext, find_backbone_joins, find_junction_stacks


def _space_group(pdb):
    group = gemmi.find_spacegroup_by_name(pdb.space_group)
    if group is None:
        raise ConversionError(f"Unknown space group: {pdb.space_group}")
    if not gemmi.UnitCell(*pdb.cell).is_compatible_with_spacegroup(group):
        raise ConversionError("Cell is incompatible with the named space-group setting")
    operations = list(group.operations())
    if pdb.symops:
        try:
            supplied = {gemmi.Op(op).wrap().triplet() for op in pdb.symops.values()}
        except (RuntimeError, ValueError) as error:
            raise ConversionError(f"Invalid input symmetry operation: {error}") from error
        expected = {op.wrap().triplet() for op in operations}
        if supplied != expected:
            raise ConversionError("Input symmetry operations disagree with the complete named space-group setting")
    return group, operations


def _check_positions(fracs, operations):
    # Systre compares points modulo its primitive lattice at 0.001 fractional
    # precision. This conservative conventional-cell guard (0.003 in each
    # coordinate) also catches near-coincidences after centering transformations.
    # Special positions need stabilizer-aware edge degrees: defer them explicitly.
    tolerance = .003
    identity = gemmi.Op('x,y,z').triplet()
    for i, point in enumerate(fracs):
        for op in operations:
            image = op.apply_to_xyz(point)
            for j in range(i + 1):
                if i == j and op.triplet() == identity:
                    continue
                delta = [a-b for a,b in zip(image,fracs[j])]
                if max(abs(d-round(d)) for d in delta) <= tolerance:
                    kind = 'special position' if i == j else 'symmetry-equivalent or near-coincident nodes'
                    raise ConversionError(f"CGD cannot preserve separate trace identities at {kind} "
                                          f"(nodes {j+1}, {i+1}; operation {op.triplet()})")


def _number(value):
    return format(value, '.12g')


def _neighbor_rows(nodes, edges, stitched=()):
    """CGD neighbor-image rows per node label, without repeated rows.

    Periodic joins use the inverse affine image for the reciprocal row, not a
    copied or negated translation. A join from a node to its own image under an
    involution (e.g. a base stacked on its 2-fold mate) has identical forward
    and reciprocal rows: it is one neighbor, written once and counted once.
    """
    neighbors = {node['label']: [] for node in nodes}
    seen = {node['label']: set() for node in nodes}
    by_label = {node['label']: node for node in nodes}

    def add(label, image):
        key = tuple(round(x, 9) for x in image)
        if key not in seen[label]:
            seen[label].add(key)
            neighbors[label].append(image)
    for left, right in edges:
        add(left, by_label[right]['fractional'])
        add(right, by_label[left]['fractional'])
    for edge in stitched:
        add(edge['source'], edge['target_image_fractional'])
        add(edge['target'], edge['reverse_image_fractional'])
    return neighbors


def _serialize(graph):
    lines = ['crystal', 'name '+graph['name'],
             'cell '+' '.join(_number(x) for x in graph['cell']),
             'group '+graph['space_group']]
    neighbors = _neighbor_rows(graph['nodes'], graph['edges'], graph.get('stitched_edges', ()))
    for node in graph['nodes']:
        label = node['label']
        lines.append(f"atom {label} {node['degree']} "+' '.join(map(_number,node['fractional'])))
        for image in neighbors[label]:
            lines.append(f"edge {label} "+' '.join(map(_number,image)))
    lines.append('end')
    return '\n'.join(lines)+'\n'


def _write_new_files(payloads):
    paths = [path.resolve() for path in payloads]
    if len(paths) != len(set(paths)):
        raise ConversionError('Output paths must be distinct')
    for path in payloads:
        if path.exists() or path.is_symlink():
            raise ConversionError(f'Refusing to overwrite output: {path}')
        if not path.parent.is_dir():
            raise ConversionError(f'Output parent does not exist: {path.parent}')
    created = []
    try:
        for path, content in payloads.items():
            with path.open('xb') as stream:
                created.append(path)
                stream.write(content)
    except OSError as error:
        for path in created:
            path.unlink()
        raise ConversionError(str(error)) from error


def export_cgd(source, output, *, cmm=None, crystal_cmm=None, receipt=None,
               name=None, model=1, altloc=None, input_format='auto',
               stitch_nicks=False, join_cutoff=3.5, o5_cutoff=6.0, include_o5_fallback=False,
               stitch_junctions=False, junction_scope='face'):
    """Write a new CGD and optional ChimeraX views/JSON receipt; return receipt.

    Every output is self-contained and refused if it exists. cmm contains the
    observed asymmetric-unit trace. crystal_cmm contains one intact copy for
    each conventional-cell symmetry operation.

    With stitch_nicks, unambiguous P...O3' backbone nick joins between segment
    termini (NAtopo rule, any symmetry image) become periodic edges. With
    stitch_junctions, base stacks without verified backbone continuation
    (NAtopo junction rule; scope 'face' or 'all') do too. Either option makes
    the receipt schema 2. Other candidates are reported, never silently used.
    """
    stitching = bool(stitch_nicks or stitch_junctions)
    source, output = Path(source), Path(output)
    targets = [Path(p) for p in (output,cmm,crystal_cmm,receipt) if p is not None]
    # Check before dict construction so identical companion paths cannot collapse.
    if len({p.resolve() for p in targets}) != len(targets):
        raise ConversionError('Output paths must be distinct')
    if any(p.resolve() == source.resolve() for p in targets):
        raise ConversionError('Input cannot be an output')
    if name is None:
        suffix = {(False, False): '_O3_chain_aware', (True, False): '_O3_nick_stitched',
                  (False, True): '_O3_junction_stitched', (True, True): '_O3_stitched'}
        name = re.sub(r'[^A-Za-z0-9_.-]', '_', source.stem)+suffix[bool(stitch_nicks), bool(stitch_junctions)]
    if not re.fullmatch(r'[A-Za-z0-9_][A-Za-z0-9_.-]*',name):
        raise ConversionError('Graph name must be a single ASCII identifier')
    raw = source.read_bytes()
    pdb = read_structure(raw.decode('utf-8'),model=model,input_format=input_format)
    atoms, pairs, _, missing = _select(pdb,["O3'"],altloc)
    transform = _fractional_transform(pdb)
    group, operations = _space_group(pdb)
    fracs = [[sum(row[k]*atom.xyz[k] for k in range(3))+row[3] for row in transform]
             for atom in atoms]
    _check_positions(fracs,operations)
    counters, labels = Counter(), []
    for atom in atoms:
        counters[atom.segment] += 1
        labels.append(f'S{atom.segment}_N{counters[atom.segment]}')
    stitched, candidates, junctions, nick_summary, junction_summary = [], [], [], None, None
    if stitch_junctions and junction_scope not in JUNCTION_SCOPES:
        raise ConversionError(f'junction_scope must be one of {", ".join(JUNCTION_SCOPES)}')
    if stitching:
        context = StitchContext(pdb, altloc=altloc, node_atoms=atoms, node_labels=labels,
                                transform=transform, operations=operations, trace_pairs=pairs)
        if stitch_nicks:
            edges, candidates, nick_summary = find_backbone_joins(
                context, join_cutoff=join_cutoff, o5_cutoff=o5_cutoff,
                include_o5_fallback=include_o5_fallback)
            stitched += edges
        if stitch_junctions:
            edges, junctions, junction_summary = find_junction_stacks(context, scope=junction_scope)
            stitched += edges
    rows = _neighbor_rows([dict(label=label, fractional=frac) for label, frac in zip(labels, fracs)],
                          [[labels[a], labels[b]] for a, b in pairs], stitched)
    nodes = [dict(label=labels[i],degree=len(rows[labels[i]]),segment=a.segment,
                  chain=a.chain,seq=a.seq,insertion=a.insertion,residue=a.residue,
                  atom=a.name,altloc=a.altloc,label_asym=a.label_asym,
                  label_seq=a.label_seq,label_comp=a.label_comp,label_atom=a.label_atom,
                  cartesian=list(a.xyz),fractional=fracs[i],description=a.description)
             for i,a in enumerate(atoms)]
    if stitching:
        kinds = Counter(e['kind'] for e in stitched)
        warnings = []
        if stitch_nicks:
            counts = Counter(c['status'] for c in candidates)
            nick_count = len(stitched) - kinds['junction_stack']
            warnings.append(f'Nick stitching: {nick_count} backbone joins stitched.')
            for status, text in (('ambiguous', 'ambiguous nick candidates (several partners) were not stitched'),
                                 ('segment_internal', 'nick candidates at segment-internal breaks were not stitched'),
                                 ('o5_fallback_disabled', "O5' fallback nick candidates reported only (use include_o5_fallback)"),
                                 ('no_graph_node', 'nick candidates without an O3-prime graph node were not stitched'),
                                 ('duplicate_edge', 'nick candidates duplicating an existing edge were not stitched')):
                if counts[status]:
                    warnings.append(f'{counts[status]} {text}.')
        if stitch_junctions:
            counts = Counter(c['status'] for c in junctions)
            warnings.append(f"Junction stitching ({junction_scope} scope): {kinds['junction_stack']} base-stack junctions stitched.")
            for status, text in (('not_face_junction', 'stacks with backbone on both faces reported only (junction_scope all stitches them)'),
                                 ('no_graph_node', 'junction stacks without an O3-prime graph node were not stitched'),
                                 ('duplicate_edge', 'junction stacks duplicating an existing edge were not stitched')):
                if counts[status]:
                    warnings.append(f'{counts[status]} {text}.')
            issues = junction_summary['plane_diagnostics']
            if issues:
                warnings.append(f'{len(issues)} nucleotide residues lack a usable base-ring plane; see stitch_parameters.')
        if not stitch_junctions:
            warnings.append('Junctions (non-backbone base stacks) are not stitched.')
        if not stitch_nicks:
            warnings.append('Backbone nicks are not stitched.')
    else:
        warnings = ['Unstitched finite strand components; no nick, junction, or symmetry-image joins.']
    if missing:
        warnings.append(f'{len(missing)} missing O3-prime selections skipped; retained nodes bridge within each segment.')
    if not pairs and not stitched:
        warnings.append('No edges: all declared degrees are zero; Systre may report its nearest-node fallback.')
    graph = dict(schema_version=2 if stitching else 1,format='Systre CGD CRYSTAL',tool_version=__version__,
                 name=name,source_name=source.name,source_sha256=hashlib.sha256(raw).hexdigest(),
                 source_format=pdb.source_format,model=model,altloc=altloc,representative_atom="O3'",
                 cell=list(pdb.cell),space_group=group.xhm().replace(' ',''),
                 symmetry_operations=[op.triplet() for op in operations],
                 node_count=len(nodes),edge_count=len(pairs)+len(stitched),
                 directed_edge_rows=2*(len(pairs)+len(stitched)),
                 segment_count=len(counters),nodes=nodes,
                 edges=[[labels[a],labels[b]] for a,b in pairs],missing_atoms=missing,
                 stitching=stitching,warnings=warnings)
    if stitching:
        graph.update(trace_edge_count=len(pairs),stitched_edge_count=len(stitched),
                     stitched_edge_counts=dict(sorted(Counter(e['kind'] for e in stitched).items())),
                     stitched_edges=stitched,join_candidates=candidates,junction_candidates=junctions,
                     stitch_parameters=dict(nicks=nick_summary,junctions=junction_summary))
    payloads = {output:_serialize(graph).encode()}
    if cmm is not None:
        payloads[Path(cmm)] = marker_xml(graph).encode()
    if crystal_cmm is not None:
        payloads[Path(crystal_cmm)] = marker_xml(graph,crystal=True).encode()
    graph['artifacts'] = [dict(name=p.name,sha256=hashlib.sha256(data).hexdigest())
                          for p,data in payloads.items()]
    if receipt is not None:
        payloads[Path(receipt)] = (json.dumps(graph,indent=2,allow_nan=False)+'\n').encode()
    _write_new_files(payloads)
    return graph


def main(argv=None):
    parser = argparse.ArgumentParser(description='O3-prime trace to Systre CRYSTAL CGD; chains stay separate unless --stitch-nicks/--stitch-junctions.')
    parser.add_argument('input',type=Path,help='PDB or native PDBx/mmCIF')
    parser.add_argument('--out','-out',required=True,type=Path,help='new .cgd file')
    parser.add_argument('--cmm',type=Path,help='ChimeraX asymmetric-unit graph view')
    parser.add_argument('--crystal-cmm',type=Path,help='ChimeraX intact symmetry copies (plus stitched joins)')
    parser.add_argument('--receipt',type=Path,help='save node identities, edges and provenance as JSON')
    parser.add_argument('--name',help='CGD identifier; default derived from source basename')
    parser.add_argument('--input-format',choices=('auto','pdb','mmcif'),default='auto')
    parser.add_argument('--model',type=int,default=1)
    parser.add_argument('--altloc')
    parser.add_argument('--stitch-nicks',action='store_true',
                        help="join unambiguous loose P...O3' backbone nicks across chains/symmetry images (NAtopo rule)")
    parser.add_argument('--join-cutoff',type=float,default=3.5,help="P...O3' cutoff in Angstrom (default 3.5)")
    parser.add_argument('--o5-cutoff',type=float,default=6.0,help="missing-phosphate O5'...O3' cutoff (default 6.0)")
    parser.add_argument('--include-o5-fallback',action='store_true',
                        help="also stitch unambiguous O5' fallback joins (reported only by default)")
    parser.add_argument('--stitch-junctions',action='store_true',
                        help='join base stacks without backbone continuation across chains/images (NAtopo junction rule)')
    parser.add_argument('--junction-scope',choices=JUNCTION_SCOPES,default='face',
                        help="face (default): only stacks on a base face with no backbone neighbor; all: every non-backbone stack")
    parser.add_argument('--json',action='store_true',help='print the receipt as JSON')
    parser.add_argument('--version',action='version',version=f'TopoStich {__version__}')
    args=parser.parse_args(argv)
    try:
        result=export_cgd(args.input,args.out,cmm=args.cmm,crystal_cmm=args.crystal_cmm,
                          receipt=args.receipt,name=args.name,model=args.model,
                          altloc=args.altloc,input_format=args.input_format,
                          stitch_nicks=args.stitch_nicks,join_cutoff=args.join_cutoff,
                          o5_cutoff=args.o5_cutoff,include_o5_fallback=args.include_o5_fallback,
                          stitch_junctions=args.stitch_junctions,junction_scope=args.junction_scope)
    except (ConversionError,OSError,UnicodeError) as error:
        parser.exit(2,f'topostich-cgd: {error}\n')
    if args.json:
        print(json.dumps(result,indent=2,allow_nan=False))
    elif result['stitching']:
        kinds = result['stitched_edge_counts']
        nicks = result['stitched_edge_count'] - kinds.get('junction_stack', 0)
        parts = ([f"{nicks} stitched nick joins"] if args.stitch_nicks else []) + \
                ([f"{kinds.get('junction_stack', 0)} stitched junction stacks"] if args.stitch_junctions else [])
        print(f"Wrote {args.out}: {result['node_count']} nodes, {result['edge_count']} edges "
              f"({result['trace_edge_count']} trace, {', '.join(parts)}), "
              f"{result['segment_count']} strand segments")
    else:
        print(f"Wrote {args.out}: {result['node_count']} nodes, {result['edge_count']} edges, "
              f"{result['segment_count']} separate strand segments")
    if not args.json:
        for warning in result['warnings']:
            print('Note: '+warning,file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
