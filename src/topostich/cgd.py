"""Chain-aware O3-prime graphs in the named-source CRYSTAL CGD dialect.

This is a separate exporter: it does not change the existing CIF/PDB/mmCIF
converter. Explicit edges represent traces, not chemical bond assertions.
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


def _serialize(graph):
    lines = ['crystal', 'name '+graph['name'],
             'cell '+' '.join(_number(x) for x in graph['cell']),
             'group '+graph['space_group']]
    neighbors = {node['label']: [] for node in graph['nodes']}
    by_label = {node['label']: node for node in graph['nodes']}
    for left, right in graph['edges']:
        neighbors[left].append(right)
        neighbors[right].append(left)
    for node in graph['nodes']:
        label = node['label']
        lines.append(f"atom {label} {node['degree']} "+' '.join(map(_number,node['fractional'])))
        for neighbor in neighbors[label]:
            lines.append(f"edge {label} "+' '.join(map(_number,by_label[neighbor]['fractional'])))
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
               name=None, model=1, altloc=None, input_format='auto'):
    """Write a new CGD and optional ChimeraX views/JSON receipt; return receipt.

    Every output is self-contained and refused if it exists. cmm contains the
    observed asymmetric-unit trace. crystal_cmm contains one intact copy for
    each conventional-cell symmetry operation, with no inter-copy links.
    """
    source, output = Path(source), Path(output)
    targets = [Path(p) for p in (output,cmm,crystal_cmm,receipt) if p is not None]
    # Check before dict construction so identical companion paths cannot collapse.
    if len({p.resolve() for p in targets}) != len(targets):
        raise ConversionError('Output paths must be distinct')
    if any(p.resolve() == source.resolve() for p in targets):
        raise ConversionError('Input cannot be an output')
    if name is None:
        name = re.sub(r'[^A-Za-z0-9_.-]', '_', source.stem)+'_O3_chain_aware'
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
    degrees = Counter(i for edge in pairs for i in edge)
    nodes = [dict(label=labels[i],degree=degrees[i],segment=a.segment,
                  chain=a.chain,seq=a.seq,insertion=a.insertion,residue=a.residue,
                  atom=a.name,altloc=a.altloc,label_asym=a.label_asym,
                  label_seq=a.label_seq,label_comp=a.label_comp,label_atom=a.label_atom,
                  cartesian=list(a.xyz),fractional=fracs[i],description=a.description)
             for i,a in enumerate(atoms)]
    warnings = ['Unstitched finite strand components; no nick, junction, or symmetry-image joins.']
    if missing:
        warnings.append(f'{len(missing)} missing O3-prime selections skipped; retained nodes bridge within each segment.')
    if not pairs:
        warnings.append('No edges: all declared degrees are zero; Systre may report its nearest-node fallback.')
    graph = dict(schema_version=1,format='Systre CGD CRYSTAL',tool_version=__version__,
                 name=name,source_name=source.name,source_sha256=hashlib.sha256(raw).hexdigest(),
                 source_format=pdb.source_format,model=model,altloc=altloc,representative_atom="O3'",
                 cell=list(pdb.cell),space_group=group.xhm().replace(' ',''),
                 symmetry_operations=[op.triplet() for op in operations],
                 node_count=len(nodes),edge_count=len(pairs),directed_edge_rows=2*len(pairs),
                 segment_count=len(counters),nodes=nodes,
                 edges=[[labels[a],labels[b]] for a,b in pairs],missing_atoms=missing,
                 stitching=False,warnings=warnings)
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
    parser = argparse.ArgumentParser(description='O3-prime trace to Systre CRYSTAL CGD; keep chains separate.')
    parser.add_argument('input',type=Path,help='PDB or native PDBx/mmCIF')
    parser.add_argument('--out','-out',required=True,type=Path,help='new .cgd file')
    parser.add_argument('--cmm',type=Path,help='ChimeraX asymmetric-unit graph view')
    parser.add_argument('--crystal-cmm',type=Path,help='ChimeraX intact symmetry copies (no joins)')
    parser.add_argument('--receipt',type=Path,help='save node identities, edges and provenance as JSON')
    parser.add_argument('--name',help='CGD identifier; default derived from source basename')
    parser.add_argument('--input-format',choices=('auto','pdb','mmcif'),default='auto')
    parser.add_argument('--model',type=int,default=1)
    parser.add_argument('--altloc')
    parser.add_argument('--json',action='store_true',help='print the receipt as JSON')
    parser.add_argument('--version',action='version',version=f'TopoStich {__version__}')
    args=parser.parse_args(argv)
    try:
        result=export_cgd(args.input,args.out,cmm=args.cmm,crystal_cmm=args.crystal_cmm,
                          receipt=args.receipt,name=args.name,model=args.model,
                          altloc=args.altloc,input_format=args.input_format)
    except (ConversionError,OSError,UnicodeError) as error:
        parser.exit(2,f'topostich-cgd: {error}\n')
    if args.json:
        print(json.dumps(result,indent=2,allow_nan=False))
    else:
        print(f"Wrote {args.out}: {result['node_count']} nodes, {result['edge_count']} edges, "
              f"{result['segment_count']} separate strand segments")
        for warning in result['warnings']:
            print('Note: '+warning,file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
