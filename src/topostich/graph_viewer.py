"""ChimeraX native marker/link views from the same explicit CGD graph.

No distance-based bonding is performed. Every chain/TER segment and symmetry
copy is a separate model so ChimeraX's Models panel can toggle it independently.
Stitched nick joins and chain-end stacks cross segments and images, so each
family gets its own model with explicit endpoint markers; every stitched edge
is drawn in each view.
"""
from collections import defaultdict
import xml.etree.ElementTree as ET

import gemmi


PALETTE = ((.16,.65,.85),(.95,.55,.18),(.50,.78,.35),(.72,.46,.85),(.93,.35,.45),(.28,.77,.70))
JOIN_FAMILIES = (('nick joins', ('backbone_join', 'backbone_join_o5', 'deposited_link'), dict(r='1', g='0.84', b='0')),
                 ('terminal stacks', ('terminal_stack',), dict(r='0.9', g='0.2', b='0.75')))


def _xyz(values):
    return dict(zip(('x','y','z'),(format(x,'.12g') for x in values)))


def _join_model(root, graph, ops, crystal, family, kinds, rgb):
    edges = [e for e in graph.get('stitched_edges') or () if e['kind'] in kinds]
    if not edges:
        return
    name = f"{graph['name']} {family}"
    if crystal:
        name += ' (all symmetry copies)'
    model = ET.SubElement(root,'marker_set',name=name)
    cell = gemmi.UnitCell(*graph['cell'])
    by_label = {node['label']: node for node in graph['nodes']}
    serial = 0
    for op_index, triplet in enumerate(ops,1):
        op = gemmi.Op(triplet)
        for edge in edges:
            source = by_label[edge['source']]
            if crystal:
                ends = [list(cell.orthogonalize(gemmi.Fractional(*op.apply_to_xyz(f))))
                        for f in (source['fractional'], edge['target_image_fractional'])]
            else:
                ends = [source['cartesian'], edge['target_image_cartesian']]
            notes = [f"{edge['source']} {edge['kind']} {edge['source_description']}",
                     f"{edge['target']} image {edge['image_operation']} {edge['target_description']}"]
            if crystal:
                notes = [n+f' symmetry {op_index}: {triplet}' for n in notes]
            ids = []
            for xyz, note in zip(ends, notes):
                serial += 1
                ids.append(str(serial))
                ET.SubElement(model,'marker',id=ids[-1],**_xyz(xyz),**rgb,radius='0.6',note=note)
            ET.SubElement(model,'link',id1=ids[0],id2=ids[1],**rgb,radius='0.3',
                          note=f"{edge['kind']} {edge['distance_A']:.3f} A")


def marker_xml(graph, *, crystal=False):
    root = ET.Element('marker_sets')
    segments = defaultdict(list)
    for node in graph['nodes']:
        segments[node['segment']].append(node)
    ops = graph['symmetry_operations'] if crystal else ['x,y,z']
    cell = gemmi.UnitCell(*graph['cell'])
    for op_index, triplet in enumerate(ops,1):
        op = gemmi.Op(triplet)
        for color_index,(segment,nodes) in enumerate(segments.items()):
            rgb = dict(zip(('r','g','b'),map(str,PALETTE[color_index % len(PALETTE)])))
            chain = nodes[0]['chain'] or '(blank)'
            name = f"{graph['name']} chain {chain} segment {segment}"
            if crystal:
                name += f' symmetry {op_index}: {triplet}'
            model = ET.SubElement(root,'marker_set',name=name)
            ids = {n['label']:str(i) for i,n in enumerate(nodes,1)}
            for node in nodes:
                xyz = node['cartesian']
                if crystal:
                    frac = op.apply_to_xyz(node['fractional'])
                    xyz = list(cell.orthogonalize(gemmi.Fractional(*frac)))
                note = node['label']+' '+node['description']+f" degree={node['degree']}"
                ET.SubElement(model,'marker',id=ids[node['label']],**_xyz(xyz),
                              **rgb,radius='0.85' if node['degree']<=1 else '0.5',note=note)
            for left,right in graph['edges']:
                if left in ids and right in ids:
                    ET.SubElement(model,'link',id1=ids[left],id2=ids[right],**rgb,radius='0.22')
    for family, kinds, rgb in JOIN_FAMILIES:
        _join_model(root, graph, ops, crystal, family, kinds, rgb)
    ET.indent(root)
    return ET.tostring(root,encoding='unicode')+'\n'
