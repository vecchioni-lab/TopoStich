"""ChimeraX native marker/link views from the same explicit CGD graph.

No distance-based bonding is performed. Every chain/TER segment and symmetry
copy is a separate model so ChimeraX's Models panel can toggle it independently.
"""
from collections import defaultdict
import xml.etree.ElementTree as ET

import gemmi


PALETTE = ((.16,.65,.85),(.95,.55,.18),(.50,.78,.35),(.72,.46,.85),(.93,.35,.45),(.28,.77,.70))


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
                ET.SubElement(model,'marker',id=ids[node['label']],
                              **dict(zip(('x','y','z'),(format(x,'.12g') for x in xyz))),
                              **rgb,radius='0.85' if node['degree']<=1 else '0.5',note=note)
            for left,right in graph['edges']:
                if left in ids and right in ids:
                    ET.SubElement(model,'link',id1=ids[left],id2=ids[right],**rgb,radius='0.22')
    ET.indent(root)
    return ET.tostring(root,encoding='unicode')+'\n'
