"""Standalone HTML views of the ASU quotient graphs (inline SVG, no libraries).

Layouts follow NAtopo's Bokeh views: ``asu_quotient.quotient_view`` (chains on
a circle, ligands outside) and ``residue_quotient.residue_quotient_view``
(residues around a circle by chain and sequence; rim arcs join sequence
neighbours, chords join everything else). Arrows point from the original ASU
node to the copy named by the symmetry code; ASU-internal edges are unlabelled.
"""
from __future__ import annotations

from collections import defaultdict
import html
import json
import math

from .asu_graphs import EDGE_LABELS

EDGE_COLORS = {'covalent': 'var(--covalent)', 'backbone_join': '#FF7A00', 'backbone_join_o5': '#FFB000',
               'terminal_stack': '#BD45E8'}
PALETTE = ('#2E86C1', '#E67E22', '#27AE60', '#8E44AD', '#C0392B', '#16A085', '#D4AC0D', '#7F8C8D',
           '#1F618D', '#AF601A')


def _add(p, q, s=1.):
    return (p[0] + s*q[0], p[1] + s*q[1])


def _scale(p, s):
    return (p[0]*s, p[1]*s)


def _norm(p):
    return math.hypot(*p)


def _loop(p0, size, angle):
    """Cubic teardrop self-loop: (path points, midpoint, tangent)."""
    d = (math.cos(angle), math.sin(angle))
    n = (-d[1], d[0])
    c0 = _add(p0, _add(_scale(d, 1.6*size), _scale(n, .9*size)))
    c1 = _add(p0, _add(_scale(d, 1.6*size), _scale(n, -.9*size)))
    mid = _add(_scale(p0, .25), _scale(_add(c0, c1), .375))
    return ('C', [p0, c0, c1, p0]), mid, (c1[0]-c0[0], c1[1]-c0[1])


def _quad(p0, control, p1):
    mid = _add(_add(_scale(p0, .25), _scale(control, .5)), _scale(p1, .25))
    return ('Q', [p0, control, p1]), mid, (p1[0]-p0[0], p1[1]-p0[1])


def _bend(p0, p1, bend):
    middle = _scale(_add(p0, p1), .5)
    normal = (-(p1[1]-p0[1]), p1[0]-p0[0])
    return _quad(p0, _add(middle, normal, bend), p1)


def _chain_layout(graph):
    nodes = [n['id'] for n in graph['nodes']]
    kind = {n['id']: n['kind'] for n in graph['nodes']}
    polymers = [n for n in nodes if kind[n] == 'polymer']
    others = [n for n in nodes if kind[n] != 'polymer']
    neighbours = defaultdict(set)
    for e in graph['edges']:
        if e['source'] != e['target']:
            neighbours[e['source']].add(e['target'])
            neighbours[e['target']].add(e['source'])
    pos = {}
    for i, n in enumerate(polymers):  # networkx.circular_layout
        theta = 2*math.pi*i/len(polymers)
        pos[n] = (math.cos(theta), math.sin(theta)) if len(polymers) > 1 else (0., 0.)
    polymer_set = set(polymers)
    attached = [n for n in others if neighbours[n] & polymer_set]
    loose = [n for n in others if n not in attached]
    slots = defaultdict(int)
    for i, n in enumerate(attached):
        partners = [pos[m] for m in neighbours[n] if m in polymer_set]
        anchor = (sum(p[0] for p in partners)/len(partners), sum(p[1] for p in partners)/len(partners))
        angle = math.atan2(anchor[1], anchor[0]) if _norm(anchor) > .15 else 2*math.pi*i/len(attached)
        slot = round(angle/.35)
        angle = slot*.35 + .45 + .12*slots[slot]
        radius = max(1.3, _norm(anchor) + .65)
        pos[n] = _scale((math.cos(angle), math.sin(angle)), radius + .18*slots[slot])
        slots[slot] += 1
    placed = polymers + attached
    center = (sum(pos[n][0] for n in placed)/len(placed), sum(pos[n][1] for n in placed)/len(placed)) \
        if placed else (0., 0.)
    bottom = min((p[1] for p in pos.values()), default=0.) - .45
    for i, n in enumerate(loose):
        in_row = min(14, len(loose) - (i // 14)*14)
        pos[n] = (((i % 14) - (in_row-1)/2)*.22, bottom - .22*(i // 14))
    return pos, center


def _residue_layout(graph):
    nodes = graph['nodes']
    polymers = [n['id'] for n in nodes if n['kind'] == 'polymer']
    chain = {n['id']: n['chain_node'] for n in nodes}
    chains = list(dict.fromkeys(chain[n] for n in polymers))
    gap = 2.5
    slots, angle_of, cursor = len(polymers) + gap*len(chains), {}, 0.
    for c in chains:
        for n in [n for n in polymers if chain[n] == c]:
            angle_of[n] = math.pi/2 - 2*math.pi*cursor/max(slots, 1)
            cursor += 1
        cursor += gap
    pos = {n: (math.cos(t), math.sin(t)) for n, t in angle_of.items()}
    partners = defaultdict(list)
    for e in graph['edges']:
        a, b = e['source'], e['target']
        if a in angle_of and b not in angle_of:
            partners[b].append(pos[a])
        elif b in angle_of and a not in angle_of:
            partners[a].append(pos[b])
    others = [n['id'] for n in nodes if n['id'] not in angle_of]
    used = defaultdict(int)
    for n in [n for n in others if partners[n]]:
        mean = (sum(p[0] for p in partners[n])/len(partners[n]), sum(p[1] for p in partners[n])/len(partners[n]))
        angle = math.atan2(mean[1], mean[0]) if _norm(mean) > 1e-6 else 0.
        slot = round(angle/.12)
        pos[n] = _scale((math.cos(slot*.12), math.sin(slot*.12)), 1.3 + .12*used[slot])
        used[slot] += 1
    loose = [n for n in others if not partners[n]]
    bottom = min([-1.] + [p[1] for p in pos.values()]) - .62
    for i, n in enumerate(loose):
        row = min(14, len(loose) - (i // 14)*14)
        pos[n] = (((i % 14) - (row-1)/2)*.16, bottom - .16*(i // 14))
    return pos, chains, angle_of


def _edge_geometry(graph, pos, resolution, angle_of=None, center=(0., 0.)):
    groups = defaultdict(list)
    for e in graph['edges']:
        groups[frozenset((e['source'], e['target']))].append(e)
    position = {n['id']: n.get('position_in_chain') for n in graph['nodes']}
    chain = {n['id']: n.get('chain_node', n['id']) for n in graph['nodes']}
    out = []
    for members in groups.values():
        members.sort(key=lambda e: (not e['intra_asu'], e['symmetry_code'], e['type']))
        for i, e in enumerate(members):
            a, b = e['source'], e['target']
            p0, p1 = pos[a], pos[b]
            spread = i - (len(members)-1)/2
            if a == b:
                outward = _add(p0, center, -1.)
                base = math.atan2(outward[1], outward[0]) if _norm(outward) > 1e-6 else math.pi/2
                if resolution == 'chain':
                    shape = _loop(p0, .28*(1 + .35*(i // 4)), base + spread*.55)
                else:
                    shape = _loop(p0, .07*(1 + .35*i), base)
            elif resolution == 'chain':
                shape = _bend(p0, p1, .18*spread)
            elif chain[a] == chain[b] and a in angle_of and b in angle_of \
                    and abs(position[a] - position[b]) == 1:
                shape = _bend(p0, p1, .22*spread)  # sequence neighbours: rim arcs
            else:
                middle = _scale(_add(p0, p1), .5)
                pull = .55 if (a in angle_of and b in angle_of) else .85
                shape = _quad(p0, _scale(middle, pull + .18*spread), p1)
            out.append((e, shape))
    return out


def _svg_point(p, s):
    return f'{p[0]*s:.2f},{-p[1]*s:.2f}'


def _info(lines):
    return html.escape('\n'.join(str(x) for x in lines if x), quote=True)


def _edge_lines(e, labels):
    lines = [EDGE_LABELS.get(e['type'], e['type']),
             f"From (original ASU): {labels[e['source']]}",
             f"To (copy): {labels[e['target']]} {e['symmetry_code']}" + (' (inside the ASU)' if e['intra_asu'] else ''),
             f"Operation: {e['operation']}   cell translation {e['lattice_translation']}"]
    if e['type'] == 'terminal_stack':
        lines.append(f"End polarity: {', '.join(e['polarities'])}" + (' (same-polarity stack)' if e['same_polarity'] else ''))
    if e['distance_min_A'] is not None:
        lo, hi = e['distance_min_A'], e['distance_max_A']
        lines.append(f'Distance: {lo:.2f} A' if abs(hi-lo) < 5e-3 else f'Distance: {lo:.2f}-{hi:.2f} A')
    if 'residue_pairs' in e:
        pairs = e['residue_pairs']
        lines.append(f"Residue pairs ({e['residue_edge_count']}): " + '; '.join(pairs[:8])
                     + (f' ... (+{len(pairs)-8})' if len(pairs) > 8 else ''))
    for text in e.get('contacts', [])[:6]:
        lines.append(text)
    return lines


def render_html(graph, *, title=None):
    """Return a self-contained HTML page for one ASU quotient graph."""
    resolution = graph['resolution']
    nodes = graph['nodes']
    labels = {n['id']: n['label'] for n in nodes}
    chain_ids = list(dict.fromkeys(n.get('chain_node', n['id']) for n in nodes))
    color = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(chain_ids)}
    angle_of, chains, center = None, [], (0., 0.)
    if resolution == 'chain':
        pos, center = _chain_layout(graph)
        scale = 260.
    else:
        pos, chains, angle_of = _residue_layout(graph)
        scale = 360.
    shapes = _edge_geometry(graph, pos, resolution, angle_of, center)
    loops = defaultdict(int)
    partners = defaultdict(set)
    for e in graph['edges']:
        if e['source'] == e['target']:
            loops[e['source']] += 1
        else:
            partners[e['source']].add(e['target'])
            partners[e['target']].add(e['source'])

    xs, ys = [], []
    edge_svg, arrow_svg, code_svg = [], [], []
    width = 3.0 if resolution == 'chain' else 2.0
    for e, ((cmd, pts), mid, tangent) in shapes:
        for p in pts + [mid]:
            xs.append(p[0]); ys.append(p[1])
        d = f'M{_svg_point(pts[0], scale)} {cmd}' + ' '.join(_svg_point(p, scale) for p in pts[1:])
        attrs = (f'data-type="{e["type"]}" data-sym="{0 if e["intra_asu"] else 1}" '
                 f'data-info="{_info(_edge_lines(e, labels))}"')
        dash = ' stroke-dasharray="7 4"' if e.get('same_polarity') else ''
        stroke = EDGE_COLORS.get(e['type'], '#888')
        edge_svg.append(f'<path class="edge" d="{d}" stroke="{stroke}" stroke-width="{width}"{dash} {attrs}/>')
        edge_svg.append(f'<path class="hit" d="{d}" {attrs}/>')
        if not e['intra_asu']:
            angle = math.degrees(math.atan2(-tangent[1], tangent[0]))
            size = 9 if resolution == 'chain' else 6
            mx, my = mid[0]*scale, -mid[1]*scale
            arrow_svg.append(f'<path class="arrow" data-type="{e["type"]}" data-sym="1" fill="{stroke}" '
                             f'd="M{size},0 L{-size*.7:.1f},{size*.6:.1f} L{-size*.7:.1f},{-size*.6:.1f} Z" '
                             f'transform="translate({mx:.2f},{my:.2f}) rotate({angle:.1f})"/>')
            tag = e['symmetry_code'] + (' ' + '/'.join(e['polarities']) if e.get('same_polarity') else '')
            code_svg.append(f'<text class="code" data-type="{e["type"]}" data-sym="1" x="{mx+6:.1f}" '
                            f'y="{my-4:.1f}" fill="{stroke}">{html.escape(tag)}</text>')
    node_svg, label_svg = [], []
    if resolution == 'chain':
        for n in nodes:
            p = pos[n['id']]
            xs.append(p[0]); ys.append(p[1])
            r = 17 if n['kind'] == 'polymer' else 9
            info = [f"Chain node {n['label']} (segment {n['segment']})", f"Residues: {n['residue_count']}",
                    f"Covalent links inside this chain copy: {n['internal_covalent_edges']}",
                    f"Self-loops (joins with own copies): {loops[n['id']]}"]
            node_svg.append(f'<circle class="node" cx="{p[0]*scale:.2f}" cy="{-p[1]*scale:.2f}" r="{r}" '
                            f'fill="{color[n["id"]]}" data-info="{_info(info)}"/>')
            label_svg.append(f'<text class="nodelabel" x="{p[0]*scale:.2f}" y="{-p[1]*scale:.2f}">'
                             f'{html.escape(n["label"])}</text>')
    else:
        r = max(2.5, min(6.5, 260/max(len(nodes), 1)))
        for n in nodes:
            p = pos[n['id']]
            xs.append(p[0]); ys.append(p[1])
            info = [f"Residue {n['label']}", f"Chain {n['chain_label']} (segment {n['segment']})",
                    f"Distinct partner residues: {len(partners[n['id']])}",
                    f"Self-loops (joins with own copies): {loops[n['id']]}"]
            on_rim = n['id'] in angle_of
            node_svg.append(f'<circle class="node" cx="{p[0]*scale:.2f}" cy="{-p[1]*scale:.2f}" '
                            f'r="{r if on_rim else r*1.25:.2f}" fill="{color[n["chain_node"]]}" data-info="{_info(info)}"/>')
            angle = math.atan2(p[1], p[0]) if on_rim else 0.
            flip = on_rim and math.cos(angle) < 0
            q = _scale(p, 1.06) if on_rim else _add(p, (0, -.06))
            rot = -math.degrees(angle + (math.pi if flip else 0.)) if on_rim else 0.
            anchor = 'end' if flip else ('start' if on_rim else 'middle')
            text = f"{n['residue']}{n['seq']}{n['insertion']}" if on_rim else n['label']
            label_svg.append(f'<text class="reslabel" text-anchor="{anchor}" transform="translate({q[0]*scale:.2f},'
                             f'{-q[1]*scale:.2f}) rotate({rot:.1f})">{html.escape(text)}</text>')
        for c in chains:
            members = [pos[n['id']] for n in nodes if n['chain_node'] == c and n['id'] in angle_of]
            mean = (sum(p[0] for p in members)/len(members), sum(p[1] for p in members)/len(members))
            point = _scale(mean, 1.42/_norm(mean)) if _norm(mean) > 1e-6 else (0., 1.42)
            xs.append(point[0]*1.08); ys.append(point[1]*1.08)
            name = next(n['chain_label'] for n in nodes if n['chain_node'] == c)
            label_svg.append(f'<text class="chainname" x="{point[0]*scale:.2f}" y="{-point[1]*scale:.2f}" '
                             f'fill="{color[c]}">{html.escape(name)}</text>')
    pad = .32
    x0, x1 = (min(xs, default=-1.) - pad)*scale, (max(xs, default=1.) + pad)*scale
    y0, y1 = (-max(ys, default=1.) - pad)*scale, (-min(ys, default=-1.) + pad)*scale
    view = f'{x0:.1f} {y0:.1f} {x1-x0:.1f} {y1-y0:.1f}'
    kinds = [k for k in EDGE_LABELS if any(e['type'] == k for e in graph['edges'])]
    toggles = ''.join(f'<label><input type="checkbox" data-kind="{k}" checked> '
                      f'<span class="swatch" style="background:{EDGE_COLORS[k]}"></span>{html.escape(EDGE_LABELS[k])}'
                      f' ({graph["edge_type_counts"][k]})</label>' for k in kinds)
    if title is None:
        what = 'chain' if resolution == 'chain' else 'residue'
        title = f"{graph['name']}: one ASU at {what} resolution ({graph['space_group']})"
    unit = 'chains/TER segments' if resolution == 'chain' else 'residues'
    summary = ''.join(f'<li>{html.escape(line)}</li>' for line in graph['summary']) or '<li>No nodes.</li>'
    data = json.dumps(dict(view=view))
    return f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(graph['name'])} ASU {resolution} graph</title>
<style>
:root {{ --bg:#ffffff; --fg:#1d2730; --muted:#5b6670; --panel:#f3f5f7; --line:#d5dbe0; --covalent:#555555; --node-stroke:#25313b; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg:#14181c; --fg:#e4e9ee; --muted:#9aa5ae; --panel:#1e2429; --line:#323a42; --covalent:#b4bcc4; --node-stroke:#e4e9ee; }} }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--fg); font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif; }}
header {{ padding:14px 16px 6px; }}
h1 {{ font-size:17px; margin:0 0 6px; }}
p, ul {{ margin:4px 0; color:var(--muted); max-width:1100px; }}
ul {{ padding-left:20px; color:var(--fg); }}
.controls {{ display:flex; flex-wrap:wrap; gap:6px 16px; padding:8px 16px; border-top:1px solid var(--line); border-bottom:1px solid var(--line); background:var(--panel); }}
.controls label {{ display:inline-flex; align-items:center; gap:5px; cursor:pointer; }}
.swatch {{ width:14px; height:4px; border-radius:2px; display:inline-block; }}
button {{ font:inherit; padding:2px 10px; border:1px solid var(--line); border-radius:4px; background:var(--bg); color:var(--fg); cursor:pointer; }}
.wrap {{ display:grid; grid-template-columns:minmax(0,1fr) 320px; min-height:70vh; }}
@media (max-width: 800px) {{ .wrap {{ grid-template-columns:1fr; }} }}
svg {{ width:100%; height:78vh; display:block; cursor:grab; touch-action:none; }}
svg.dragging {{ cursor:grabbing; }}
#info {{ border-left:1px solid var(--line); padding:12px 14px; white-space:pre-wrap; font-size:13px; background:var(--panel); overflow:auto; max-height:78vh; }}
.edge {{ fill:none; stroke-linecap:round; opacity:.88; vector-effect:non-scaling-stroke; }}
.hit {{ fill:none; stroke:transparent; stroke-width:12; vector-effect:non-scaling-stroke; cursor:pointer; }}
.node {{ stroke:var(--node-stroke); stroke-width:1; vector-effect:non-scaling-stroke; cursor:pointer; }}
.nodelabel {{ fill:#fff; font-weight:700; font-size:12px; text-anchor:middle; dominant-baseline:central; pointer-events:none; }}
.reslabel {{ fill:var(--muted); font-size:8px; dominant-baseline:central; pointer-events:none; }}
.chainname {{ font-weight:700; font-size:15px; text-anchor:middle; dominant-baseline:central; pointer-events:none; }}
.code {{ font-size:{10 if resolution == 'chain' else 8}px; paint-order:stroke; stroke:var(--bg); stroke-width:3px; pointer-events:none; }}
.hidden {{ display:none; }}
.nolabels .reslabel, .nolabels .code {{ display:none; }}
</style></head><body>
<header><h1>{html.escape(title)}</h1>
<p>Nodes are the {unit} of the <b>original asymmetric unit</b>. Every edge is one symmetry orbit of covalent bonds or stitched joins: X &rarr; Y labelled g means X of the original ASU is joined to the copy of Y made by <b>symmetry code</b> g (op_klm; 555 = no cell shift; unlabelled = inside the ASU). Self-loops join a node to another copy of itself. Composing the codes around any cycle gives the operation the real lattice path ends on; a non-identity result means the graph continues through the crystal. Dashed purple edges are same-polarity (5&prime;-5&prime; or 3&prime;-3&prime;) chain-end stacks. Hover or tap an edge or node for details; drag to pan, scroll to zoom.</p>
<p><b>Periodicity</b> (all edges):</p><ul>{summary}</ul></header>
<div class="controls">{toggles}<label><input type="checkbox" id="symonly"> only symmetry joins</label><label><input type="checkbox" id="labels" checked> labels</label><button id="reset" type="button">Reset view</button></div>
<div class="wrap"><svg id="graph" viewBox="{view}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="{html.escape(title)}">
<g id="edges">{''.join(edge_svg)}</g><g id="arrows">{''.join(arrow_svg)}</g><g id="nodes">{''.join(node_svg)}</g><g id="labels-g">{''.join(label_svg)}{''.join(code_svg)}</g>
</svg><div id="info">Hover or tap an edge or node.</div></div>
<script>
(function () {{
  const cfg = {data};
  const svg = document.getElementById('graph'), info = document.getElementById('info');
  const initial = cfg.view.split(' ').map(Number);
  let box = initial.slice();
  const apply = () => svg.setAttribute('viewBox', box.join(' '));
  function filter() {{
    const on = new Set([...document.querySelectorAll('input[data-kind]')].filter(i => i.checked).map(i => i.dataset.kind));
    const symOnly = document.getElementById('symonly').checked;
    svg.querySelectorAll('[data-type]').forEach(el => {{
      el.classList.toggle('hidden', !on.has(el.dataset.type) || (symOnly && el.dataset.sym === '0'));
    }});
    svg.classList.toggle('nolabels', !document.getElementById('labels').checked);
  }}
  document.querySelectorAll('.controls input').forEach(i => i.addEventListener('change', filter));
  document.getElementById('reset').addEventListener('click', () => {{ box = initial.slice(); apply(); }});
  svg.querySelectorAll('[data-info]').forEach(el => {{
    const show = () => {{ info.textContent = el.dataset.info; }};
    el.addEventListener('mouseenter', show); el.addEventListener('click', show);
  }});
  svg.addEventListener('wheel', ev => {{
    ev.preventDefault();
    const r = svg.getBoundingClientRect(), f = ev.deltaY > 0 ? 1.15 : 1 / 1.15;
    const s = Math.max(box[2] / r.width, box[3] / r.height);
    const cx = box[0] + (ev.clientX - r.left - (r.width - box[2] / s) / 2) * s;
    const cy = box[1] + (ev.clientY - r.top - (r.height - box[3] / s) / 2) * s;
    box = [cx - (cx - box[0]) * f, cy - (cy - box[1]) * f, box[2] * f, box[3] * f]; apply();
  }}, {{ passive: false }});
  let drag = null;
  svg.addEventListener('pointerdown', ev => {{ drag = [ev.clientX, ev.clientY, box[0], box[1]]; svg.classList.add('dragging'); }});
  window.addEventListener('pointerup', () => {{ drag = null; svg.classList.remove('dragging'); }});
  window.addEventListener('pointermove', ev => {{
    if (!drag) return;
    const r = svg.getBoundingClientRect(), s = Math.max(box[2] / r.width, box[3] / r.height);
    box[0] = drag[2] - (ev.clientX - drag[0]) * s; box[1] = drag[3] - (ev.clientY - drag[1]) * s; apply();
  }});
}})();
</script></body></html>
'''
