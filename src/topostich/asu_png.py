"""PNG images of the ASU quotient graphs (optional; needs matplotlib).

Same layouts, edge geometry, colors and arrowheads as the HTML views in
``asu_view.py``, drawn with matplotlib's Agg backend for papers and slides.
Install with ``pip install "topostich[png]"``. Gemmi stays the only required
dependency; without matplotlib the PNG options fail before any file is written.
"""
from __future__ import annotations

import io
import math

from .asu_graphs import EDGE_LABELS
from .asu_view import (ARROW_SIZE, EDGE_COLORS, PALETTE, _add, _chain_layout, _edge_geometry, _norm,
                       _residue_layout, _scale)
from .converter import ConversionError

COLORS = dict(EDGE_COLORS, covalent='#555555')


def _pyplot():
    try:
        import matplotlib
        matplotlib.use('Agg', force=True)
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.path import Path as MPath
        from matplotlib.patches import PathPatch
    except ImportError as error:
        raise ConversionError('PNG output needs matplotlib: pip install "topostich[png]"') from error
    return plt, Line2D, MPath, PathPatch


def render_png(graph, *, title=None, dpi=200):
    """Return PNG bytes for one ASU quotient graph."""
    plt, Line2D, MPath, PathPatch = _pyplot()
    resolution = graph['resolution']
    nodes = graph['nodes']
    chain_ids = list(dict.fromkeys(n.get('chain_node', n['id']) for n in nodes))
    color = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(chain_ids)}
    angle_of, chains, center = None, [], (0., 0.)
    if resolution == 'chain':
        pos, center = _chain_layout(graph)
    else:
        pos, chains, angle_of = _residue_layout(graph)
    shapes = _edge_geometry(graph, pos, resolution, angle_of, center)
    summary = graph['summary'] or ['No nodes.']
    fig = plt.figure(figsize=(11, 11 + 0.22*len(summary)), dpi=dpi, facecolor='white')
    top = 1 - (0.6 + 0.22*len(summary)) / (11 + 0.22*len(summary))
    ax = fig.add_axes((0.03, 0.03, 0.94, top - 0.04))
    ax.set_aspect('equal')
    ax.axis('off')
    xs, ys = [], []
    width = 2.4 if resolution == 'chain' else 1.6
    arrow = ARROW_SIZE[resolution]
    codes_font = 9 if resolution == 'chain' else 7
    curve = {'Q': MPath.CURVE3, 'C': MPath.CURVE4}
    for e, ((cmd, pts), mid, tangent) in shapes:
        xs.extend(p[0] for p in pts + [mid]); ys.extend(p[1] for p in pts + [mid])
        stroke = COLORS.get(e['type'], '#888888')
        path = MPath(pts, [MPath.MOVETO] + [curve[cmd]]*(len(pts)-1))
        ax.add_patch(PathPatch(path, fill=False, edgecolor=stroke, linewidth=width, alpha=.9,
                               linestyle=(0, (5, 3)) if e.get('same_polarity') else 'solid', zorder=1))
        if not e['intra_asu']:
            heading = math.degrees(math.atan2(tangent[1], tangent[0]))
            ax.plot(*mid, marker=(3, 0, heading - 90), markersize=arrow, color=stroke,
                    markeredgecolor=stroke, zorder=2)
            tag = e['symmetry_code'] + (' ' + '/'.join(e['polarities']) if e.get('same_polarity') else '')
            ax.annotate(tag, mid, xytext=(.6*arrow, .4*arrow), textcoords='offset points', fontsize=codes_font,
                        color=stroke, zorder=4,
                        bbox=dict(boxstyle='round,pad=0.1', facecolor='white', edgecolor='none', alpha=.75))
    if resolution == 'chain':
        for n in nodes:
            p = pos[n['id']]
            xs.append(p[0]); ys.append(p[1])
            big = n['kind'] == 'polymer'
            ax.scatter(*p, s=900 if big else 220, color=color[n['id']], edgecolors='#25313b', linewidths=1, zorder=3)
            ax.text(*p, n['label'], ha='center', va='center', color='white', fontsize=11 if big else 7,
                    fontweight='bold', zorder=5)
    else:
        size = max(10., min(60., 2600/max(len(nodes), 1)))
        label_font = max(4., min(7., 600/max(len(nodes), 1)))
        for n in nodes:
            p = pos[n['id']]
            xs.append(p[0]); ys.append(p[1])
            on_rim = n['id'] in angle_of
            ax.scatter(*p, s=size if on_rim else size*1.5, color=color[n['chain_node']], edgecolors='#25313b',
                       linewidths=.5, zorder=3)
            angle = math.atan2(p[1], p[0]) if on_rim else 0.
            flip = on_rim and math.cos(angle) < 0
            q = _scale(p, 1.06) if on_rim else _add(p, (0, -.06))
            text = f"{n['residue']}{n['seq']}{n['insertion']}" if on_rim else n['label']
            ax.text(*q, text, rotation=math.degrees(angle + (math.pi if flip else 0.)) if on_rim else 0.,
                    rotation_mode='anchor', ha='right' if flip else ('left' if on_rim else 'center'),
                    va='center', fontsize=label_font, color='#4a5560', zorder=5)
        for c in chains:
            members = [pos[n['id']] for n in nodes if n['chain_node'] == c and n['id'] in angle_of]
            mean = (sum(p[0] for p in members)/len(members), sum(p[1] for p in members)/len(members))
            point = _scale(mean, 1.42/_norm(mean)) if _norm(mean) > 1e-6 else (0., 1.42)
            xs.append(point[0]*1.08); ys.append(point[1]*1.08)
            name = next(n['chain_label'] for n in nodes if n['chain_node'] == c)
            ax.text(*point, name, ha='center', va='center', fontsize=14, fontweight='bold', color=color[c])
    pad = .32
    ax.set_xlim(min(xs, default=-1.) - pad, max(xs, default=1.) + pad)
    ax.set_ylim(min(ys, default=-1.) - pad, max(ys, default=1.) + pad)
    kinds = [k for k in EDGE_LABELS if any(e['type'] == k for e in graph['edges'])]
    handles = [Line2D([], [], color=COLORS[k], linewidth=width, label=f"{EDGE_LABELS[k]} ({graph['edge_type_counts'][k]})")
               for k in kinds]
    if any(e.get('same_polarity') for e in graph['edges']):
        handles.append(Line2D([], [], color=COLORS['terminal_stack'], linewidth=width, linestyle=(0, (5, 3)),
                              label="Same-polarity (5'-5' or 3'-3') chain-end stack"))
    if handles:
        ax.legend(handles=handles, loc='lower left', fontsize=8, frameon=False)
    if title is None:
        what = 'chain' if resolution == 'chain' else 'residue'
        title = f"{graph['name']}: one ASU at {what} resolution ({graph['space_group']})"
    fig.text(0.03, 1 - 0.25/fig.get_figheight(), title, fontsize=14, fontweight='bold', va='top')
    fig.text(0.03, 1 - 0.6/fig.get_figheight(),
             'Periodicity: ' + '\n'.join(summary) + '\nArrows point to the copy named by the symmetry code '
             'op_(k,l,m): operation, then signed cell shifts along a, b, c; unlabelled edges are inside the ASU.', fontsize=9, va='top', color='#333333')
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', dpi=dpi, facecolor='white', metadata={'Software': None})
    plt.close(fig)
    return buffer.getvalue()
