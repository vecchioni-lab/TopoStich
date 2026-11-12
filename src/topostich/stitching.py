"""Backbone nick and junction stitching for the O3-prime CGD graph.

Both rules are ported from the lab's NAtopo project, without its NetworkX,
NumPy or SciPy dependencies. Gemmi supplies the symmetry operator algebra.

**Nicks** (``nucleotide_interactions.find_backbone_join_candidates`` with the
endpoint rules of ``prepare_nucleotides`` and the polymer links of
``build_atom_graph``):

* A **loose P** belongs to a residue that also has O5' and has no covalent link
  to another residue. A **loose O3'** has its own C3' and no such link.
* Residues are covalently linked within a chain/TER segment when consecutive
  residues have O3'(i)-P(i+1) at 0.8-2.2 A, or when an explicit input
  CONECT/_struct_conn joins the atoms at 0.8-2.3 A.
* A **backbone join** pairs a loose P with a loose O3' closer than
  ``join_cutoff`` (3.5 A) on a *different chain copy*: another segment, or the
  same segment in another crystal image. The **O5' fallback** pairs a
  phosphate-less residue's O5' (own C5') with a loose O3' closer than
  ``o5_cutoff`` (6.0 A); it is stitched only on request.
* Only one-to-one joins between segment termini become edges.

**Junctions** (``find_pi_interactions``/``fit_base_planes`` and
``junctions.find_junction_candidates``, method ``nonbackbone_pi_stack_v2``):

* Base planes are SVD fits of complete ring heavy atoms (9 purine, 6
  pyrimidine). Ring atom names come from TopoStich's NARestraints-derived
  residue catalog, so modified residues map canonical ring positions to their
  actual names; unmapped or ambiguous residues are reported, not guessed.
* A **pi stack** has centroids < 5.5 A, normals within 30 deg, both normal
  separations 2.0-4.5 A and both lateral shifts <= 3.0 A.
* A **junction candidate** is a stack whose two residues are not joined by a
  verified O3'-P bond. Its reason is ``stack_on_face_without_backbone_continuation``
  when, for either endpoint, every verified backbone neighbor lies on the
  opposite face (0.5 A side tolerance); otherwise ``non_backbone_stack``.
* ``junction_scope='face'`` (default) stitches face junctions only: the helix
  continues by stacking past a strand end. ``'all'`` stitches every candidate.

Every symmetry operation and every lattice translation that can satisfy a
cutoff is enumerated exactly (no finite cell window). Candidates that are not
stitched are reported with a status. No atoms are moved, added or repaired.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import math

import gemmi

from .converter import ConversionError, normalized_name

NICK_METHOD = 'natopo_backbone_join_v1'
JUNCTION_METHOD = 'natopo_nonbackbone_pi_stack_v2'
POLYMER_LINK_A = (0.8, 2.2)
EXPLICIT_COVALENT_A = (0.8, 2.3)
ENDPOINT_ATOMS = ('P', "O5'", "C5'", "C3'", "O3'")
DEN = 24
PURINE_RING = ('N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4')
PYRIMIDINE_RING = ('N1', 'C2', 'N3', 'C4', 'C5', 'C6')
RING_ELEMENTS = {'C', 'N', 'O', 'S', 'SE'}
PI_DEFAULTS = dict(centroid_cutoff_A=5.5, max_angle_deg=30.0, max_lateral_shift_A=3.0,
                   min_plane_separation_A=2.0, max_plane_separation_A=4.5, min_side_projection_A=0.5)
JUNCTION_SCOPES = ('face', 'all')


def _positive(value, name):
    if not isinstance(value, (int, float)) or isinstance(value, bool) or not math.isfinite(value) or value <= 0:
        raise ConversionError(f'{name} must be a finite positive distance in Angstrom')
    return float(value)


def _invert3(m):
    a, b, c = m
    det = (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0]))
    if det == 0:
        raise ConversionError('Singular fractional transform')
    return [[(b[1]*c[2]-b[2]*c[1])/det, (a[2]*c[1]-a[1]*c[2])/det, (a[1]*b[2]-a[2]*b[1])/det],
            [(b[2]*c[0]-b[0]*c[2])/det, (a[0]*c[2]-a[2]*c[0])/det, (a[2]*b[0]-a[0]*b[2])/det],
            [(b[0]*c[1]-b[1]*c[0])/det, (a[1]*c[0]-a[0]*c[1])/det, (a[0]*b[1]-a[1]*b[0])/det]]


def _matvec(m, v):
    return [sum(m[i][k]*v[k] for k in range(3)) for i in range(3)]


def _dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def _sub(a, b):
    return [x-y for x, y in zip(a, b)]


def _norm(v):
    return math.sqrt(_dot(v, v))


def _symmetric_eigen(m):
    """Cyclic Jacobi eigen-decomposition of a symmetric 3x3 matrix."""
    a = [row[:] for row in m]
    v = [[1. if i == j else 0. for j in range(3)] for i in range(3)]
    for _ in range(100):
        off = sum(a[i][j]**2 for i in range(3) for j in range(3) if i != j)
        if off < 1e-30:
            break
        for p, q in ((0, 1), (0, 2), (1, 2)):
            if abs(a[p][q]) < 1e-300:
                continue
            theta = (a[q][q]-a[p][p]) / (2*a[p][q])
            t = math.copysign(1., theta) / (abs(theta) + math.sqrt(theta*theta+1))
            c = 1/math.sqrt(t*t+1)
            s = t*c
            for k in range(3):
                akp, akq = a[k][p], a[k][q]
                a[k][p], a[k][q] = c*akp - s*akq, s*akp + c*akq
            for k in range(3):
                apk, aqk = a[p][k], a[q][k]
                a[p][k], a[q][k] = c*apk - s*aqk, s*apk + c*aqk
            for k in range(3):
                vkp, vkq = v[k][p], v[k][q]
                v[k][p], v[k][q] = c*vkp - s*vkq, s*vkp + c*vkq
    values = [a[i][i] for i in range(3)]
    order = sorted(range(3), key=lambda i: values[i])
    return [values[i] for i in order], [[v[k][i] for k in range(3)] for i in order]


def _atom_identity(atom):
    return dict(chain=atom.chain, seq=atom.seq, insertion=atom.insertion, residue=atom.residue,
                atom=atom.name, altloc=atom.altloc, segment=atom.segment, serial=atom.serial,
                label_asym=atom.label_asym, label_seq=atom.label_seq, description=atom.description)


def _residue_identity(atom):
    return dict(chain=atom.chain, seq=atom.seq, insertion=atom.insertion, residue=atom.residue,
                segment=atom.segment, label_asym=atom.label_asym, label_seq=atom.label_seq,
                description=f"segment={atom.segment} chain={atom.chain or '.'} "
                            f"residue={atom.residue}:{atom.seq}{atom.insertion}")


class _Frame:
    """Fractional <-> Cartesian frame of the input, including any SCALE shift."""

    def __init__(self, transform):
        self.F = [list(row[:3]) for row in transform]
        self.u = [row[3] for row in transform]
        self.A = _invert3(self.F)

    def fractional(self, xyz):
        return [v+s for v, s in zip(_matvec(self.F, xyz), self.u)]

    def cartesian(self, frac):
        return _matvec(self.A, [f-s for f, s in zip(frac, self.u)])

    def length(self, delta_frac):
        return _norm(_matvec(self.A, delta_frac))

    def translation_ranges(self, delta, cutoff):
        # |frac_k| <= |F_k| * |cart| bounds every lattice shift that can satisfy the cutoff.
        ranges = []
        for k in range(3):
            reach = cutoff*_norm(self.F[k]) + 1e-9
            ranges.append(range(math.ceil(-delta[k]-reach), math.floor(-delta[k]+reach)+1))
        return ranges

    def images(self, operations, a, b, cutoff):
        """Yield (op, shift, distance) with |op(b)+shift - a| < cutoff (fractional a, b)."""
        for op in operations:
            delta = _sub(op.apply_to_xyz(b), a)
            rx, ry, rz = self.translation_ranges(delta, cutoff)
            for shift in ((i, j, k) for i in rx for j in ry for k in rz):
                d = self.length([x+n for x, n in zip(delta, shift)])
                if d < cutoff:
                    yield op, shift, d


def _image_op(op, shift):
    g = gemmi.Op(op.triplet())
    g.tran = [t + DEN*n for t, n in zip(g.tran, shift)]
    return g


IDENTITY_OP = gemmi.Op('x,y,z')
IDENTITY = IDENTITY_OP.triplet()


class StitchContext:
    """Residues, verified backbone links, nodes and existing edges of one graph."""

    def __init__(self, pdb, *, altloc, node_atoms, node_labels, transform, operations, trace_pairs):
        self.pdb, self.frame, self.operations = pdb, _Frame(transform), operations
        self.node_atoms, self.node_labels = node_atoms, node_labels
        self.node_of = {atom.residue_key: i for i, atom in enumerate(node_atoms)}
        self._read_residues(altloc)
        self._read_links()
        # Existing neighbors per node as (neighbor index, rounded image coordinates).
        self.occupied = defaultdict(set)
        for i, j in trace_pairs:
            self._occupy(i, j, IDENTITY_OP)

    def _read_residues(self, altloc):
        residues, order, all_atoms = {}, [], defaultdict(dict)
        for atom in self.pdb.atoms:
            if altloc is not None and atom.altloc not in ('', altloc):
                continue
            key = atom.residue_key
            if key not in residues:
                residues[key] = {}
                order.append(key)
            name = normalized_name(atom.name)
            if name in all_atoms[key]:
                if name in ENDPOINT_ATOMS:
                    raise ConversionError(f'Duplicate backbone atom {atom.description}; choose --altloc')
                all_atoms[key][name] = None  # ambiguous alternate; resolved per use
                continue
            all_atoms[key][name] = atom
            if name in ENDPOINT_ATOMS:
                residues[key][name] = atom
        self.residues, self.all_atoms = residues, all_atoms
        self.nucleotide = [key for key in order if any(n in residues[key] for n in ENDPOINT_ATOMS)]
        self.order = order
        self.by_segment = defaultdict(list)
        for key in self.nucleotide:
            self.by_segment[key[0]].append(key)

    def _read_links(self):
        residues = self.residues
        linked = defaultdict(set)   # (residue key, atom name) -> {(other key, other atom name)}
        for keys in self.by_segment.values():
            for left, right in zip(keys, keys[1:]):
                o3, p = residues[left].get("O3'"), residues[right].get('P')
                if o3 is not None and p is not None and POLYMER_LINK_A[0] <= math.dist(o3.xyz, p.xyz) <= POLYMER_LINK_A[1]:
                    linked[left, "O3'"].add((right, 'P'))
                    linked[right, 'P'].add((left, "O3'"))
        serials = defaultdict(list)
        for key, atoms in residues.items():
            for name, atom in atoms.items():
                serials[atom.serial].append((key, name, atom))
        by_serial = {atom.serial: atom for atom in self.pdb.atoms}
        explicit = []
        for a, b in self.pdb.conect:
            for here, there in ((a, b), (b, a)):
                other = by_serial.get(there)
                for key, name, atom in serials.get(here, ()):
                    if other is None or other.residue_key == key:
                        continue
                    if not EXPLICIT_COVALENT_A[0] <= math.dist(atom.xyz, other.xyz) <= EXPLICIT_COVALENT_A[1]:
                        continue
                    linked[key, name].add((other.residue_key, normalized_name(other.name)))
                    if here == a:
                        explicit.extend(((key, name, atom), right) for right in serials.get(there, ())
                                        if right[0] != key)
        self.linked, self.explicit = linked, explicit
        # Verified O3'-P backbone adjacency between residues (same image).
        self.bonded = defaultdict(set)
        for (key, name), partners in linked.items():
            for other, other_name in partners:
                if {name, other_name} == {"O3'", 'P'}:
                    self.bonded[key].add(other)
                    self.bonded[other].add(key)

    def _occupy(self, source, target, g):
        f_src = self.frame.fractional(self.node_atoms[source].xyz)
        f_tgt = self.frame.fractional(self.node_atoms[target].xyz)
        forward = g.apply_to_xyz(f_tgt)
        reverse = g.inverse().apply_to_xyz(f_src)
        key_f = (target, tuple(round(x, 6) for x in forward))
        key_r = (source, tuple(round(x, 6) for x in reverse))
        if (source == target and g.triplet() == IDENTITY) or key_f in self.occupied[source] \
                or key_r in self.occupied[target]:
            return None
        self.occupied[source].add(key_f)
        self.occupied[target].add(key_r)
        return forward, reverse

    def add_edge(self, source, target, g, *, kind, distance, source_atom, target_atom, evidence):
        """Register a periodic edge; return its record, or None if it duplicates one."""
        images = self._occupy(source, target, g)
        if images is None:
            return None
        forward, reverse = images
        return dict(source=self.node_labels[source], target=self.node_labels[target], kind=kind,
                    operation=_wrapped(g).triplet(), lattice_translation=_shift(g),
                    image_operation=g.triplet(), reverse_image_operation=g.inverse().triplet(),
                    target_image_fractional=list(forward), reverse_image_fractional=list(reverse),
                    target_image_cartesian=self.frame.cartesian(forward),
                    reverse_image_cartesian=self.frame.cartesian(reverse), distance_A=distance,
                    source_description=_residue_identity(source_atom)['description'],
                    target_description=_residue_identity(target_atom)['description'],
                    evidence=evidence)


def _wrapped(g):
    w = gemmi.Op(g.triplet())
    w.tran = [t % DEN for t in w.tran]
    return w


def _shift(g):
    return [t // DEN for t in g.tran]


# ----------------------------------------------------------------------------- nicks

def _nick_endpoints(ctx):
    loose_p, loose_o3, loose_o5 = [], [], []
    for key in ctx.nucleotide:
        atoms = ctx.residues[key]
        if 'P' in atoms and "O5'" in atoms and not ctx.linked.get((key, 'P')):
            loose_p.append(key)
        if "O3'" in atoms and "C3'" in atoms and not ctx.linked.get((key, "O3'")):
            loose_o3.append(key)
        if "O5'" in atoms and 'P' not in atoms and "C5'" in atoms and not ctx.linked.get((key, "O5'")):
            loose_o5.append(key)
    return loose_p, loose_o3, loose_o5


def find_backbone_joins(ctx, *, join_cutoff=3.5, o5_cutoff=6.0, include_o5_fallback=False):
    """Return (stitched edges, candidate records, summary) for nicks.

    Each stitched edge joins the upstream 3' end (source, stored position) to the
    downstream strand's first O3' node (target) in image g = operation + shift.
    """
    join_cutoff = _positive(join_cutoff, 'join_cutoff')
    o5_cutoff = _positive(o5_cutoff, 'o5_cutoff')
    residues, frame = ctx.residues, ctx.frame
    loose_p, loose_o3, loose_o5 = _nick_endpoints(ctx)

    def frac(key, name):
        return frame.fractional(residues[key][name].xyz)

    candidates = []
    searches = (('P', loose_p, join_cutoff, 'backbone_join'),
                ("O5'", loose_o5, o5_cutoff, 'backbone_join_o5'))
    for o3_key in loose_o3:
        a = frac(o3_key, "O3'")
        for atom_name, five_primes, cutoff, kind in searches:
            for p_key in five_primes:
                for op, shift, d in frame.images(ctx.operations, a, frac(p_key, atom_name), cutoff):
                    same_image = op.triplet() == IDENTITY and shift == (0, 0, 0)
                    if same_image and o3_key[0] == p_key[0]:
                        continue  # same chain copy, as in NAtopo's interchain-only rule
                    candidates.append(dict(kind=kind, o3_key=o3_key, p_key=p_key, p_name=atom_name,
                                           op=op, shift=list(shift), distance=d))
    for (o3_key, o3_name, o3_atom), (p_key, p_name, p_atom) in ctx.explicit:
        if (o3_name, p_name) == ('P', "O3'"):
            o3_key, p_key, o3_atom, p_atom = p_key, o3_key, p_atom, o3_atom
        elif (o3_name, p_name) != ("O3'", 'P'):
            continue
        if o3_key[0] == p_key[0] or any(c['kind'] == 'deposited_link' and c['o3_key'] == o3_key
                                        and c['p_key'] == p_key for c in candidates):
            continue
        candidates.append(dict(kind='deposited_link', o3_key=o3_key, p_key=p_key, p_name='P',
                               op=IDENTITY_OP, shift=[0, 0, 0], distance=math.dist(o3_atom.xyz, p_atom.xyz)))

    enabled = {'backbone_join', 'deposited_link'} | ({'backbone_join_o5'} if include_o5_fallback else set())
    partners_3, partners_5 = Counter(), Counter()
    for c in candidates:
        if c['kind'] in enabled:
            partners_3[c['o3_key']] += 1
            partners_5[c['p_key']] += 1

    def first_node_from(key):
        keys = ctx.by_segment[key[0]]
        for later in keys[keys.index(key):]:
            if later in ctx.node_of:
                return ctx.node_of[later]
        return None

    stitched, records = [], []
    for c in sorted(candidates, key=lambda c: (c['distance'], str(c['o3_key']), str(c['p_key']))):
        o3_key, p_key = c['o3_key'], c['p_key']
        g = _image_op(c['op'], c['shift'])
        source, target = ctx.node_of.get(o3_key), first_node_from(p_key)
        status = 'stitched'
        if c['kind'] not in enabled:
            status = 'o5_fallback_disabled'
        elif ctx.by_segment[o3_key[0]][-1] != o3_key or ctx.by_segment[p_key[0]][0] != p_key:
            status = 'segment_internal'
        elif partners_3[o3_key] > 1 or partners_5[p_key] > 1:
            status = 'ambiguous'
        elif source is None or target is None:
            status = 'no_graph_node'
        o3_atom, p_atom = residues[o3_key]["O3'"], residues[p_key][c['p_name']]
        record = dict(kind=c['kind'], status=status, method=NICK_METHOD, distance_A=c['distance'],
                      cutoff_A=(o5_cutoff if c['kind'] == 'backbone_join_o5'
                                else join_cutoff if c['kind'] == 'backbone_join' else None),
                      operation=c['op'].triplet(), lattice_translation=c['shift'], image_operation=g.triplet(),
                      o3_atom=_atom_identity(o3_atom), p_atom=_atom_identity(p_atom),
                      source=ctx.node_labels[source] if source is not None else None,
                      target=ctx.node_labels[target] if target is not None else None)
        if status == 'stitched':
            edge = ctx.add_edge(source, target, g, kind=c['kind'], distance=c['distance'],
                                source_atom=o3_atom, target_atom=p_atom,
                                evidence=dict(o3_atom=record['o3_atom'], p_atom=record['p_atom']))
            if edge is None:
                record['status'] = 'duplicate_edge'
            else:
                stitched.append(edge)
        records.append(record)
    summary = dict(method=NICK_METHOD,
                   source='Ported from NAtopo nucleotide_interactions.find_backbone_join_candidates',
                   join_cutoff_A=join_cutoff, o5_cutoff_A=o5_cutoff, include_o5_fallback=include_o5_fallback,
                   polymer_link_A=list(POLYMER_LINK_A), explicit_covalent_A=list(EXPLICIT_COVALENT_A),
                   image_search='every symmetry operation and every lattice translation within the cutoff',
                   loose_endpoints=dict(P=len(loose_p), O3=len(loose_o3), O5_missing_phosphate=len(loose_o5)),
                   candidate_status_counts=dict(sorted(Counter(r['status'] for r in records).items())))
    return stitched, records, summary


# ------------------------------------------------------------------------- junctions

def _ring_names(code, observed):
    """Canonical ring -> actual atom names for one residue, or an issue string."""
    from .nucleotides import _matching_definitions
    records, status = _matching_definitions(code, set(observed))
    if status == 'unmapped':
        return None, 'unmapped_component'
    if status == 'ambiguous':
        return None, 'ambiguous_ring_atom_mapping'
    positions = dict(records[0].positions)
    analog = records[0].base_analog.upper()
    if analog in ('A', 'G') or (analog not in ('C', 'T', 'U') and all(k in positions for k in PURINE_RING)):
        ring = PURINE_RING
    else:
        ring = PYRIMIDINE_RING
    names = [positions.get(k, '') for k in ring]
    if any(not n or n == '/' for n in names):
        return None, 'missing_ring_atom_mapping'
    if len(set(names)) != len(names):
        return None, 'ambiguous_ring_atom_mapping'
    return names, None


def base_planes(ctx):
    """(planes keyed by residue key, diagnostics) for nucleotide-like residues."""
    planes, diagnostics = {}, []
    for key in ctx.order:
        atoms = ctx.all_atoms[key]
        if key not in ctx.nucleotide and "C1'" not in atoms:
            continue
        sample = next(a for a in atoms.values() if a is not None) if any(atoms.values()) else None
        if sample is None:
            continue
        code = sample.residue
        detail = dict(residue=code, description=_residue_identity(sample)['description'])
        names, issue = _ring_names(code, atoms)
        if issue is None:
            missing = [n for n in names if n not in atoms]
            ambiguous = [n for n in names if n in atoms and atoms[n] is None]
            if ambiguous:
                raise ConversionError(f"Base ring atoms {', '.join(ambiguous)} have alternate locations in "
                                      f"{detail['description']}; choose --altloc")
            if missing:
                issue = 'incomplete_base_ring'
                detail['atoms'] = missing
            elif any(atoms[n].element.upper() not in RING_ELEMENTS for n in names):
                issue = 'unsupported_ring_atom_element'
        if issue is not None:
            diagnostics.append(dict(detail, issue=issue))
            continue
        points = [list(atoms[n].xyz) for n in names]
        center = [sum(p[k] for p in points)/len(points) for k in range(3)]
        centered = [_sub(p, center) for p in points]
        scatter = [[sum(p[i]*p[j] for p in centered) for j in range(3)] for i in range(3)]
        values, vectors = _symmetric_eigen(scatter)
        if math.sqrt(max(values[1], 0.)) <= 1e-8:
            diagnostics.append(dict(detail, issue='degenerate_base_ring'))
            continue
        normal = vectors[0]
        planes[key] = dict(centroid=center, normal=normal, ring_atoms=names, atom=sample,
                           rms_deviation_A=math.sqrt(sum(_dot(p, normal)**2 for p in centered)/len(points)))
    return planes, diagnostics


def _transform_plane(frame, g, plane):
    """Cartesian centroid and normal of a plane moved by fractional operation g."""
    c = plane['centroid']
    moved = frame.cartesian(g.apply_to_xyz(frame.fractional(c)))
    tip = frame.cartesian(g.apply_to_xyz(frame.fractional([x+y for x, y in zip(c, plane['normal'])])))
    return moved, _sub(tip, moved)


def find_junction_stacks(ctx, *, scope='face'):
    """Return (stitched edges, junction candidate records, summary)."""
    if scope not in JUNCTION_SCOPES:
        raise ConversionError(f'junction_scope must be one of {", ".join(JUNCTION_SCOPES)}')
    p = PI_DEFAULTS
    frame = ctx.frame
    planes, diagnostics = base_planes(ctx)
    keys = list(planes)
    candidates, excluded = [], 0
    for ia, a in enumerate(keys):
        fa = frame.fractional(planes[a]['centroid'])
        for b in keys[ia:]:
            for op, shift, d in frame.images(ctx.operations, fa, frame.fractional(planes[b]['centroid']),
                                             p['centroid_cutoff_A']):
                g = _image_op(op, shift)
                if a == b:
                    if g.triplet() == IDENTITY or g.triplet() > g.inverse().triplet():
                        continue  # self, or the reverse description of the same orbit
                cb, nb = _transform_plane(frame, g, planes[b])
                ca, na = planes[a]['centroid'], planes[a]['normal']
                delta = _sub(cb, ca)
                angle = math.degrees(math.acos(min(1., abs(_dot(na, nb)))))
                heights = [abs(_dot(delta, n)) for n in (na, nb)]
                shifts = [_norm(_sub(delta, [_dot(delta, n)*x for x in n])) for n in (na, nb)]
                if (angle > p['max_angle_deg'] or max(shifts) > p['max_lateral_shift_A']
                        or min(heights) < p['min_plane_separation_A'] or max(heights) > p['max_plane_separation_A']):
                    continue
                if g.triplet() == IDENTITY and b in ctx.bonded[a]:
                    excluded += 1  # verified backbone neighbors stack normally
                    continue
                candidates.append(dict(a=a, b=b, g=g, distance=d, angle=angle, heights=heights, shifts=shifts))

    def side(height):
        tol = p['min_side_projection_A']
        return None if height is None or abs(height) <= tol else (1 if height > 0 else -1)

    def context(anchor, partner_centroid):
        """Face context of one endpoint; neighbors share the anchor's image."""
        plane = planes[anchor]
        pi_height = _dot(_sub(partner_centroid, plane['centroid']), plane['normal'])
        pi_side = side(pi_height)
        neighbors = []
        for other in sorted(ctx.bonded[anchor], key=str):
            other_plane = planes.get(other)
            height = None if other_plane is None else _dot(_sub(other_plane['centroid'], plane['centroid']), plane['normal'])
            nside = side(height)
            neighbors.append(dict(residue=_residue_identity(next(x for x in ctx.all_atoms[other].values() if x))['description'],
                                  signed_projection_A=height,
                                  same_side_as_pi=None if pi_side is None or nside is None else pi_side == nside))
        empty = bool(neighbors) and all(n['same_side_as_pi'] is False for n in neighbors)
        return dict(anchor=_residue_identity(plane['atom'])['description'], pi_projection_A=pi_height,
                    verified_backbone_neighbors=neighbors, stack_face_without_verified_backbone=empty)

    stitched, records = [], []
    for c in sorted(candidates, key=lambda c: (c['distance'], str(c['a']), str(c['b']), c['g'].triplet())):
        a, b, g = c['a'], c['b'], c['g']
        cb, _ = _transform_plane(frame, g, planes[b])
        ca_in_b, _ = _transform_plane(frame, g.inverse(), planes[a])
        contexts = [context(a, cb), context(b, ca_in_b)]
        reason = ('stack_on_face_without_backbone_continuation'
                  if any(x['stack_face_without_verified_backbone'] for x in contexts) else 'non_backbone_stack')
        source, target = ctx.node_of.get(a), ctx.node_of.get(b)
        status = 'stitched'
        if scope == 'face' and reason == 'non_backbone_stack':
            status = 'not_face_junction'
        elif source is None or target is None:
            status = 'no_graph_node'
        evidence = dict(reason=reason, centroid_distance_A=c['distance'], plane_angle_deg=c['angle'],
                        plane_separation_A=c['heights'], lateral_shift_A=c['shifts'])
        record = dict(kind='junction_stack', status=status, method=JUNCTION_METHOD, reason=reason, scope=scope,
                      distance_A=c['distance'], operation=_wrapped(g).triplet(), lattice_translation=_shift(g),
                      image_operation=g.triplet(), residue1=_residue_identity(planes[a]['atom']),
                      residue2=_residue_identity(planes[b]['atom']),
                      source=ctx.node_labels[source] if source is not None else None,
                      target=ctx.node_labels[target] if target is not None else None,
                      measurements={k: v for k, v in evidence.items() if k != 'reason'}, backbone_context=contexts)
        if status == 'stitched':
            edge = ctx.add_edge(source, target, g, kind='junction_stack', distance=c['distance'],
                                source_atom=planes[a]['atom'], target_atom=planes[b]['atom'], evidence=evidence)
            if edge is None:
                record['status'] = 'duplicate_edge'
            else:
                stitched.append(edge)
        records.append(record)
    summary = dict(method=JUNCTION_METHOD,
                   source='Ported from NAtopo nucleotide_interactions.find_pi_interactions and junctions.find_junction_candidates',
                   scope=scope, **p, ring_atom_source='TopoStich NARestraints-derived residue catalog',
                   image_search='every symmetry operation and every lattice translation within the centroid cutoff',
                   base_planes=len(planes), verified_backbone_stacks_excluded=excluded,
                   plane_diagnostics=diagnostics,
                   candidate_status_counts=dict(sorted(Counter(r['status'] for r in records).items())),
                   reason_counts=dict(sorted(Counter(r['reason'] for r in records).items())))
    return stitched, records, summary
