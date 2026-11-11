"""PDB to core CIF, PDBx/mmCIF, or PDB without chemistry/distance inference.

Coordinates stay unwrapped. Trace edges are generic links, not restraints.
Only the Python standard library is required.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import math
from pathlib import Path
import re
from typing import Sequence

from . import __version__

BACKBONE_ORDER = ("P", "O5'", "C5'", "C4'", "C3'", "O3'")
FORMATS = ("cif", "pdb", "mmcif")


class ConversionError(ValueError):
    """An input cannot be converted without guessing or losing identity."""


def normalized_name(name: str) -> str:
    return name.strip().replace("*", "'").replace("\u2032", "'")


@dataclass(frozen=True)
class Atom:
    serial: str
    name: str
    altloc: str
    residue: str
    chain: str
    seq: str
    insertion: str
    segment: int
    xyz: tuple[float, float, float]
    occupancy: float | None
    b_iso: float | None
    element: str
    raw: str
    label_asym: str = ""
    label_seq: str = ""
    label_comp: str = ""
    label_atom: str = ""
    record_type: str = "ATOM"

    @property
    def residue_key(self) -> tuple:
        return self.segment, self.chain, self.seq, self.insertion, self.residue, self.label_seq

    @property
    def description(self) -> str:
        return (f"segment={self.segment} chain={self.chain or '.'} "
                f"residue={self.residue}:{self.seq}{self.insertion} "
                f"atom={self.name} altloc={self.altloc or '.'}")


@dataclass
class Pdb:
    atoms: list[Atom]
    cell: tuple[float, ...] | None
    space_group: str
    scale: tuple[tuple[float, ...], ...] | None
    symops: dict[int, str]
    conect: list[tuple[str, str]]
    headers: list[str]
    model_count: int
    source_format: str = "pdb"
    connection_warnings: list[str] = field(default_factory=list)
    scale_is_pdb_rounded: bool = True


@dataclass(frozen=True)
class ConversionResult:
    schema_version: int
    version: str
    source_name: str
    source_sha256: str
    output_name: str
    output_sha256: str
    output_format: str
    model: int
    atom_selection: list[str]
    atom_count: int
    bond_count: int
    segment_count: int
    missing_count: int
    missing: list[str]
    warnings: list[str]


def _number(value: str, context: str, optional: bool = False) -> float | None:
    if optional and not value.strip():
        return None
    try:
        result = float(value)
    except ValueError as exc:
        raise ConversionError(f"Invalid {context}: {value!r}") from exc
    if not math.isfinite(result):
        raise ConversionError(f"Non-finite {context}")
    return result


def _element(line: str) -> str:
    explicit = line[76:78].strip()
    if explicit:
        if not re.fullmatch(r"[A-Za-z]{1,2}", explicit):
            raise ConversionError(f"Invalid element: {explicit!r}")
        return explicit.capitalize()
    field = line[12:16]
    # PDB alignment distinguishes e.g. alpha carbon ' CA ' from calcium 'CA  '.
    if field[:1].isspace() or field[:1].isdigit():
        letters = re.sub(r"[^A-Za-z]", "", field)
        return letters[0].upper() if letters else "X"
    letters = re.match(r"[A-Za-z]+", field)
    return letters[0][:2].capitalize() if letters else "X"


def read_pdb(text: str, model: int = 1) -> Pdb:
    lines = text.splitlines()
    model_count = sum(line[:6].strip() == "MODEL" for line in lines) or 1
    explicit_models = any(line[:6].strip() == "MODEL" for line in lines)
    if not 1 <= model <= model_count:
        raise ConversionError(f"Model must be in 1..{model_count} (one-based file order)")
    current_model = 0 if explicit_models else 1
    in_model = not explicit_models
    segment = 0
    last_chain = None
    atoms: list[Atom] = []
    cell = None
    group = ""
    scales = {}
    symops = {}
    conect = []
    headers = []
    for line_number, raw in enumerate(lines, 1):
        line = raw.ljust(80)
        record = line[:6].strip()
        if record == "CRYST1":
            new_cell = tuple(_number(line[a:b], "CRYST1") for a, b in
                             ((6, 15), (15, 24), (24, 33), (33, 40), (40, 47), (47, 54)))
            new_group = line[55:66].strip()
            if cell is not None and (cell != new_cell or group != new_group):
                raise ConversionError("Conflicting CRYST1 records")
            cell, group = new_cell, new_group
            headers.append(raw)
        elif record in ("SCALE1", "SCALE2", "SCALE3"):
            row = tuple(_number(line[a:b], record) for a, b in ((10, 20), (20, 30), (30, 40), (45, 55)))
            index = int(record[-1]) - 1
            if index in scales and scales[index] != row:
                raise ConversionError("Conflicting SCALE records")
            scales[index] = row
            headers.append(raw)
        elif raw.startswith("REMARK 290"):
            match = re.fullmatch(r"REMARK 290\s+(\d+)555\s+([xyzXYZ0-9,+/ .-]+)\s*", raw)
            if match:
                index = int(match[1])
                operation = match[2].replace(" ", "").lower()
                if operation.count(",") != 2:
                    raise ConversionError("Malformed fractional symmetry operator in REMARK 290")
                if index in symops and symops[index] != operation:
                    raise ConversionError("Conflicting fractional symmetry operators")
                symops[index] = operation
            headers.append(raw)
        elif record == "MODEL":
            current_model += 1
            in_model = True
            last_chain = None
        elif record == "ENDMDL":
            in_model = False
            last_chain = None
        elif record == "TER" and current_model == model:
            last_chain = None
        elif record == "CONECT":
            # Archive CONECT often follows ENDMDL; model-local records must match.
            if not explicit_models or not in_model or current_model == model:
                source = line[6:11].strip()
                conect.extend((source, line[pos:pos+5].strip())
                              for pos in range(11, len(raw), 5) if line[pos:pos+5].strip())
        elif record in ("ATOM", "HETATM") and in_model and current_model == model:
            if len(raw) < 54:
                raise ConversionError(f"Truncated coordinate record on line {line_number}")
            chain = line[21].strip()
            if chain != last_chain:
                segment += 1
            last_chain = chain
            coords = tuple(_number(line[a:b], f"coordinate on line {line_number}")
                           for a, b in ((30, 38), (38, 46), (46, 54)))
            name, residue, seq = line[12:16].strip(), line[17:20].strip(), line[22:26].strip()
            if not name or not residue or not seq:
                raise ConversionError(f"Missing atom/residue identity on line {line_number}")
            atoms.append(Atom(line[6:11].strip(), name, line[16].strip(), residue, chain,
                              seq, line[26].strip(), segment, coords,
                              _number(line[54:60], "occupancy", True),
                              _number(line[60:66], "B factor", True), _element(line), raw))
    if scales and len(scales) != 3:
        raise ConversionError("Incomplete SCALE records: all three rows are required")
    if not atoms:
        raise ConversionError("No atoms in selected PDB model")
    return Pdb(atoms, cell, group, tuple(scales[i] for i in range(3)) if scales else None,
               symops, conect, headers, model_count)


def _select(pdb: Pdb, names: Sequence[str] | None, altloc: str | None):
    if altloc is not None and (len(altloc) != 1 or not altloc.isalnum()):
        raise ConversionError("--altloc must be one letter or digit")
    requested = {normalized_name(name) for name in names or []}
    unknown = requested.difference(BACKBONE_ORDER)
    if unknown:
        from .nucleotides import SUGAR_ATOMS, classify_atom, header_base_atoms
        residues = {}
        for atom in pdb.atoms:
            residues.setdefault(atom.residue_key, []).append(atom.name)
        base_names = unknown.intersection(header_base_atoms())
        for atom in pdb.atoms:
            if normalized_name(atom.name) in unknown:
                role = classify_atom(atom.residue, atom.name, observed_atoms=residues[atom.residue_key])
                if 'base' in role.candidate_roles:
                    base_names.add(normalized_name(atom.name))
        if base_names:
            raise ConversionError(f"Rejected nucleobase atom selection: {', '.join(sorted(base_names))}. "
                                  "The trace accepts the sugar/backbone path only; omit -atom for all-atom conversion.")
        sugars = unknown.intersection(SUGAR_ATOMS)
        if sugars:
            raise ConversionError(f"Recognized sugar atoms {', '.join(sorted(sugars))} are outside the linear backbone path. "
                                  "Sugar-ring/ribbon connections are not implemented.")
        raise ConversionError(f"Unsupported backbone selections: {', '.join(sorted(unknown))}. "
                              f"Use {', '.join(BACKBONE_ORDER)}; omit -atom to convert all atoms.")
    order = [name for name in BACKBONE_ORDER if name in requested]
    retained = [a for a in pdb.atoms if altloc is None or a.altloc in ("", altloc)]
    missing = []
    if requested:
        residues = {}
        original_residues = {}
        for atom in pdb.atoms:
            residues.setdefault(atom.residue_key, [])
            original_residues.setdefault(atom.residue_key, []).append(atom)
        for atom in retained:
            if normalized_name(atom.name) in requested:
                residues[atom.residue_key].append(atom)
        selected = []
        for key, residue_atoms in residues.items():
            by_name = {}
            for atom in residue_atoms:
                if atom.altloc and altloc is None:
                    raise ConversionError("Selected atoms have alternate locations; choose --altloc A (or another ID)")
                name = normalized_name(atom.name)
                if name in by_name:
                    raise ConversionError(f"Duplicate selected atom: {atom.description}")
                by_name[name] = atom
            # Non-polymer waters/ions lacking any backbone atom are not missing residues.
            original = original_residues[key]
            is_backbone_residue = any(normalized_name(a.name) in BACKBONE_ORDER
                                     or normalized_name(a.name) in ("C1'", "N9", "N1") for a in original)
            for name in order:
                if name in by_name:
                    selected.append(by_name[name])
                elif is_backbone_residue:
                    missing.append(f"segment={key[0]} chain={key[1] or '.'} "
                                   f"residue={key[4]}:{key[2]}{key[3]} atom={name}")
        edges = [(i-1, i) for i in range(1, len(selected))
                 if selected[i-1].segment == selected[i].segment]
    else:
        selected = retained
        serials = {}
        for i, atom in enumerate(selected):
            serials.setdefault(atom.serial, []).append(i)
        pairs = set()
        for left, right in pdb.conect:
            if left not in serials or right not in serials:
                continue
            if len(serials[left]) != 1 or len(serials[right]) != 1:
                raise ConversionError("Ambiguous CONECT endpoint: repeated PDB serial")
            pair = tuple(sorted((serials[left][0], serials[right][0])))
            if pair[0] != pair[1]:
                pairs.add(pair)
        edges = sorted(pairs)
    if not selected:
        raise ConversionError("No atoms matched the selection")
    return selected, edges, order, missing


def _fractional_transform(pdb: Pdb):
    if pdb.cell is None or not pdb.space_group:
        raise ConversionError("Core CIF requires cell and space group (PDB CRYST1 or mmCIF cell/symmetry)")
    a, b, c, alpha, beta, gamma = pdb.cell
    if min(a, b, c) <= 0 or any(not 0 < angle < 180 for angle in (alpha, beta, gamma)):
        raise ConversionError("Invalid CRYST1 cell dimensions/angles")
    if pdb.cell == (1., 1., 1., 90., 90., 90.) and pdb.space_group.replace(" ", "") == "P1":
        raise ConversionError("CRYST1 is the PDB non-crystallographic placeholder; a real cell is required")
    ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
    sg = math.sin(math.radians(gamma))
    volume_factor = 1 - ca*ca - cb*cb - cg*cg + 2*ca*cb*cg
    if volume_factor <= 1e-12:
        raise ConversionError("Degenerate CRYST1 cell")
    # Invert the conventional upper-triangular fractional-to-Cartesian matrix.
    by, cx, cy = b*sg, c*cb, c*(ca-cb*cg)/sg
    cz = c*math.sqrt(volume_factor)/sg
    canonical = ((1/a, -b*cg/(a*by), (b*cg*cy/by-cx)/(a*cz), 0.),
                 (0., 1/by, -cy/(by*cz), 0.), (0., 0., 1/cz, 0.))
    if pdb.scale is None:
        return canonical
    scale = pdb.scale
    det = (scale[0][0]*(scale[1][1]*scale[2][2]-scale[1][2]*scale[2][1])
           - scale[0][1]*(scale[1][0]*scale[2][2]-scale[1][2]*scale[2][0])
           + scale[0][2]*(scale[1][0]*scale[2][1]-scale[1][1]*scale[2][0]))
    if det == 0:
        raise ConversionError("Singular SCALE matrix")
    if not pdb.scale_is_pdb_rounded:
        return scale
    # PDB SCALE has only six decimals. Use precise cell-derived values when
    # it is simply the rounded conventional transform; retain rotated/shifted frames.
    if all(abs(scale[i][j]-canonical[i][j]) <= (5.01e-7 if j < 3 else 5.01e-6)
           for i in range(3) for j in range(4)):
        return canonical
    return scale


def _quote(value: str) -> str:
    if not value or "\n" in value or "\r" in value:
        raise ConversionError("Invalid scalar CIF text")
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_+.,:/-]*", value) and not value.lower().startswith(
            ("data_", "save_", "loop_", "stop_", "global_")):
        return value
    if "'" not in value:
        return "'" + value + "'"
    if '"' not in value:
        return '"' + value + '"'
    raise ConversionError("CIF identifier contains both quote characters")


def _num(value: float | None) -> str:
    return "?" if value is None else f"{value:.10g}"


def _loop(tags, rows):
    return "loop_\n" + "\n".join(tags) + "\n" + "\n".join(" ".join(row) for row in rows) + "\n\n"


def _group(pdb: Pdb) -> str:
    # wwPDB's H lattice denotes rhombohedral symmetry in hexagonal axes.
    return re.sub(r"^H(?=\s|[0-9-])", "R", pdb.space_group)


def _cif_header(pdb: Pdb, mmcif: bool) -> str:
    text = "data_topostich\n\n"
    if pdb.cell:
        for name, value in zip(("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma"), pdb.cell):
            text += f"_cell{'.' if mmcif else '_'}{name} {_num(value)}\n"
    if pdb.space_group:
        text += f"_symmetry{'.' if mmcif else '_'}space_group_name_H-M {_quote(_group(pdb))}\n"
    text += "\n"
    symops = pdb.symops or ({1: "x,y,z"} if pdb.space_group.replace(" ", "") == "P1" else {})
    if symops:
        sep = "." if mmcif else "_"
        text += _loop([f"_space_group_symop{sep}id", f"_space_group_symop{sep}operation_xyz"],
                      [[str(i), _quote(op)] for i, op in sorted(symops.items())])
    return text


def _core_cif(pdb: Pdb, atoms: list[Atom], edges, trace: bool) -> str:
    transform = _fractional_transform(pdb)
    text = _cif_header(pdb, False)
    labels = [f"{a.element}{i+1}" for i, a in enumerate(atoms)]
    text += "".join(f"# {label}: {atom.description}\n" for label, atom in zip(labels, atoms)) + "\n"
    alt_ids = {alt: str(i+1) for i, alt in enumerate(sorted({a.altloc for a in atoms if a.altloc}))}
    residues = {key: str(i+1) for i, key in enumerate(dict.fromkeys(a.residue_key for a in atoms))}
    rows = []
    for label, atom in zip(labels, atoms):
        frac = [sum(row[i]*atom.xyz[i] for i in range(3)) + row[3] for row in transform]
        rows.append([label, atom.element, *(_num(v) for v in frac), _num(atom.occupancy),
                     _num(None if atom.b_iso is None else atom.b_iso/(8*math.pi**2)),
                     residues[atom.residue_key] if atom.altloc else ".", alt_ids.get(atom.altloc, ".")])
    text += _loop(["_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
                   "_atom_site_fract_z", "_atom_site_occupancy", "_atom_site_U_iso_or_equiv",
                   "_atom_site_disorder_assembly", "_atom_site_disorder_group"], rows)
    if edges:
        text += _loop(["_geom_bond_atom_site_label_1", "_geom_bond_atom_site_label_2", "_geom_bond_distance",
                       "_geom_bond_site_symmetry_1", "_geom_bond_site_symmetry_2"],
                      [[labels[a], labels[b], _num(math.dist(atoms[a].xyz, atoms[b].xyz)), ".", "."] for a, b in edges])
        # Also retain the same graph in the IUCr topology dictionary. Generic
        # links express coarse traces without claiming a chemical bond order.
        text += _loop(["_topol_node.id", "_topol_node.label"],
                      [[str(i+1), label] for i, label in enumerate(labels)])
        text += _loop(["_topol_atom.id", "_topol_atom.node_id", "_topol_atom.atom_label", "_topol_atom.element_symbol"],
                      [[str(i+1), str(i+1), label, atoms[i].element] for i, label in enumerate(labels)])
        text += _loop(["_topol_link.id", "_topol_link.node_id_1", "_topol_link.node_id_2", "_topol_link.distance", "_topol_link.type"],
                      [[str(i+1), str(a+1), str(b+1), _num(math.dist(atoms[a].xyz, atoms[b].xyz)), "gl"]
                       for i, (a, b) in enumerate(edges)])
    return text


def _mmcif(pdb: Pdb, atoms: list[Atom], edges, trace: bool) -> str:
    text = _cif_header(pdb, True)
    residues = {}
    for atom in atoms:
        if atom.residue_key not in residues:
            residues[atom.residue_key] = 1 + sum(key[0] == atom.segment for key in residues)
    if pdb.scale:
        for i, row in enumerate(pdb.scale, 1):
            for j, value in enumerate(row[:3], 1):
                text += f"_atom_sites.fract_transf_matrix[{i}][{j}] {_num(value)}\n"
            text += f"_atom_sites.fract_transf_vector[{i}] {_num(row[3])}\n"
        text += "\n"
    def identity(atom):
        return [f"S{atom.segment}", _quote(atom.residue), str(residues[atom.residue_key]), _quote(atom.name)]
    rows = [[a.raw[:6].strip() if a.raw else a.record_type, str(i+1), a.element, _quote(a.name), _quote(a.altloc) if a.altloc else ".",
             _quote(a.residue), f"S{a.segment}", "?", str(residues[a.residue_key]),
             _quote(a.insertion) if a.insertion else "?", *(_num(v) for v in a.xyz),
             _num(a.occupancy), _num(a.b_iso), _quote(a.seq), _quote(a.residue),
             _quote(a.chain) if a.chain else ".", _quote(a.name), "1"] for i, a in enumerate(atoms)]
    tags = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
            "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z",
            "occupancy", "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"]
    text += _loop(["_atom_site." + tag for tag in tags], rows)
    if edges:
        tags = ["id", "conn_type_id", "ptnr1_label_asym_id", "ptnr1_label_comp_id", "ptnr1_label_seq_id", "ptnr1_label_atom_id",
                "ptnr2_label_asym_id", "ptnr2_label_comp_id", "ptnr2_label_seq_id", "ptnr2_label_atom_id",
                "pdbx_ptnr1_label_alt_id", "pdbx_ptnr2_label_alt_id", "ptnr1_symmetry", "ptnr2_symmetry", "pdbx_dist_value", "details"]
        text += _loop(["_struct_conn." + tag for tag in tags],
                      [[str(i+1), "?", *identity(atoms[a]), *identity(atoms[b]),
                        atoms[a].altloc or ".", atoms[b].altloc or ".", "1_555", "1_555",
                        _num(math.dist(atoms[a].xyz, atoms[b].xyz)),
                        _quote("TopoStich virtual backbone link" if trace else "PDB CONECT; bond order unspecified")]
                       for i, (a, b) in enumerate(edges)])
    return text


def _pdb_field(value, width, precision):
    text = " " * width if value is None else f"{value:{width}.{precision}f}"
    if len(text) != width:
        raise ConversionError("Value exceeds PDB field width; choose CIF or mmCIF output")
    return text


def _pdb_atom_record(atom: Atom, serial: int) -> str:
    if (len(atom.name) > 4 or len(atom.residue) > 3 or len(atom.chain) > 1
            or len(atom.altloc) > 1 or len(atom.insertion) > 1
            or not re.fullmatch(r"-?\d+", atom.seq) or len(f"{int(atom.seq):4d}") != 4):
        raise ConversionError("Atom/residue identity cannot fit PDB fields; choose CIF or mmCIF output")
    name = (f" {atom.name:<3}" if len(atom.element) == 1 and len(atom.name) < 4
            and not atom.name[:1].isdigit() else f"{atom.name:<4}")
    return (f"{atom.record_type:<6}{serial:5d} {name}{atom.altloc:1}{atom.residue:>3} {atom.chain:1}"
            f"{int(atom.seq):4d}{atom.insertion:1}   "
            + "".join(_pdb_field(v, 8, 3) for v in atom.xyz)
            + _pdb_field(atom.occupancy, 6, 2) + _pdb_field(atom.b_iso, 6, 2)
            + f"          {atom.element:>2}")


def _pdb_headers(pdb: Pdb):
    if pdb.source_format == "pdb":
        return pdb.headers
    headers = []
    if pdb.cell:
        group = pdb.space_group
        # PDB uses H for rhombohedral symmetry expressed in hexagonal axes.
        if group.startswith('R') and pdb.cell[3:] == (90., 90., 120.) and ':R' not in group:
            group = 'H' + group[1:].removesuffix(':H')
        if len(group) > 11:
            raise ConversionError("Space-group name exceeds PDB field width")
        headers.append('CRYST1' + ''.join(_pdb_field(v, 9, 3) for v in pdb.cell[:3])
                       + ''.join(_pdb_field(v, 7, 2) for v in pdb.cell[3:]) + f' {group:<11}')
    if pdb.scale:
        for i, row in enumerate(pdb.scale, 1):
            headers.append(f'SCALE{i}    ' + ''.join(_pdb_field(v, 10, 6) for v in row[:3])
                           + '     ' + _pdb_field(row[3], 10, 5))
    headers.extend(f'REMARK 290 {i:7d}555   {op.upper()}' for i,op in pdb.symops.items())
    return headers


def _pdb(pdb: Pdb, atoms: list[Atom], edges) -> str:
    lines = ["REMARK 999 TOPOSTICH COORDINATE DERIVATIVE; LINKS ARE NOT RESTRAINTS"] + _pdb_headers(pdb)
    serials = []
    serial = 0
    for i, atom in enumerate(atoms):
        serial += 1
        if serial > 99999:
            raise ConversionError("PDB serial limit exceeded; choose CIF or mmCIF")
        serials.append(serial)
        lines.append(atom.raw[:6] + f"{serial:5d}" + atom.raw[11:] if atom.raw
                     else _pdb_atom_record(atom, serial))
        if i == len(atoms)-1 or atom.segment != atoms[i+1].segment:
            serial += 1
            if serial > 99999:
                raise ConversionError("PDB serial limit exceeded; choose CIF or mmCIF")
            lines.append(f"TER   {serial:5d}      {atom.residue:>3} {atom.chain:1}{atom.seq:>4}{atom.insertion:1}")
    adjacency = {}
    for a, b in edges:
        adjacency.setdefault(serials[a], []).append(serials[b])
        adjacency.setdefault(serials[b], []).append(serials[a])
    for atom, neighbours in sorted(adjacency.items()):
        for offset in range(0, len(neighbours), 4):
            lines.append(f"CONECT{atom:5d}" + "".join(f"{n:5d}" for n in sorted(neighbours)[offset:offset+4]))
    return "\n".join(lines + ["END", ""])


def read_structure(text: str, model: int = 1, input_format: str = "auto") -> Pdb:
    """Read a coordinate snapshot by content, never by filename extension alone."""
    fmt = input_format.lower()
    if fmt not in ("auto", "pdb", "mmcif"):
        raise ConversionError(f"Unknown input format: {input_format}")
    text = text.lstrip('\ufeff')
    first = next((line.strip() for line in text.splitlines() if line.strip() and not line.lstrip().startswith('#')), '')
    if fmt == "auto":
        fmt = "mmcif" if first.lower().startswith('data_') else "pdb"
    if fmt == "mmcif":
        from .mmcif_input import read_mmcif
        return read_mmcif(text, model)
    return read_pdb(text, model)


def convert_file(source: Path | str, output: Path | str, *, atoms: Sequence[str] | None = None,
                 output_format: str = "cif", model: int = 1, altloc: str | None = None,
                 input_format: str = "auto") -> ConversionResult:
    source, output = Path(source), Path(output)
    fmt = output_format.lower()
    if fmt not in FORMATS:
        raise ConversionError(f"Unknown output format: {output_format}")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Output already exists: {output}")
    raw = source.read_bytes()
    pdb = read_structure(raw.decode("utf-8-sig"), model, input_format)
    selected, edges, order, missing = _select(pdb, atoms, altloc)
    source_hash = hashlib.sha256(raw).hexdigest()
    warnings = []
    if not order:
        warnings.extend(pdb.connection_warnings)
    if pdb.model_count > 1:
        warnings.append(f"Selected model {model} of {pdb.model_count}; other models excluded")
    if missing:
        warnings.append(f"Skipped {len(missing)} missing selected atoms; bridged within each chain/TER segment")
    if not pdb.symops and pdb.space_group.replace(" ", "") not in ("", "P1"):
        warnings.append("No fractional symmetry operators in input; space-group symbol retained for the reader")
    if fmt == "pdb":
        if pdb.source_format == "mmcif":
            warnings.append("PDB output rounds coordinates to 0.001 A and occupancy/B values to 0.01; use CIF/mmCIF to retain precision")
        text = _pdb(pdb, selected, edges)
        text = f"REMARK 999 SOURCE SHA256 {source_hash}\n" + text
    else:
        text = _core_cif(pdb, selected, edges, bool(order)) if fmt == "cif" else _mmcif(pdb, selected, edges, bool(order))
        text = (f"# TopoStich {__version__}; source SHA256 {source_hash}\n"
                f"# Model {model}; selection {','.join(order) or 'all atoms'}; no inferred stitching\n"
                "# Coordinate/graph derivative, not a refinement or deposition dictionary.\n" + text)
    encoded = text.encode("ascii")
    # Exclusive creation prevents overwriting an input, historical result, or symlink.
    with output.open("xb") as stream:
        stream.write(encoded)
    return ConversionResult(1, __version__, source.name, source_hash, output.name,
                            hashlib.sha256(encoded).hexdigest(), fmt, model, order, len(selected), len(edges),
                            len({a.segment for a in selected}), len(missing), missing, warnings)
