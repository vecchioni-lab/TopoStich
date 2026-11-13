"""Native PDBx/mmCIF input, retaining row order, precision and chain identities."""
import re

from .converter import Atom, Pdb, ConversionError, _number


def read_mmcif(text: str, model: int = 1) -> Pdb:
    try:
        import gemmi
    except ImportError as exc:
        raise ConversionError('mmCIF input requires Gemmi; reinstall TopoTrace with its declared dependencies') from exc
    try:
        document = gemmi.cif.read_string(text)
    except (ValueError, RuntimeError) as exc:
        raise ConversionError(f'Invalid mmCIF syntax: {exc}') from exc
    if len(document) != 1:
        raise ConversionError('Input must contain exactly one mmCIF data block')
    block = document.sole_block()

    def clean(value):
        return '' if value in (None, '.', '?') else gemmi.cif.as_string(value)

    def scalar(tag):
        return clean(block.find_value(tag))

    def category(prefix, preserve_unknown=False):
        table = block.find_mmcif_category(prefix)
        names = [tag[len(prefix):].lower() for tag in table.tags]
        return [dict(zip(names, ('?' if preserve_unknown and value == '?' else clean(value) for value in row)))
                for row in table]

    rows = category('_atom_site.')
    if not rows or any(not all(key in row for key in ('cartn_x','cartn_y','cartn_z')) for row in rows):
        raise ConversionError('Expected macromolecular mmCIF _atom_site.cartn_x/y/z coordinates; core CIF/CCD input is not supported')
    model_ids = list(dict.fromkeys(row.get('pdbx_pdb_model_num') or '1' for row in rows))
    if not 1 <= model <= len(model_ids):
        raise ConversionError(f'Model must be in 1..{len(model_ids)} (one-based file order)')
    model_id = model_ids[model-1]
    chosen = [row for row in rows if (row.get('pdbx_pdb_model_num') or '1') == model_id]
    atoms, atom_rows = [], []
    previous = None
    segment = 0
    for row in chosen:
        def identity(name):
            return row.get('auth_'+name) or row.get('label_'+name) or ''
        name, residue, seq, chain = (identity(k) for k in ('atom_id','comp_id','seq_id','asym_id'))
        if not name or not residue or not seq:
            raise ConversionError('Missing mmCIF atom/residue identity (auth and label fields)')
        label = row.get('label_asym_id') or chain
        boundary = label, chain
        if boundary != previous:
            segment += 1
        previous = boundary
        record = row.get('group_pdb') or 'ATOM'
        if record not in ('ATOM', 'HETATM'):
            raise ConversionError(f'Invalid mmCIF group_pdb: {record!r}')
        element = row.get('type_symbol')
        if not element or not re.fullmatch(r'[A-Za-z]{1,2}', element):
            raise ConversionError('Missing or invalid mmCIF atom element')
        serial = row.get('id') or str(len(atoms)+1)
        coords = tuple(_number(row[k], f'mmCIF {k}') for k in ('cartn_x','cartn_y','cartn_z'))
        atoms.append(Atom(serial, name, row.get('label_alt_id',''), residue, chain, seq,
            row.get('pdbx_pdb_ins_code',''), segment, coords,
            _number(row.get('occupancy',''), 'occupancy', True),
            _number(row.get('b_iso_or_equiv',''), 'B factor', True), element.capitalize(), '',
            label, row.get('label_seq_id',''), row.get('label_comp_id',''), row.get('label_atom_id',''), record))
        atom_rows.append(row)
    values = [scalar('_cell.'+key) for key in ('length_a','length_b','length_c','angle_alpha','angle_beta','angle_gamma')]
    if any(values) and not all(values):
        raise ConversionError('Incomplete mmCIF cell: all six parameters are required')
    cell = tuple(_number(v,'mmCIF cell') for v in values) if all(values) else None
    group = scalar('_space_group.name_H-M_alt') or scalar('_symmetry.space_group_name_H-M')
    if not group and scalar('_space_group.IT_number'):
        try:
            group = gemmi.find_spacegroup_by_number(int(scalar('_space_group.IT_number'))).hm
        except (ValueError, AttributeError) as exc:
            raise ConversionError('Invalid mmCIF space-group number') from exc
    scale_values = [[scalar(f'_atom_sites.fract_transf_matrix[{i}][{j}]') for j in (1,2,3)]
                    + [scalar(f'_atom_sites.fract_transf_vector[{i}]')] for i in (1,2,3)]
    flat = [v for row in scale_values for v in row]
    if any(flat) and not all(flat):
        raise ConversionError('Incomplete mmCIF fractional transformation: matrix and vector are required')
    scale = tuple(tuple(_number(v,'mmCIF fractional transform') for v in row) for row in scale_values) if all(flat) else None
    symops = {}
    operations = category('_space_group_symop.') or category('_symmetry_equiv.')
    for ordinal, row in enumerate(operations, 1):
        op = row.get('operation_xyz') or row.get('pos_as_xyz')
        if not op:
            raise ConversionError('Missing mmCIF symmetry operator')
        try:
            index = int(row.get('id') or ordinal)
            gemmi.Op(op)
        except (ValueError, RuntimeError) as exc:
            raise ConversionError('Invalid mmCIF symmetry operator/id') from exc
        if index <= 0 or (index in symops and symops[index] != op):
            raise ConversionError('Conflicting mmCIF symmetry operator IDs')
        symops[index] = op.replace(' ','').lower()
    connections, warnings = [], []
    for row in category('_struct_conn.', preserve_unknown=True):
        # Do not invent an image for a connection involving crystallographic symmetry.
        identity_images = True
        for partner in (1,2):
            image = row.get(f'ptnr{partner}_symmetry','')
            if image in ('', '.'):
                continue
            match = re.fullmatch(r'(\d+)_555',image)
            index = int(match[1]) if match else None
            op = symops.get(index) if index is not None else None
            if not match or not ((op and gemmi.Op(op) == gemmi.Op('x,y,z')) or (index == 1 and not symops)):
                identity_images = False
        if not identity_images:
            warnings.append('Omitted a symmetry-image _struct_conn link; image/stitching links are not supported')
            continue
        endpoints = []
        out_of_model = False
        for partner in (1,2):
            prefix = f'ptnr{partner}_'
            partner_model = row.get(f'pdbx_ptnr{partner}_pdb_model_num')
            if partner_model == '?':
                partner_model = ''
            if partner_model and partner_model != model_id:
                out_of_model = True
                break
            criteria = [(key, row[prefix+key]) for key in
                        ('label_asym_id','label_seq_id','label_comp_id','label_atom_id',
                         'auth_asym_id','auth_seq_id','auth_comp_id','auth_atom_id') if row.get(prefix+key) not in (None,'','?')]
            if not any(key.endswith('atom_id') for key,_ in criteria):
                raise ConversionError('mmCIF connection endpoint lacks an atom identifier')
            alt = row.get(f'pdbx_ptnr{partner}_label_alt_id','')
            insertion = row.get(f'pdbx_ptnr{partner}_pdb_ins_code','')
            alt = '' if alt == '?' else alt
            insertion = '' if insertion == '?' else insertion
            candidates = [a for a, atom_row in zip(atoms,atom_rows)
                          if all(atom_row.get(key) == value for key,value in criteria)
                          and (not alt or a.altloc == alt) and (not insertion or a.insertion == insertion)]
            if len(candidates) > 1:
                raise ConversionError('Ambiguous mmCIF connection endpoint; alternate/residue identity is incomplete')
            if not candidates:
                warnings.append('Omitted an unresolved _struct_conn endpoint in the selected model')
                out_of_model = True
                break
            endpoints.append(candidates[0].serial)
        if not out_of_model:
            connections.append(tuple(endpoints))
    # PDB-derived mmCIF often carries a six-decimal SCALE approximation. Full
    # precision native matrices and every nonzero origin shift must stay exact.
    rounded = scale is None or (all(row[3] == 0 for row in scale) and
        all(abs(v-round(v,6)) < 1e-14 for row in scale for v in row[:3]))
    return Pdb(atoms, cell, group, scale, symops, connections, [], len(model_ids), 'mmcif', warnings, rounded)
