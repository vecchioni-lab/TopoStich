"""Residue-aware nucleotide atom roles, separate from graph construction.

The NARestraints workbook maps canonical base roles to actual names. It is not
a bond dictionary, a universal blacklist, or a list of atoms guaranteed present.
"""
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
import json
from typing import Iterable

SUGAR_ATOMS = frozenset(("C1'", "C2'", "C3'", "C4'", "C5'", "O2'", "O3'", "O4'", "O5'"))
PHOSPHATE_ATOMS = frozenset(('P', 'OP1', 'OP2', 'OP3'))
PHOSPHATE_ALIASES = {'O1P': 'OP1', 'O2P': 'OP2', 'O3P': 'OP3'}

def normalize_atom_name(name: str) -> str:
    name = name.strip().replace('*', "'").replace('′', "'")
    return PHOSPHATE_ALIASES.get(name, name)

@lru_cache(maxsize=1)
def _catalog():
    return json.loads(files('topostich').joinpath('data/nucleotide_atoms.json').read_text())

@dataclass(frozen=True)
class ResidueDefinition:
    code: str
    source_sheet: str
    source_row: int
    base_analog: str
    sugar_type: str
    positions: tuple[tuple[str, str], ...]
    substituents: tuple[tuple[str, str], ...]

    @property
    def atom_names(self) -> frozenset[str]:
        return frozenset(value for _, value in self.positions + self.substituents)

@lru_cache(maxsize=1)
def _definitions():
    result = {}
    for row in _catalog()['records']:
        definition = ResidueDefinition(**{k:v for k,v in row.items() if k not in ('positions', 'substituents')},
            positions=tuple((k, normalize_atom_name(v)) for k,v in row['positions'].items()),
            substituents=tuple((k, normalize_atom_name(v)) for k,v in row['substituents'].items()))
        result.setdefault(definition.code, []).append(definition)
    return {k:tuple(v) for k,v in result.items()}

def lookup_residue(code: str) -> tuple[ResidueDefinition, ...]:
    """Return all source rows; no canonical-base substitution or first-row guess."""
    return _definitions().get(code.strip().upper(), ())

def header_base_atoms() -> frozenset[str]:
    """Canonical base-role headers, useful for selector diagnostics only."""
    return frozenset(_catalog()['header_base_atoms'])

def _matching_definitions(code, observed):
    records = lookup_residue(code)
    if not records:
        return (), 'unmapped'
    scores = [len(record.atom_names.intersection(observed)) for record in records]
    best = max(scores)
    selected = tuple(record for record, score in zip(records, scores) if score == best)
    signatures = {(record.positions, record.substituents) for record in selected}
    if len(signatures) > 1:
        status = 'ambiguous'
    elif len(records) == 1:
        status = 'unique'
    elif len(selected) == 1:
        status = 'matched'
    else:
        status = 'equivalent-duplicates'
    return selected, status

@dataclass(frozen=True)
class AtomRecognition:
    residue_code: str
    normalized_name: str
    role: str
    canonical_positions: tuple[str, ...]
    candidate_roles: tuple[str, ...]
    mapping_status: str

def classify_atom(residue_code: str, atom_name: str, *, observed_atoms: Iterable[str] = ()) -> AtomRecognition:
    """Classify an actual name using its residue row and conventional sugar names.

    Unknown rows/names remain unknown. Competing mappings or sugar/base name
    collisions remain ambiguous. No coordinates or chemical bonds are invented.
    """
    name = normalize_atom_name(atom_name)
    observed = {normalize_atom_name(v) for v in observed_atoms} or {name}
    records, status = _matching_definitions(residue_code, observed)
    roles = set()
    if name in SUGAR_ATOMS:
        roles.add('sugar')
    if name in PHOSPHATE_ATOMS:
        roles.add('phosphate')
    positions = {key for record in records for key,value in record.positions if value == name}
    mapped = [name in record.atom_names for record in records]
    if any(mapped):
        roles.add('base')
    conflict = any(mapped) and not all(mapped)
    role = ('ambiguous' if len(roles) > 1 or conflict else next(iter(roles)) if roles else 'unknown')
    return AtomRecognition(residue_code.strip().upper(), name, role, tuple(sorted(positions)),
                           tuple(sorted(roles)), status)

@dataclass(frozen=True)
class ResidueRecognition:
    residue_code: str
    nucleotide_like: bool
    mapping_status: str
    definitions: tuple[ResidueDefinition, ...]
    atoms: tuple[AtomRecognition, ...]
    present_sugar_atoms: tuple[str, ...]
    present_phosphate_atoms: tuple[str, ...]
    present_base_atoms: tuple[str, ...]
    unknown_atoms: tuple[str, ...]
    ambiguous_atoms: tuple[str, ...]
    missing_mapped_base_atoms: tuple[str, ...]

def recognize_residue(residue_code: str, atom_names: Iterable[str]) -> ResidueRecognition:
    """Recognize one observed residue; call per chain/residue/insertion identity."""
    observed = tuple(dict.fromkeys(normalize_atom_name(name) for name in atom_names))
    records, status = _matching_definitions(residue_code, set(observed))
    atoms = tuple(classify_atom(residue_code, name, observed_atoms=observed) for name in observed)
    def names(role):
        return tuple(a.normalized_name for a in atoms if a.role == role)
    # Report only unanimously mapped missing names when competing rows remain.
    common = set.intersection(*(set(r.atom_names) for r in records)) if records else set()
    return ResidueRecognition(residue_code.strip().upper(), bool(records or SUGAR_ATOMS.intersection(observed)),
        status, records, atoms, names('sugar'), names('phosphate'), names('base'), names('unknown'),
        names('ambiguous'), tuple(sorted(common.difference(observed))))
