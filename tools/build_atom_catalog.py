"""Extract an auditable atom-role lookup; the input workbook is never modified.

Development-only dependency: openpyxl. Runtime uses the packaged JSON snapshot.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

SHEETS = ('adenine', 'thymine', 'guanine', 'cytosine', 'D', 'B', 'S', 'Z', 'P', 'I', 'X', 'K', 'unique')
ATOM = re.compile(r"(?:C|N|O|S|P|F|CL|BR|I|SE|TE)[A-Z0-9']{0,3}$")
POSITION = re.compile(r'[CNOS][0-9]+$')

def build(workbook):
    import openpyxl
    raw = Path(workbook).read_bytes()
    book = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    records, excluded, headers = [], [], set()
    try:
        for name in SHEETS:
            sheet = book[name]
            rows = sheet.iter_rows(values_only=True)
            columns = list(next(rows))
            start = columns.index('Function') + 1
            for number, row in enumerate(rows, 2):
                if row[0] is None:
                    continue
                record = {'code': str(row[0]).strip().upper(), 'source_sheet': name,
                    'source_row': number, 'base_analog': str(row[3] or '').strip(),
                    'sugar_type': str(row[5] or '').strip(), 'positions': {}, 'substituents': {}}
                for index in range(start, len(columns)):
                    header = columns[index]
                    if not header:
                        continue
                    if POSITION.fullmatch(header):
                        headers.add(header)
                    value = row[index]
                    if value is None:
                        continue
                    value = str(value).strip().replace('*', "'").replace('′', "'")
                    if not ATOM.fullmatch(value):
                        excluded.append({'sheet': name, 'row': number, 'column': index+1, 'header': header})
                        continue
                    key = 'positions' if POSITION.fullmatch(header) else 'substituents'
                    record[key][header] = value
                records.append(record)
    finally:
        book.close()
    counts = Counter(r['code'] for r in records)
    return {'schema_version': 1,
        'source': {'repository': 'NARestraints', 'path': 'restraints/data/Ligands.xlsx',
                   'sha256': hashlib.sha256(raw).hexdigest(), 'sheets': list(SHEETS)},
        'semantics': 'Residue code -> canonical role -> actual atom name. Only observed atoms are present.',
        'normalization': 'Trim mapped names; normalize star/Unicode prime. Preserve row alternatives and role collisions.',
        'header_base_atoms': sorted(headers), 'record_count': len(records), 'residue_count': len(counts),
        'duplicate_codes': sorted(k for k,v in counts.items() if v > 1),
        'excluded_non_atom_cells': excluded, 'records': records}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workbook', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    catalog = build(args.workbook)
    with args.output.open('x') as stream:
        json.dump(catalog, stream, indent=2)
        stream.write('\n')
    print(f"Wrote {catalog['residue_count']} residue codes / {catalog['record_count']} rows")

if __name__ == '__main__':
    main()
