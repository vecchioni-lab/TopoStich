"""Audit a TopoStich CGD against its receipt using Systre's own parser.

    python tools/check_stitched_cgd.py graph.cgd graph.json --systre Systre-19.6.0.jar

Compiles tools/CheckCgd.java into a temporary directory, runs it in --stitched
mode (Systre parse + declared degrees), then compares Systre's primitive-cell
edges with the receipt (schema 1 or 2):

* every Systre edge joins a node pair listed in the receipt, and vice versa;
* for each node pair, Systre's edge count equals the expected orbit count:
  one edge per primitive-cell copy, or half that for an edge reversed by its
  own image (same node, involutive image operation);
* periodic dimensions of Systre's connected components are reported.

Needs a JDK (javac/java) and a separately obtained Systre jar; neither is a
TopoStich runtime dependency. Exit status 0 = pass, 1 = mismatch, 2 = usage.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def systre_records(cgd, jar, java='java', javac='javac'):
    with tempfile.TemporaryDirectory(prefix='checkcgd-') as build:
        subprocess.run([javac, '-cp', str(jar), '-d', build, str(HERE/'CheckCgd.java')], check=True)
        run = subprocess.run([java, '-cp', f'{jar}{":" if sys.platform != "win32" else ";"}{build}',
                              'CheckCgd', '--stitched', str(cgd)], text=True, capture_output=True)
    if run.returncode:
        raise RuntimeError(f'Systre audit failed:\n{run.stderr.strip()}')
    records = [line.split('\t') for line in run.stdout.splitlines() if line]
    return records


def expected_pairs(receipt):
    """Unordered node pair -> number of edge orbits whose image is self-reversing or not."""
    import gemmi
    identity = gemmi.Op('x,y,z').triplet()
    pairs = Counter()
    half = Counter()
    for left, right in receipt['edges']:
        pairs[frozenset((left, right))] += 1
    for edge in receipt.get('stitched_edges', ()):
        key = frozenset((edge['source'], edge['target']))
        op = gemmi.Op(edge['image_operation'])
        if edge['source'] == edge['target'] and op.combine(op).triplet() == identity:
            half[key] += 1
        else:
            pairs[key] += 1
    return pairs, half


def compare(records, receipt):
    problems = []
    nodes = next(int(r[1]) for r in records if r[0] == 'NODES')
    if nodes % receipt['node_count']:
        problems.append(f'Systre node count {nodes} is not a multiple of {receipt["node_count"]} receipt nodes')
        return problems, {}
    copies = nodes // receipt['node_count']
    full, half = expected_pairs(receipt)
    expected = Counter()
    for key, n in full.items():
        expected[key] += n*copies
    for key, n in half.items():
        if copies % 2:
            problems.append(f'self-reversing edge {sorted(key)} with odd primitive multiplicity {copies}')
        expected[key] += n*copies//2
    observed = Counter(frozenset(r[1:3]) for r in records if r[0] == 'EDGE')
    for key in sorted(set(expected) | set(observed), key=sorted):
        if expected[key] != observed[key]:
            problems.append(f'{"-".join(sorted(key))}: Systre {observed[key]} edges, receipt implies {expected[key]}')
    components = Counter((int(r[1]), int(r[2])) for r in records if r[0] == 'COMPONENT')
    summary = dict(systre_nodes=nodes, primitive_copies=copies,
                   systre_edges=sum(observed.values()), expected_edges=sum(expected.values()),
                   trace_edges=len(receipt['edges']), stitched_edges=len(receipt.get('stitched_edges', ())),
                   components={f'dimension {d}, {n} nodes': c for (d, n), c in sorted(components.items())},
                   warnings=[r[1] for r in records if r[0] == 'WARNING'])
    return problems, summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cgd', type=Path)
    parser.add_argument('receipt', type=Path)
    parser.add_argument('--systre', type=Path, required=True, help='Systre jar (e.g. Systre-19.6.0.jar)')
    parser.add_argument('--java', default='java')
    parser.add_argument('--javac', default='javac')
    args = parser.parse_args(argv)
    receipt = json.loads(args.receipt.read_text())
    try:
        records = systre_records(args.cgd, args.systre.resolve(), args.java, args.javac)
    except (RuntimeError, subprocess.CalledProcessError, OSError) as error:
        print(f'FAIL: {error}', file=sys.stderr)
        return 1
    problems, summary = compare(records, receipt)
    print(json.dumps(summary, indent=2))
    if problems:
        for problem in problems:
            print('FAIL: '+problem, file=sys.stderr)
        return 1
    print(f'PASS: Systre adjacency matches {args.receipt.name}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
