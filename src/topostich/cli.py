"""Installed topostich command and machine-readable conversion receipts."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

from . import __version__
from .converter import BACKBONE_ORDER, ConversionError, convert_file


def _normalize_all_option(argv):
    """Let argparse accept the requested dash-prefixed -atom -all value."""
    tokens = list(sys.argv[1:] if argv is None else argv)
    index = 0
    while index < len(tokens) - 1:
        if tokens[index] == "--":
            break
        if tokens[index] in ("-atom", "--atom") and tokens[index + 1] == "-all":
            tokens[index:index + 2] = [tokens[index] + "=-all"]
        index += 1
    return tokens


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="topostich", description="Convert PDB/mmCIF coordinates to core CIF, with optional backbone traces.")
    parser.add_argument("input", type=Path, help="input PDB or PDBx/mmCIF path")
    parser.add_argument("--input-format", choices=("auto", "pdb", "mmcif"), default="auto",
                        help="input syntax (default: detect from content)")
    parser.add_argument("-atom", "--atom", action="append", help="use -atom -all for all-atom conversion, or repeat backbone atom names "
                        "and quote primes, e.g. -atom \"O3'\". "
                        "Order is always " + " -> ".join(BACKBONE_ORDER))
    parser.add_argument("-out", "--out", required=True, type=Path, help="new output file (never overwritten)")
    formats = parser.add_mutually_exclusive_group()
    formats.add_argument("--format", type=str.lower, choices=("cif", "pdb", "mmcif"), dest="format")
    for flag, value in (("-cif", "cif"), ("-pdb", "pdb"), ("-mmCif", "mmcif"), ("-mmcif", "mmcif")):
        formats.add_argument(flag, action="store_const", const=value, dest="format")
    parser.set_defaults(format="cif")
    parser.add_argument("--model", type=int, default=1, help="one-based model in file order (default: first)")
    parser.add_argument("--altloc", help="retain blank alternate plus this alternate ID; required for an ambiguous trace")
    parser.add_argument("--json", action="store_true", help="print schema-1 receipt with counts, SHA-256 and missing atoms")
    parser.add_argument("--version", action="version", version=f"TopoStich {__version__}")
    args = parser.parse_args(_normalize_all_option(argv))
    all_values = {"all", "-all"}
    if args.atom and any(name in all_values for name in args.atom):
        if any(name not in all_values for name in args.atom):
            parser.error("-atom -all cannot be combined with named atom selections")
        args.atom = None
    try:
        result = convert_file(args.input, args.out, atoms=args.atom, output_format=args.format,
                              model=args.model, altloc=args.altloc, input_format=args.input_format)
    except (ConversionError, OSError, UnicodeError) as exc:
        print(f"topostich: {exc}", file=sys.stderr)
        return 2
    if args.json:
        print(json.dumps(asdict(result), indent=2))
    else:
        print(f"Wrote {args.out}: {result.atom_count} atoms, {result.bond_count} links, "
              f"{result.segment_count} chain/TER segments ({result.output_format}).")
        for warning in result.warnings:
            print(f"Note: {warning}", file=sys.stderr)
    return 0
