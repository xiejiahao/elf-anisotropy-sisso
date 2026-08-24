from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from .analysis import analyze_elfcar
from .errors import DescriptorError
from .output import format_json, format_text, write_outputs
from .profiles import DEFAULT_SAMPLES


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="elf-anisotropy",
        description="Calculate static ELF anisotropy for six-coordinate Sn/Pb iodides.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    analyze = subparsers.add_parser(
        "analyze",
        help="calculate A_ELF and A_ELF^2 from one ELFCAR",
    )
    analyze.add_argument("elfcar", help="path to a VASP ELFCAR")
    analyze.add_argument(
        "--output",
        metavar="DIRECTORY",
        help="write summary.json, per_site.csv, and per_bond.csv",
    )
    analyze.add_argument(
        "--format",
        choices=("text", "json"),
        default="text",
        help="standard-output format (default: text)",
    )
    analyze.add_argument(
        "--samples",
        type=int,
        default=DEFAULT_SAMPLES,
        help=f"points sampled along each metal-I direction (default: {DEFAULT_SAMPLES})",
    )
    return parser


def run(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        result = analyze_elfcar(args.elfcar, samples=args.samples)
        if args.output:
            write_outputs(result, args.output)
        print(format_json(result) if args.format == "json" else format_text(result))
    except (DescriptorError, OSError, ValueError) as exc:
        print(f"elf-anisotropy: error: {exc}", file=sys.stderr)
        return 2
    return 0


def main() -> None:
    raise SystemExit(run())


if __name__ == "__main__":
    main()
