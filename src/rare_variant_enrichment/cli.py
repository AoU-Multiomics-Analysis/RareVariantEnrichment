import argparse
import math
from pathlib import Path

from rare_variant_enrichment.phenotypes import prepare_phenotypes


COMMANDS = ("prepare-phenotypes", "classify-chromosome", "gather", "calculate")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rare-variant-enrichment")
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare-phenotypes")
    prepare_parser.add_argument("--phenotype-bed", required=True, type=Path)
    prepare_parser.add_argument("--vcf-samples", required=True, type=Path)
    prepare_parser.add_argument("--chromosomes", required=True, type=parse_csv_strings)
    prepare_parser.add_argument("--z-thresholds", required=True, type=parse_csv_floats)
    prepare_parser.add_argument("--tail", required=True, choices=("absolute", "positive", "negative"))
    prepare_parser.add_argument("--feature-output", required=True, type=Path)
    prepare_parser.add_argument("--sample-output", required=True, type=Path)
    prepare_parser.add_argument("--qc-output", required=True, type=Path)
    for command in COMMANDS[1:]:
        subparsers.add_parser(command)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "prepare-phenotypes":
        prepare_phenotypes(
            args.phenotype_bed,
            args.vcf_samples,
            args.chromosomes,
            args.z_thresholds,
            args.tail,
            args.feature_output,
            args.sample_output,
            args.qc_output,
        )
    return 0


def parse_csv_strings(value: str) -> list[str]:
    values = [item.strip() for item in value.split(",")]
    if not values or any(not item for item in values):
        raise argparse.ArgumentTypeError("Expected a non-empty comma-separated list")
    return values


def parse_csv_ints(value: str) -> list[int]:
    try:
        return [int(item) for item in parse_csv_strings(value)]
    except ValueError as error:
        raise argparse.ArgumentTypeError("Expected comma-separated integers") from error


def parse_csv_floats(value: str) -> list[float]:
    try:
        values = [float(item) for item in parse_csv_strings(value)]
    except ValueError as error:
        raise argparse.ArgumentTypeError("Expected comma-separated numbers") from error
    if not all(math.isfinite(item) for item in values):
        raise argparse.ArgumentTypeError("Expected finite comma-separated numbers")
    return values


if __name__ == "__main__":
    raise SystemExit(main())
