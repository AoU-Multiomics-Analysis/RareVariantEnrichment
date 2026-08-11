import argparse
import math
from pathlib import Path

from rare_variant_enrichment.aggregation import gather_outputs
from rare_variant_enrichment.phenotypes import prepare_phenotypes
from rare_variant_enrichment.statistics import calculate_enrichment
from rare_variant_enrichment.vat import prepare_vat
from rare_variant_enrichment.variants import classify_chromosome


COMMANDS = ("prepare-phenotypes", "prepare-vat", "classify-chromosome", "gather", "calculate")


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
    vat_parser = subparsers.add_parser("prepare-vat")
    vat_parser.add_argument("--vat", required=True, type=Path)
    vat_parser.add_argument("--chromosomes", required=True, type=parse_csv_strings)
    vat_parser.add_argument("--schema-output", required=True, type=Path)
    vat_parser.add_argument("--loftee-enabled-output", required=True, type=Path)
    classify_parser = subparsers.add_parser("classify-chromosome")
    classify_parser.add_argument("--vcf", required=True, type=Path)
    classify_parser.add_argument("--vat", required=True, type=Path)
    classify_parser.add_argument("--vat-schema", required=True, type=Path)
    classify_parser.add_argument("--features", required=True, type=Path)
    classify_parser.add_argument("--shared-samples", required=True, type=Path)
    classify_parser.add_argument("--chromosome", required=True)
    classify_parser.add_argument("--exact-ac", required=True, type=parse_csv_ints)
    classify_parser.add_argument("--cumulative-ac-max", required=True, type=parse_csv_ints)
    classify_parser.add_argument("--consequence-classes", required=True, type=parse_csv_strings)
    classify_parser.add_argument("--maximum-gvs-maf", required=True, type=float)
    classify_parser.add_argument("--max-distance", required=True, type=int)
    classify_parser.add_argument("--annotation-chunk-size-bp", required=True, type=int)
    classify_parser.add_argument("--carrier-output", required=True, type=Path)
    classify_parser.add_argument("--regions-output", required=True, type=Path)
    classify_parser.add_argument("--qc-output", required=True, type=Path)
    gather_parser = subparsers.add_parser("gather")
    gather_parser.add_argument("--carrier-input", action="append", required=True, type=Path)
    gather_parser.add_argument("--qc-input", action="append", required=True, type=Path)
    gather_parser.add_argument("--carrier-output", required=True, type=Path)
    gather_parser.add_argument("--qc-output", required=True, type=Path)
    calculate_parser = subparsers.add_parser("calculate")
    calculate_parser.add_argument("--phenotype-bed", required=True, type=Path)
    calculate_parser.add_argument("--shared-samples", required=True, type=Path)
    calculate_parser.add_argument("--carriers", required=True, type=Path)
    calculate_parser.add_argument("--features", required=True, type=Path)
    calculate_parser.add_argument("--exact-ac", required=True, type=parse_csv_ints)
    calculate_parser.add_argument("--cumulative-ac-max", required=True, type=parse_csv_ints)
    calculate_parser.add_argument("--z-thresholds", required=True, type=parse_csv_floats)
    calculate_parser.add_argument("--distance-thresholds", required=True, type=parse_csv_ints)
    calculate_parser.add_argument(
        "--tail", required=True, choices=("absolute", "positive", "negative")
    )
    calculate_parser.add_argument("--output-tsv", required=True, type=Path)
    calculate_parser.add_argument("--output-json", required=True, type=Path)
    calculate_parser.add_argument("--phenotype-qc", type=Path)
    calculate_parser.add_argument("--chromosome-qc", type=Path)
    calculate_parser.add_argument("--selected-chromosomes", type=parse_csv_strings)
    calculate_parser.add_argument("--container-image")
    calculate_parser.add_argument("--workflow-version", default="unknown")
    calculate_parser.add_argument("--max-retries", type=int, default=0)
    calculate_parser.add_argument(
        "--index-provenance",
        choices=("generated", "supplied", "unknown"),
        default="unknown",
    )
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
    elif args.command == "prepare-vat":
        schema = prepare_vat(args.vat, args.chromosomes, args.schema_output)
        args.loftee_enabled_output.write_text("true\n" if schema.lof is not None else "false\n")
    elif args.command == "classify-chromosome":
        classify_chromosome(
            args.vcf,
            args.vat,
            args.vat_schema,
            args.features,
            args.shared_samples,
            args.chromosome,
            args.exact_ac,
            args.cumulative_ac_max,
            args.consequence_classes,
            args.maximum_gvs_maf,
            args.max_distance,
            args.annotation_chunk_size_bp,
            args.carrier_output,
            args.regions_output,
            args.qc_output,
        )
    elif args.command == "gather":
        gather_outputs(args.carrier_input, args.qc_input, args.carrier_output, args.qc_output)
    elif args.command == "calculate":
        calculate_enrichment(
            args.phenotype_bed,
            args.shared_samples,
            args.carriers,
            args.features,
            args.exact_ac,
            args.cumulative_ac_max,
            args.z_thresholds,
            args.distance_thresholds,
            args.tail,
            args.output_tsv,
            args.output_json,
            phenotype_qc_path=args.phenotype_qc,
            chromosome_qc_path=args.chromosome_qc,
            selected_chromosomes=args.selected_chromosomes,
            container_image=args.container_image,
            workflow_version=args.workflow_version,
            max_retries=args.max_retries,
            index_provenance=args.index_provenance,
        )
    return 0


def parse_csv_strings(value: str) -> list[str]:
    values = [item.strip() for item in value.split(",")]
    if not values or any(not item for item in values):
        raise argparse.ArgumentTypeError("Expected a non-empty comma-separated list")
    return values


def parse_csv_ints(value: str) -> list[int]:
    if value == "":
        return []
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
