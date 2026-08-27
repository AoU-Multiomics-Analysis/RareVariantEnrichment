import argparse
import logging
import math
from pathlib import Path

from rare_variant_enrichment.aggregation import gather_outputs
from rare_variant_enrichment.carrier_aggregation import gather_variant_carriers
from rare_variant_enrichment.carrier_definitions import build_carrier_definitions
from rare_variant_enrichment.carrier_extraction import (
    extract_chromosome_carriers,
    prepare_carrier_inputs,
)
from rare_variant_enrichment.io import read_nonempty_lines, write_json
from rare_variant_enrichment.lof_pc import (
    build_pc_chunks,
    calculate_carrier_pc_enrichment,
    calculate_lof_pc_enrichment,
    merge_lof_pc_enrichment,
    prepare_protein_coding_genes,
    read_principal_component_header,
)
from rare_variant_enrichment.pc_selection import (
    DEFAULT_PLATEAU_FRACTION,
    DEFAULT_SELECTION_Z_THRESHOLDS,
    analyze_lof_pc_enrichment,
)
from rare_variant_enrichment.phenotypes import prepare_phenotypes
from rare_variant_enrichment.statistics import calculate_enrichment
from rare_variant_enrichment.vat import prepare_vat
from rare_variant_enrichment.variants import classify_chromosome


COMMANDS = (
    "prepare-phenotypes",
    "prepare-vat",
    "classify-chromosome",
    "gather",
    "calculate",
    "prepare-protein-coding-genes",
    "lof-pc-enrichment",
    "carrier-pc-enrichment",
    "merge-lof-pc-enrichment",
    "analyze-lof-pc-enrichment",
    "pc-chunks",
    "prepare-carrier-inputs",
    "extract-gene-carriers",
    "gather-gene-carriers",
    "build-carrier-definitions",
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rare-variant-enrichment")
    subparsers = parser.add_subparsers(dest="command", required=True)
    coding_gene_parser = subparsers.add_parser("prepare-protein-coding-genes")
    coding_gene_parser.add_argument("--gtf", required=True, type=Path)
    coding_gene_parser.add_argument("--genes-output", required=True, type=Path)
    coding_gene_parser.add_argument("--qc-output", required=True, type=Path)
    lof_pc_parser = subparsers.add_parser("lof-pc-enrichment")
    lof_pc_parser.add_argument("--phenotype-bed", required=True, type=Path)
    lof_pc_parser.add_argument("--lof-carriers", required=True, type=Path)
    lof_pc_parser.add_argument("--principal-components", required=True, type=Path)
    lof_pc_parser.add_argument("--additional-covariates", type=Path)
    lof_pc_parser.add_argument("--protein-coding-genes", required=True, type=Path)
    lof_pc_parser.add_argument(
        "--negative-z-thresholds", required=True, type=parse_csv_floats
    )
    lof_pc_parser.add_argument("--pc-counts", required=True, type=parse_csv_ints)
    lof_pc_parser.add_argument("--results-output", required=True, type=Path)
    lof_pc_parser.add_argument("--summary-output", required=True, type=Path)
    lof_pc_parser.add_argument("--gene-pc-qc-output", required=True, type=Path)
    lof_pc_parser.add_argument("--analysis-qc-output", required=True, type=Path)
    lof_pc_parser.add_argument("--pc-grid-mode", choices=("adaptive", "explicit"))
    carrier_pc_parser = subparsers.add_parser("carrier-pc-enrichment")
    carrier_pc_parser.add_argument("--phenotype-bed", required=True, type=Path)
    carrier_pc_parser.add_argument("--carrier-table", required=True, type=Path)
    carrier_pc_parser.add_argument("--carrier-manifest", required=True, type=Path)
    carrier_pc_parser.add_argument("--principal-components", required=True, type=Path)
    carrier_pc_parser.add_argument("--additional-covariates", type=Path)
    carrier_pc_parser.add_argument("--protein-coding-genes", required=True, type=Path)
    carrier_pc_parser.add_argument(
        "--negative-z-thresholds", required=True, type=parse_csv_floats
    )
    carrier_pc_parser.add_argument("--pc-counts", required=True, type=parse_csv_ints)
    carrier_pc_parser.add_argument("--results-output", required=True, type=Path)
    carrier_pc_parser.add_argument("--summary-output", required=True, type=Path)
    carrier_pc_parser.add_argument("--gene-pc-qc-output", required=True, type=Path)
    carrier_pc_parser.add_argument("--analysis-qc-output", required=True, type=Path)
    carrier_pc_parser.add_argument("--pc-grid-mode", choices=("adaptive", "explicit"))
    merge_lof_pc_parser = subparsers.add_parser("merge-lof-pc-enrichment")
    merge_lof_pc_parser.add_argument("--results-input", action="append", required=True, type=Path)
    merge_lof_pc_parser.add_argument("--summary-input", action="append", required=True, type=Path)
    merge_lof_pc_parser.add_argument("--gene-pc-qc-input", action="append", required=True, type=Path)
    merge_lof_pc_parser.add_argument("--analysis-qc-input", action="append", required=True, type=Path)
    merge_lof_pc_parser.add_argument("--results-output", required=True, type=Path)
    merge_lof_pc_parser.add_argument("--summary-output", required=True, type=Path)
    merge_lof_pc_parser.add_argument("--gene-pc-qc-output", required=True, type=Path)
    merge_lof_pc_parser.add_argument("--analysis-qc-output", required=True, type=Path)
    analyze_lof_pc_parser = subparsers.add_parser("analyze-lof-pc-enrichment")
    analyze_lof_pc_parser.add_argument("--results-input", required=True, type=Path)
    analyze_lof_pc_parser.add_argument("--selection-output", required=True, type=Path)
    analyze_lof_pc_parser.add_argument("--plot-output", required=True, type=Path)
    analyze_lof_pc_parser.add_argument(
        "--selection-z-thresholds",
        type=parse_csv_floats,
        default=list(DEFAULT_SELECTION_Z_THRESHOLDS),
    )
    analyze_lof_pc_parser.add_argument(
        "--plateau-fraction", type=float, default=DEFAULT_PLATEAU_FRACTION
    )
    pc_chunks_parser = subparsers.add_parser("pc-chunks")
    pc_chunks_parser.add_argument("--principal-components", required=True, type=Path)
    pc_chunks_parser.add_argument("--pc-counts", required=True, type=parse_csv_ints)
    pc_chunks_parser.add_argument("--pc-counts-per-job", required=True, type=int)
    pc_chunks_parser.add_argument("--output", required=True, type=Path)
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
    classify_parser.add_argument(
        "--consequence-classes", required=True, type=parse_csv_consequences
    )
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
    calculate_parser.add_argument("--vat-schema", type=Path)
    calculate_parser.add_argument("--selected-chromosomes", type=parse_csv_strings)
    calculate_parser.add_argument("--container-image")
    calculate_parser.add_argument("--workflow-version", default="unknown")
    calculate_parser.add_argument("--max-retries", type=int, default=0)
    calculate_parser.add_argument(
        "--index-provenance",
        choices=("generated", "supplied", "unknown"),
        default="unknown",
    )
    calculate_parser.add_argument(
        "--consequence-classes", type=parse_csv_consequences, default=[]
    )
    calculate_parser.add_argument("--loftee-enabled", type=parse_boolean, default=False)
    calculate_parser.add_argument(
        "--vat-index-provenance",
        choices=("generated", "supplied", "unknown"),
        default="unknown",
    )
    calculate_parser.add_argument("--maximum-gvs-maf", type=float, default=0.01)
    calculate_parser.add_argument("--annotation-chunk-size-bp", type=int, default=10_000_000)

    prepare_carrier_parser = subparsers.add_parser("prepare-carrier-inputs")
    prepare_carrier_parser.add_argument("--vcf", required=True, type=Path)
    prepare_carrier_parser.add_argument("--annotations", required=True, type=Path)
    prepare_carrier_parser.add_argument(
        "--chromosomes", required=True, type=parse_csv_strings
    )
    prepare_carrier_parser.add_argument(
        "--vcf-index-provenance", required=True, choices=("supplied", "generated")
    )
    prepare_carrier_parser.add_argument(
        "--transcript-index-provenance",
        required=True,
        choices=("supplied", "generated"),
    )
    prepare_carrier_parser.add_argument("--schema-output", required=True, type=Path)
    prepare_carrier_parser.add_argument("--qc-output", required=True, type=Path)

    extract_carrier_parser = subparsers.add_parser("extract-gene-carriers")
    extract_carrier_parser.add_argument("--vcf", required=True, type=Path)
    extract_carrier_parser.add_argument("--annotations", required=True, type=Path)
    extract_carrier_parser.add_argument("--schema", required=True, type=Path)
    extract_carrier_parser.add_argument("--chromosome", required=True)
    extract_carrier_parser.add_argument(
        "--chunk-size-bp", required=True, type=parse_positive_int
    )
    extract_carrier_parser.add_argument("--audit-output", required=True, type=Path)
    extract_carrier_parser.add_argument("--qc-output", required=True, type=Path)

    gather_carrier_parser = subparsers.add_parser("gather-gene-carriers")
    gather_carrier_parser.add_argument(
        "--audit-input", action="append", required=True, type=Path
    )
    gather_carrier_parser.add_argument(
        "--qc-input", action="append", required=True, type=Path
    )
    gather_carrier_parser.add_argument("--preparation-qc", required=True, type=Path)
    gather_carrier_parser.add_argument("--audit-output", required=True, type=Path)
    gather_carrier_parser.add_argument("--carrier-output", required=True, type=Path)
    gather_carrier_parser.add_argument("--qc-output", required=True, type=Path)

    definition_parser = subparsers.add_parser("build-carrier-definitions")
    definition_parser.add_argument("--audit", required=True, type=Path)
    definition_parser.add_argument("--extraction-qc", required=True, type=Path)
    definition_parser.add_argument("--definitions", required=True, type=Path)
    definition_parser.add_argument("--container-image", required=True)
    definition_parser.add_argument("--output", required=True, type=Path)
    definition_parser.add_argument("--qc-output", required=True, type=Path)
    return parser


def main() -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    args = build_parser().parse_args()
    if args.command == "prepare-protein-coding-genes":
        prepare_protein_coding_genes(args.gtf, args.genes_output, args.qc_output)
    elif args.command == "lof-pc-enrichment":
        calculate_arguments = (
            args.phenotype_bed,
            args.lof_carriers,
            args.principal_components,
            args.protein_coding_genes,
            args.negative_z_thresholds,
            args.pc_counts,
            args.results_output,
            args.summary_output,
            args.gene_pc_qc_output,
            args.analysis_qc_output,
        )
        calculate_options = {}
        if args.pc_grid_mode is not None:
            calculate_options["pc_grid_mode"] = args.pc_grid_mode
        if args.additional_covariates is not None:
            calculate_options["additional_covariates_path"] = args.additional_covariates
        calculate_lof_pc_enrichment(*calculate_arguments, **calculate_options)
    elif args.command == "carrier-pc-enrichment":
        calculate_arguments = (
            args.phenotype_bed,
            args.carrier_table,
            args.carrier_manifest,
            args.principal_components,
            args.protein_coding_genes,
            args.negative_z_thresholds,
            args.pc_counts,
            args.results_output,
            args.summary_output,
            args.gene_pc_qc_output,
            args.analysis_qc_output,
        )
        calculate_options = {}
        if args.pc_grid_mode is not None:
            calculate_options["pc_grid_mode"] = args.pc_grid_mode
        if args.additional_covariates is not None:
            calculate_options["additional_covariates_path"] = (
                args.additional_covariates
            )
        calculate_carrier_pc_enrichment(*calculate_arguments, **calculate_options)
    elif args.command == "pc-chunks":
        available_pc_count = read_principal_component_header(args.principal_components)
        write_json(
            args.output,
            build_pc_chunks(
                args.pc_counts,
                available_pc_count,
                args.pc_counts_per_job,
            ),
        )
    elif args.command == "merge-lof-pc-enrichment":
        merge_lof_pc_enrichment(
            args.results_input,
            args.summary_input,
            args.gene_pc_qc_input,
            args.analysis_qc_input,
            args.results_output,
            args.summary_output,
            args.gene_pc_qc_output,
            args.analysis_qc_output,
        )
    elif args.command == "analyze-lof-pc-enrichment":
        analyze_lof_pc_enrichment(
            args.results_input,
            args.selection_output,
            args.plot_output,
            selection_z_thresholds=args.selection_z_thresholds,
            plateau_fraction=args.plateau_fraction,
        )
    elif args.command == "prepare-phenotypes":
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
            consequence_classes=args.consequence_classes,
            loftee_enabled=args.loftee_enabled,
            vat_index_provenance=args.vat_index_provenance,
            maximum_gvs_maf=args.maximum_gvs_maf,
            annotation_chunk_size_bp=args.annotation_chunk_size_bp,
            vat_schema_path=args.vat_schema,
        )
    elif args.command == "prepare-carrier-inputs":
        prepare_carrier_inputs(
            args.vcf,
            args.annotations,
            args.chromosomes,
            args.vcf_index_provenance,
            args.transcript_index_provenance,
            args.schema_output,
            args.qc_output,
        )
    elif args.command == "extract-gene-carriers":
        extract_chromosome_carriers(
            args.vcf,
            args.annotations,
            args.schema,
            args.chromosome,
            args.chunk_size_bp,
            args.audit_output,
            args.qc_output,
        )
    elif args.command == "gather-gene-carriers":
        gather_variant_carriers(
            args.audit_input,
            args.qc_input,
            args.preparation_qc,
            args.audit_output,
            args.carrier_output,
            args.qc_output,
        )
    elif args.command == "build-carrier-definitions":
        build_carrier_definitions(
            args.audit,
            args.extraction_qc,
            args.definitions,
            args.output,
            args.qc_output,
            container_image=args.container_image,
        )
    return 0


def parse_csv_strings(value: str) -> list[str]:
    values = [item.strip() for item in value.split(",")]
    if not values or any(not item for item in values):
        raise argparse.ArgumentTypeError("Expected a non-empty comma-separated list")
    return values


def parse_csv_consequences(value: str) -> list[str]:
    if value == "":
        return []
    return parse_csv_strings(value)


def parse_csv_ints(value: str) -> list[int]:
    if value == "":
        return []
    try:
        return [int(item) for item in parse_csv_strings(value)]
    except ValueError as error:
        raise argparse.ArgumentTypeError("Expected comma-separated integers") from error


def parse_positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("Expected a positive integer") from error
    if parsed < 1:
        raise argparse.ArgumentTypeError("Expected a positive integer")
    return parsed


def parse_boolean(value: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise argparse.ArgumentTypeError("Expected true or false")


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
