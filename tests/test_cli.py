import subprocess
import sys
from pathlib import Path

import pytest

from rare_variant_enrichment.annotations import VatSchema
from rare_variant_enrichment import cli
from rare_variant_enrichment.cli import build_parser, parse_csv_ints


def test_cli_lists_workflow_subcommands():
    result = subprocess.run(
        [sys.executable, "-m", "rare_variant_enrichment.cli", "--help"],
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0
    for command in (
        "prepare-phenotypes",
        "prepare-vat",
        "classify-chromosome",
        "gather",
        "calculate",
        "prepare-protein-coding-genes",
        "lof-pc-enrichment",
        "merge-lof-pc-enrichment",
        "analyze-lof-pc-enrichment",
        "prepare-carrier-inputs",
        "extract-gene-carriers",
        "gather-gene-carriers",
        "build-carrier-definitions",
        "carrier-pc-enrichment",
        "merge-carrier-pc-enrichment",
        "analyze-carrier-pc-enrichment",
    ):
        assert command in result.stdout


def test_prepare_carrier_inputs_cli_dispatches_exact_contract(monkeypatch):
    received: list[object] = []
    monkeypatch.setattr(cli, "prepare_carrier_inputs", lambda *args: received.extend(args))
    monkeypatch.setattr(sys, "argv", [
        "rare-variant-enrichment", "prepare-carrier-inputs",
        "--vcf", "variants.vcf.gz",
        "--annotations", "transcript.tsv.bgz",
        "--chromosomes", "chr1,chr2",
        "--vcf-index-provenance", "supplied",
        "--transcript-index-provenance", "generated",
        "--schema-output", "transcript.schema.json",
        "--qc-output", "transcript.prepare.qc.json",
    ])

    assert cli.main() == 0
    assert received == [
        Path("variants.vcf.gz"), Path("transcript.tsv.bgz"), ["chr1", "chr2"],
        "supplied", "generated", Path("transcript.schema.json"),
        Path("transcript.prepare.qc.json"),
    ]


def test_extract_gene_carriers_cli_dispatches_exact_contract(monkeypatch):
    received: list[object] = []
    monkeypatch.setattr(cli, "extract_chromosome_carriers", lambda *args: received.extend(args))
    monkeypatch.setattr(sys, "argv", [
        "rare-variant-enrichment", "extract-gene-carriers",
        "--vcf", "variants.vcf.gz", "--annotations", "transcript.tsv.bgz",
        "--schema", "transcript.schema.json", "--chromosome", "chr1",
        "--chunk-size-bp", "10000000", "--audit-output", "chr1.audit.tsv.gz",
        "--qc-output", "chr1.qc.json",
    ])

    assert cli.main() == 0
    assert received == [
        Path("variants.vcf.gz"), Path("transcript.tsv.bgz"),
        Path("transcript.schema.json"), "chr1", 10_000_000,
        Path("chr1.audit.tsv.gz"), Path("chr1.qc.json"),
    ]


def test_gather_gene_carriers_cli_dispatches_exact_contract(monkeypatch):
    received: list[object] = []
    monkeypatch.setattr(cli, "gather_variant_carriers", lambda *args: received.extend(args))
    monkeypatch.setattr(sys, "argv", [
        "rare-variant-enrichment", "gather-gene-carriers",
        "--audit-input", "chr1.audit.tsv.gz", "--qc-input", "chr1.qc.json",
        "--preparation-qc", "transcript.prepare.qc.json",
        "--audit-output", "variant_carrier_audit.tsv.gz",
        "--carrier-output", "variant_carriers.tsv.gz",
        "--qc-output", "variant_carriers.qc.json",
    ])

    assert cli.main() == 0
    assert received == [
        [Path("chr1.audit.tsv.gz")], [Path("chr1.qc.json")],
        Path("transcript.prepare.qc.json"), Path("variant_carrier_audit.tsv.gz"),
        Path("variant_carriers.tsv.gz"), Path("variant_carriers.qc.json"),
    ]


def test_build_carrier_definitions_cli_dispatches_exact_contract(monkeypatch):
    received: list[object] = []
    keywords: dict[str, object] = {}

    def fake_build(*arguments: object, **keyword_arguments: object) -> None:
        received.extend(arguments)
        keywords.update(keyword_arguments)

    monkeypatch.setattr(cli, "build_carrier_definitions", fake_build)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "build-carrier-definitions",
            "--audit",
            "variant_carrier_audit.tsv.gz",
            "--extraction-qc",
            "variant_carriers.qc.json",
            "--definitions",
            "definitions.json",
            "--container-image",
            "image@sha256:abc",
            "--output",
            "carrier_definitions.tsv.gz",
            "--qc-output",
            "carrier_definitions.qc.json",
        ],
    )

    assert cli.main() == 0
    assert received == [
        Path("variant_carrier_audit.tsv.gz"),
        Path("variant_carriers.qc.json"),
        Path("definitions.json"),
        Path("carrier_definitions.tsv.gz"),
        Path("carrier_definitions.qc.json"),
    ]
    assert keywords == {"container_image": "image@sha256:abc"}


def test_prepare_protein_coding_genes_cli_dispatches_exact_paths(monkeypatch):
    received: list[Path] = []

    def fake_prepare(gtf: Path, genes_output: Path, qc_output: Path) -> None:
        received.extend([gtf, genes_output, qc_output])

    monkeypatch.setattr(cli, "prepare_protein_coding_genes", fake_prepare)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "prepare-protein-coding-genes",
            "--gtf",
            "annotation.gtf.gz",
            "--genes-output",
            "genes.tsv",
            "--qc-output",
            "genes.qc.json",
        ],
    )

    assert cli.main() == 0
    assert received == [
        Path("annotation.gtf.gz"),
        Path("genes.tsv"),
        Path("genes.qc.json"),
    ]


def test_lof_pc_enrichment_cli_dispatches_exact_contract(monkeypatch):
    received: list[object] = []

    def fake_calculate(*arguments: object) -> None:
        received.extend(arguments)

    monkeypatch.setattr(cli, "calculate_lof_pc_enrichment", fake_calculate)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "lof-pc-enrichment",
            "--phenotype-bed",
            "phenotypes.bed.gz",
            "--lof-carriers",
            "lof.tsv",
            "--principal-components",
            "pcs.tsv",
            "--protein-coding-genes",
            "genes.tsv",
            "--negative-z-thresholds=-2,-3",
            "--pc-counts",
            "0,10",
            "--results-output",
            "results.tsv",
            "--summary-output",
            "summary.json",
            "--gene-pc-qc-output",
            "gene-pc-qc.tsv.gz",
            "--analysis-qc-output",
            "analysis-qc.json",
        ],
    )

    assert cli.main() == 0
    assert received == [
        Path("phenotypes.bed.gz"),
        Path("lof.tsv"),
        Path("pcs.tsv"),
        Path("genes.tsv"),
        [-2.0, -3.0],
        [0, 10],
        Path("results.tsv"),
        Path("summary.json"),
        Path("gene-pc-qc.tsv.gz"),
        Path("analysis-qc.json"),
    ]


def test_carrier_pc_enrichment_cli_dispatches_exact_contract(monkeypatch):
    received: list[object] = []
    keywords: dict[str, object] = {}

    def fake_calculate(*arguments: object, **keyword_arguments: object) -> None:
        received.extend(arguments)
        keywords.update(keyword_arguments)

    monkeypatch.setattr(cli, "calculate_carrier_pc_enrichment", fake_calculate)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "carrier-pc-enrichment",
            "--phenotype-bed",
            "phenotypes.bed.gz",
            "--carrier-table",
            "carrier_definitions.tsv.gz",
            "--carrier-manifest",
            "carrier_definitions.qc.json",
            "--principal-components",
            "pcs.tsv",
            "--additional-covariates",
            "covariates.tsv",
            "--protein-coding-genes",
            "genes.tsv",
            "--negative-z-thresholds=-2,-3",
            "--pc-counts",
            "0,10",
            "--pc-grid-mode",
            "adaptive",
            "--results-output",
            "results.tsv",
            "--summary-output",
            "summary.json",
            "--gene-pc-qc-output",
            "gene-pc-qc.tsv.gz",
            "--analysis-qc-output",
            "analysis-qc.json",
        ],
    )

    assert cli.main() == 0
    assert received == [
        Path("phenotypes.bed.gz"),
        Path("carrier_definitions.tsv.gz"),
        Path("carrier_definitions.qc.json"),
        Path("pcs.tsv"),
        Path("genes.tsv"),
        [-2.0, -3.0],
        [0, 10],
        Path("results.tsv"),
        Path("summary.json"),
        Path("gene-pc-qc.tsv.gz"),
        Path("analysis-qc.json"),
    ]
    assert keywords == {
        "pc_grid_mode": "adaptive",
        "additional_covariates_path": Path("covariates.tsv"),
    }


def test_lof_pc_enrichment_cli_dispatches_optional_additional_covariates(
    monkeypatch,
):
    received: list[object] = []
    keywords: dict[str, object] = {}

    def fake_calculate(*arguments: object, **keyword_arguments: object) -> None:
        received.extend(arguments)
        keywords.update(keyword_arguments)

    monkeypatch.setattr(cli, "calculate_lof_pc_enrichment", fake_calculate)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "lof-pc-enrichment",
            "--phenotype-bed",
            "phenotypes.bed.gz",
            "--lof-carriers",
            "lof.tsv",
            "--principal-components",
            "pcs.tsv",
            "--additional-covariates",
            "genetic-pcs.tsv",
            "--protein-coding-genes",
            "genes.tsv",
            "--negative-z-thresholds=-2",
            "--pc-counts",
            "0",
            "--results-output",
            "results.tsv",
            "--summary-output",
            "summary.json",
            "--gene-pc-qc-output",
            "gene-pc-qc.tsv.gz",
            "--analysis-qc-output",
            "analysis-qc.json",
        ],
    )

    assert cli.main() == 0
    assert received[:4] == [
        Path("phenotypes.bed.gz"),
        Path("lof.tsv"),
        Path("pcs.tsv"),
        Path("genes.tsv"),
    ]
    assert keywords == {
        "additional_covariates_path": Path("genetic-pcs.tsv"),
    }


def test_pc_chunks_cli_dispatches_header_count_chunking_and_json_output(monkeypatch):
    received: list[object] = []

    def fake_read_header(principal_components: Path) -> int:
        received.append(principal_components)
        return 30

    def fake_build_chunks(
        pc_counts: list[int], available_pc_count: int, pc_counts_per_job: int
    ) -> list[list[int]]:
        received.extend([pc_counts, available_pc_count, pc_counts_per_job])
        return [[0, 1], [10]]

    def fake_write_json(output: Path, chunks: list[list[int]]) -> None:
        received.extend([output, chunks])

    monkeypatch.setattr(cli, "read_principal_component_header", fake_read_header)
    monkeypatch.setattr(cli, "build_pc_chunks", fake_build_chunks)
    monkeypatch.setattr(cli, "write_json", fake_write_json)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "pc-chunks",
            "--principal-components",
            "pcs.tsv",
            "--pc-counts",
            "0,1,10",
            "--pc-counts-per-job",
            "2",
            "--output",
            "chunks.json",
        ],
    )

    assert cli.main() == 0
    assert received == [Path("pcs.tsv"), [0, 1, 10], 30, 2, Path("chunks.json"), [[0, 1], [10]]]


def test_merge_lof_pc_enrichment_cli_dispatches_direct_inputs(monkeypatch):
    received: list[object] = []

    def fake_merge(*arguments: object) -> None:
        received.extend(arguments)

    monkeypatch.setattr(cli, "merge_lof_pc_enrichment", fake_merge)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "merge-lof-pc-enrichment",
            "--results-input",
            "/cromwell/results-0.tsv",
            "--results-input",
            "/cromwell/results-1.tsv",
            "--summary-input",
            "/cromwell/summary-0.json",
            "--summary-input",
            "/cromwell/summary-1.json",
            "--gene-pc-qc-input",
            "/cromwell/gene-0.tsv.gz",
            "--gene-pc-qc-input",
            "/cromwell/gene-1.tsv.gz",
            "--analysis-qc-input",
            "/cromwell/analysis-0.json",
            "--analysis-qc-input",
            "/cromwell/analysis-1.json",
            "--results-output",
            "results.tsv",
            "--summary-output",
            "summary.json",
            "--gene-pc-qc-output",
            "gene-pc-qc.tsv.gz",
            "--analysis-qc-output",
            "analysis-qc.json",
        ],
    )

    assert cli.main() == 0
    assert received == [
        [Path("/cromwell/results-0.tsv"), Path("/cromwell/results-1.tsv")],
        [Path("/cromwell/summary-0.json"), Path("/cromwell/summary-1.json")],
        [Path("/cromwell/gene-0.tsv.gz"), Path("/cromwell/gene-1.tsv.gz")],
        [Path("/cromwell/analysis-0.json"), Path("/cromwell/analysis-1.json")],
        Path("results.tsv"),
        Path("summary.json"),
        Path("gene-pc-qc.tsv.gz"),
        Path("analysis-qc.json"),
    ]


def test_merge_carrier_pc_enrichment_cli_dispatches_direct_inputs(monkeypatch):
    received: list[object] = []
    monkeypatch.setattr(
        cli, "merge_carrier_pc_enrichment", lambda *args: received.extend(args)
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "merge-carrier-pc-enrichment",
            "--results-input",
            "results-0.tsv",
            "--summary-input",
            "summary-0.json",
            "--gene-pc-qc-input",
            "gene-0.tsv.gz",
            "--analysis-qc-input",
            "analysis-0.json",
            "--results-output",
            "results.tsv",
            "--summary-output",
            "summary.json",
            "--gene-pc-qc-output",
            "gene.tsv.gz",
            "--analysis-qc-output",
            "analysis.json",
        ],
    )

    assert cli.main() == 0
    assert received == [
        [Path("results-0.tsv")],
        [Path("summary-0.json")],
        [Path("gene-0.tsv.gz")],
        [Path("analysis-0.json")],
        Path("results.tsv"),
        Path("summary.json"),
        Path("gene.tsv.gz"),
        Path("analysis.json"),
    ]


def test_analyze_lof_pc_enrichment_cli_dispatches_selection_options(monkeypatch):
    received: list[object] = []

    def fake_analyze(*arguments: object, **keyword_arguments: object) -> None:
        received.extend([arguments, keyword_arguments])

    monkeypatch.setattr(cli, "analyze_lof_pc_enrichment", fake_analyze)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "analyze-lof-pc-enrichment",
            "--results-input",
            "results.tsv",
            "--selection-output",
            "selection.json",
            "--plot-output",
            "plot.svg",
            "--selection-z-thresholds=-3,-4,-5,-6",
            "--plateau-fraction",
            "0.95",
        ],
    )

    assert cli.main() == 0
    assert received == [
        (Path("results.tsv"), Path("selection.json"), Path("plot.svg")),
        {"selection_z_thresholds": [-3.0, -4.0, -5.0, -6.0], "plateau_fraction": 0.95},
    ]


def test_analyze_carrier_pc_enrichment_cli_dispatches_dynamic_definitions(monkeypatch):
    received: list[object] = []

    def fake_analyze(*arguments: object, **keyword_arguments: object) -> None:
        received.extend([arguments, keyword_arguments])

    monkeypatch.setattr(cli, "analyze_carrier_pc_enrichment", fake_analyze)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "analyze-carrier-pc-enrichment",
            "--results-input",
            "results.tsv",
            "--selection-output",
            "selection.json",
            "--plot-output",
            "plot.svg",
            "--carrier-definitions",
            "lof_hc,missense,splice_any",
            "--selection-z-thresholds=-3,-4",
            "--plateau-fraction",
            "0.95",
        ],
    )

    assert cli.main() == 0
    assert received == [
        (Path("results.tsv"), Path("selection.json"), Path("plot.svg")),
        {
            "carrier_definitions": ["lof_hc", "missense", "splice_any"],
            "selection_z_thresholds": [-3.0, -4.0],
            "plateau_fraction": 0.95,
        },
    ]


def test_prepare_vat_cli_dispatches_paths_and_chromosomes_unchanged(tmp_path: Path, monkeypatch):
    vat = tmp_path / "annotations.tsv.bgz"
    schema_output = tmp_path / "schema.json"
    loftee_enabled_output = tmp_path / "loftee_enabled.txt"
    chromosomes = ["chr1", "chrX"]
    schema = VatSchema.from_header(
        ["chrom", "pos", "ref", "alt", "gene_id", "consequence", "gvs_max_af", "LoF"]
    )
    received: list[object] = []

    def fake_prepare_vat(
        received_vat: Path, received_chromosomes: list[str], received_schema_output: Path
    ) -> VatSchema:
        received.extend([received_vat, received_chromosomes, received_schema_output])
        return schema

    monkeypatch.setattr(cli, "prepare_vat", fake_prepare_vat)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "prepare-vat",
            "--vat",
            str(vat),
            "--chromosomes",
            "chr1,chrX",
            "--schema-output",
            str(schema_output),
            "--loftee-enabled-output",
            str(loftee_enabled_output),
        ],
    )

    assert cli.main() == 0
    assert received == [vat, chromosomes, schema_output]
    assert loftee_enabled_output.read_text() == "true\n"


def test_gather_cli_writes_aggregated_outputs(tmp_path: Path):
    carrier = tmp_path / "chr1.tsv"
    carrier.write_text(
        "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
        "S1\tGENE1\tAC=1\tbaseline\tall_rare_variants\t4\n"
    )
    qc = tmp_path / "chr1.json"
    qc.write_text('{"chromosome":"chr1","extracted_records":1}')
    carrier_output = tmp_path / "all.tsv"
    qc_output = tmp_path / "qc.tsv"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "rare_variant_enrichment.cli",
            "gather",
            "--carrier-input",
            str(carrier),
            "--qc-input",
            str(qc),
            "--carrier-output",
            str(carrier_output),
            "--qc-output",
            str(qc_output),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "S1\tGENE1\tAC=1\tbaseline\tall_rare_variants\t4" in carrier_output.read_text().splitlines()
    assert qc_output.read_text().splitlines() == ["chromosome\textracted_records", "chr1\t1"]


def test_classify_cli_dispatches_vat_annotation_arguments(tmp_path: Path, monkeypatch):
    received: list[object] = []

    def fake_classify(*arguments):
        received.extend(arguments)

    monkeypatch.setattr(cli, "classify_chromosome", fake_classify)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "classify-chromosome",
            "--vcf", "variants.vcf.gz",
            "--vat", "annotations.tsv.bgz",
            "--vat-schema", "vat-schema.json",
            "--features", "features.tsv",
            "--shared-samples", "shared.txt",
            "--chromosome", "chr1",
            "--exact-ac", "1,2",
            "--cumulative-ac-max", "1,2",
            "--consequence-classes", "stop_gained,missense_variant",
            "--maximum-gvs-maf", "0.01",
            "--max-distance", "100",
            "--annotation-chunk-size-bp", "25",
            "--carrier-output", "carriers.tsv",
            "--regions-output", "regions.bed",
            "--qc-output", "qc.json",
        ],
    )

    assert cli.main() == 0
    assert received == [
        Path("variants.vcf.gz"),
        Path("annotations.tsv.bgz"),
        Path("vat-schema.json"),
        Path("features.tsv"),
        Path("shared.txt"),
        "chr1",
        [1, 2],
        [1, 2],
        ["stop_gained", "missense_variant"],
        0.01,
        100,
        25,
        Path("carriers.tsv"),
        Path("regions.bed"),
        Path("qc.json"),
    ]


def test_calculate_cli_dispatches_annotation_configuration_as_literal_data(
    tmp_path: Path, monkeypatch
):
    """Replacing direct dispatch with shell evaluation must fail this contract."""
    received: dict[str, object] = {}
    marker = tmp_path / "must-not-exist"
    consequence = f"stop_gained;touch {marker}"

    def fake_calculate(*arguments: object, **keywords: object) -> None:
        received["arguments"] = arguments
        received["keywords"] = keywords

    monkeypatch.setattr(cli, "calculate_enrichment", fake_calculate)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "rare-variant-enrichment",
            "calculate",
            "--phenotype-bed", "phenotypes.bed",
            "--shared-samples", "shared_samples.txt",
            "--carriers", "carriers.tsv",
            "--features", "features.tsv",
            "--exact-ac", "1",
            "--cumulative-ac-max", "",
            "--z-thresholds", "2.0",
            "--distance-thresholds", "100",
            "--tail", "absolute",
            "--output-tsv", "enrichment.tsv",
            "--output-json", "enrichment.json",
            "--vat-schema", "vat-schema.json",
            "--consequence-classes", consequence,
            "--loftee-enabled", "true",
            "--vat-index-provenance", "supplied",
            "--maximum-gvs-maf", "0.01",
            "--annotation-chunk-size-bp", "500",
        ],
    )

    assert cli.main() == 0
    assert received["keywords"] == {
        "phenotype_qc_path": None,
        "chromosome_qc_path": None,
        "selected_chromosomes": None,
        "container_image": None,
        "workflow_version": "unknown",
        "max_retries": 0,
        "index_provenance": "unknown",
        "consequence_classes": [consequence],
        "loftee_enabled": True,
        "vat_index_provenance": "supplied",
        "maximum_gvs_maf": 0.01,
        "annotation_chunk_size_bp": 500,
        "vat_schema_path": Path("vat-schema.json"),
    }
    assert not marker.exists()


def test_parse_csv_ints_decodes_empty_wdl_array():
    assert parse_csv_ints("") == []


@pytest.mark.parametrize("command", ["classify-chromosome", "calculate"])
def test_consequence_cli_options_decode_empty_wdl_array(command: str):
    if command == "classify-chromosome":
        arguments = [
            command,
            "--vcf", "variants.vcf.gz",
            "--vat", "annotations.tsv.bgz",
            "--vat-schema", "vat-schema.json",
            "--features", "features.tsv",
            "--shared-samples", "shared_samples.txt",
            "--chromosome", "chr1",
            "--exact-ac", "1",
            "--cumulative-ac-max", "",
            "--consequence-classes", "",
            "--maximum-gvs-maf", "0.01",
            "--max-distance", "1000",
            "--annotation-chunk-size-bp", "10000000",
            "--carrier-output", "carriers.tsv",
            "--regions-output", "regions.bed",
            "--qc-output", "qc.json",
        ]
    else:
        arguments = [
            command,
            "--phenotype-bed", "phenotypes.bed",
            "--shared-samples", "shared_samples.txt",
            "--carriers", "carriers.tsv",
            "--features", "features.tsv",
            "--exact-ac", "1",
            "--cumulative-ac-max", "",
            "--z-thresholds", "2.0",
            "--distance-thresholds", "1000",
            "--tail", "absolute",
            "--output-tsv", "enrichment.tsv",
            "--output-json", "enrichment.json",
            "--consequence-classes", "",
        ]

    parsed = build_parser().parse_args(arguments)

    assert parsed.consequence_classes == []


@pytest.mark.parametrize("value", [",stop_gained", "stop_gained,", "stop_gained,,missense_variant"])
def test_consequence_cli_options_keep_nonempty_csv_validation(value: str):
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "calculate",
                "--phenotype-bed", "phenotypes.bed",
                "--shared-samples", "shared_samples.txt",
                "--carriers", "carriers.tsv",
                "--features", "features.tsv",
                "--exact-ac", "1",
                "--cumulative-ac-max", "",
                "--z-thresholds", "2.0",
                "--distance-thresholds", "1000",
                "--tail", "absolute",
                "--output-tsv", "enrichment.tsv",
                "--output-json", "enrichment.json",
                "--consequence-classes", value,
            ]
        )


@pytest.mark.parametrize(
    ("command", "exact_ac", "cumulative_ac_max"),
    [
        ("classify-chromosome", "", "1"),
        ("classify-chromosome", "1", ""),
        ("calculate", "", "1"),
        ("calculate", "1", ""),
    ],
)
def test_ac_cli_options_allow_one_empty_family(
    command: str, exact_ac: str, cumulative_ac_max: str
):
    common = [
        "--exact-ac",
        exact_ac,
        "--cumulative-ac-max",
        cumulative_ac_max,
    ]
    if command == "classify-chromosome":
        arguments = [
            command,
            "--vcf",
            "variants.vcf.gz",
            "--vat",
            "annotations.tsv.bgz",
            "--vat-schema",
            "vat-schema.json",
            "--features",
            "features.tsv",
            "--shared-samples",
            "shared_samples.txt",
            "--chromosome",
            "chr1",
            *common,
            "--consequence-classes",
            "stop_gained",
            "--maximum-gvs-maf",
            "0.01",
            "--max-distance",
            "1000",
            "--annotation-chunk-size-bp",
            "10000000",
            "--carrier-output",
            "carriers.tsv",
            "--regions-output",
            "regions.bed",
            "--qc-output",
            "qc.json",
        ]
    else:
        arguments = [
            command,
            "--phenotype-bed",
            "phenotypes.bed",
            "--shared-samples",
            "shared_samples.txt",
            "--carriers",
            "carriers.tsv",
            "--features",
            "features.tsv",
            *common,
            "--z-thresholds",
            "2.0",
            "--distance-thresholds",
            "1000",
            "--tail",
            "absolute",
            "--output-tsv",
            "enrichment.tsv",
            "--output-json",
            "enrichment.json",
        ]

    parsed = build_parser().parse_args(arguments)

    assert parsed.exact_ac == ([] if exact_ac == "" else [1])
    assert parsed.cumulative_ac_max == ([] if cumulative_ac_max == "" else [1])
