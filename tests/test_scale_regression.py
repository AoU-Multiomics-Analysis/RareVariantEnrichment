import tracemalloc
from pathlib import Path

from rare_variant_enrichment.aggregation import gather_outputs
from rare_variant_enrichment.annotation_storage import VatChunkStore
from rare_variant_enrichment.annotations import VatSchema, VariantKey
from rare_variant_enrichment.statistics import calculate_enrichment, fisher_exact_two_sided


MTTOVCF_HEADER = [
    "chrom", "pos", "ref", "alt", "rsid", "gene_id", "gene_symbol",
    "transcript", "is_canonical_transcript", "consequence", "aa_change",
    "LoF", "LoF_filter", "LoF_flags", "LoF_info", "gvs_max_af",
    "gvs_max_subpop",
]


def test_fisher_wide_support_uses_constant_python_memory():
    tracemalloc.start()
    try:
        assert fisher_exact_two_sided(125_000, 125_000, 125_000, 125_000) == 1.0
        _, peak_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    print(f"fisher_wide_support_peak_mib={peak_bytes / 1024 / 1024:.3f}")
    assert peak_bytes < 1 * 1024 * 1024


def test_gather_reduction_has_bounded_python_heap(tmp_path: Path):
    carrier = tmp_path / "chr1.tsv"
    with carrier.open("w", encoding="utf-8") as handle:
        handle.write(
            "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
        )
        for index in range(40_000):
            handle.write(
                f"S{index % 1000}\tG{index}\tAC=1\tbaseline\tall_rare_variants\t{index % 100}\n"
            )
    qc = tmp_path / "chr1.json"
    qc.write_text('{"chromosome":"chr1","emitted_keys":40000}')

    tracemalloc.start()
    try:
        gather_outputs(
            [carrier],
            [qc],
            tmp_path / "gathered.tsv",
            tmp_path / "chromosome_qc.tsv",
        )
        _, peak_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    print(f"gather_40000_keys_peak_mib={peak_bytes / 1024 / 1024:.3f}")
    assert (tmp_path / "gathered.tsv").read_text().count("\n") == 40_001
    assert peak_bytes < 4 * 1024 * 1024


def test_calculation_streams_thirty_thousand_carrier_keys_with_bounded_heap(
    tmp_path: Path,
):
    sample_ids = [f"S{index}" for index in range(1000)]
    features = tmp_path / "features.tsv"
    bed = tmp_path / "phenotypes.bed"
    carriers = tmp_path / "carriers.tsv"
    shared = tmp_path / "shared.txt"
    shared.write_text("\n".join(sample_ids) + "\n")
    with features.open("w", encoding="utf-8") as feature_handle, bed.open(
        "w", encoding="utf-8"
    ) as bed_handle, carriers.open("w", encoding="utf-8") as carrier_handle:
        feature_handle.write("chrom\ttss\tfeature_id\n")
        bed_handle.write("#chr\tstart\tend\tgene_id\t" + "\t".join(sample_ids) + "\n")
        carrier_handle.write(
            "sample_id\tfeature_id\tac_class\tannotation_family\tannotation_class\tminimum_distance_bp\n"
        )
        values = ["3", *(["0"] * 999)]
        for feature_index in range(30):
            feature_id = f"G{feature_index}"
            tss = 1000 + feature_index
            feature_handle.write(f"chr1\t{tss}\t{feature_id}\n")
            bed_handle.write(
                f"chr1\t{tss - 1}\t{tss}\t{feature_id}\t" + "\t".join(values) + "\n"
            )
            for sample_id in sample_ids:
                carrier_handle.write(
                    f"{sample_id}\t{feature_id}\tAC=1\tbaseline\tall_rare_variants\t0\n"
                )

    tracemalloc.start()
    try:
        calculate_enrichment(
            bed,
            shared,
            carriers,
            features,
            [1],
            [],
            [2.0],
            [0],
            "absolute",
            tmp_path / "enrichment.tsv",
            tmp_path / "summary.json",
        )
        _, peak_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    print(f"calculate_30000_keys_peak_mib={peak_bytes / 1024 / 1024:.3f}")
    assert peak_bytes < 8 * 1024 * 1024


def test_vat_chunk_store_streams_two_hundred_thousand_transcripts_with_bounded_heap(
    tmp_path: Path,
):
    """A chromosome-sized VAT chunk must remain disk-backed during transcript collapse."""
    transcript_rows = tmp_path / "vat-chunk.tsv"
    with transcript_rows.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(MTTOVCF_HEADER) + "\n")
        for allele_gene_index in range(25_000):
            position = allele_gene_index + 1
            for transcript_index in range(8):
                handle.write(
                    "\t".join(
                        [
                            "chr1",
                            str(position),
                            "A",
                            "C",
                            f"rs{allele_gene_index}",
                            f"ENSG{allele_gene_index:011d}.1",
                            f"GENE{allele_gene_index}",
                            f"ENST{allele_gene_index:011d}{transcript_index}",
                            "true",
                            "missense_variant",
                            "p.Ala1Val",
                            "HC",
                            ".",
                            ".",
                            ".",
                            "0.001",
                            "afr",
                        ]
                    )
                    + "\n"
                )

    schema = VatSchema.from_header(MTTOVCF_HEADER)
    tracemalloc.start()
    try:
        with VatChunkStore(
            tmp_path, schema, 0.01, ["missense_variant"]
        ) as chunk_store:
            with transcript_rows.open(encoding="utf-8") as handle:
                next(handle)
                for raw_row in handle:
                    chunk_store.ingest(raw_row.rstrip("\r\n").split("\t"))
            qc = chunk_store.finalize()
            assert qc["vat_rows"] == 200_000
            assert qc["unique_vat_allele_gene_pairs"] == 25_000
            for allele_gene_index in (0, 12_499, 24_999):
                key = VariantKey("chr1", allele_gene_index + 1, "A", "C")
                assert chunk_store.qualifying_maf(key) == 0.001
                assert chunk_store.gene_annotation(
                    key, f"ENSG{allele_gene_index:011d}"
                ).consequence == "missense_variant"
            _, peak_bytes = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    print(f"vat_chunk_200000_rows_peak_mib={peak_bytes / 1024 / 1024:.3f}")
    assert peak_bytes < 12 * 1024 * 1024
