import tracemalloc
from pathlib import Path

from rare_variant_enrichment.aggregation import gather_outputs
from rare_variant_enrichment.statistics import calculate_enrichment, fisher_exact_two_sided


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
        handle.write("sample_id\tfeature_id\tac_class\tminimum_distance_bp\n")
        for index in range(40_000):
            handle.write(f"S{index % 1000}\tG{index}\tAC=1\t{index % 100}\n")
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
        carrier_handle.write("sample_id\tfeature_id\tac_class\tminimum_distance_bp\n")
        values = ["3", *(["0"] * 999)]
        for feature_index in range(30):
            feature_id = f"G{feature_index}"
            tss = 1000 + feature_index
            feature_handle.write(f"chr1\t{tss}\t{feature_id}\n")
            bed_handle.write(
                f"chr1\t{tss - 1}\t{tss}\t{feature_id}\t" + "\t".join(values) + "\n"
            )
            for sample_id in sample_ids:
                carrier_handle.write(f"{sample_id}\t{feature_id}\tAC=1\t0\n")

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
