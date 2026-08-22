"""Prepare and extract gene-matched carrier records from indexed inputs."""

from __future__ import annotations

import csv
import gzip
from io import TextIOWrapper
from pathlib import Path
import re
import subprocess
from typing import Callable, Iterator, Sequence

from rare_variant_enrichment.annotations import VariantKey
from rare_variant_enrichment.carrier_annotation_storage import CarrierAnnotationChunkStore
from rare_variant_enrichment.carrier_annotations import (
    TranscriptCarrierSchema,
    initial_variant_classes,
)
from rare_variant_enrichment.io import open_text, write_json
from rare_variant_enrichment.variants import QueryChunk, parse_variant_alleles


AUDIT_HEADER = (
    "sample_id",
    "gene_id",
    "gene_symbol",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_id",
    "variant_ac",
    "variant_af",
    "sample_alt_allele_count",
    "most_severe_consequence",
    "all_consequences",
    "unknown_consequences",
    "loftee",
    "revel",
    "gvs_max_af",
    "variant_classes",
)


def build_chromosome_chunks(
    chromosome: str, chromosome_length: int, chunk_size_bp: int
) -> tuple[QueryChunk, ...]:
    if (
        isinstance(chunk_size_bp, bool)
        or not isinstance(chunk_size_bp, int)
        or chunk_size_bp < 1
    ):
        raise ValueError("chunk_size_bp must be a positive integer")
    if chromosome_length < 1:
        raise ValueError("VCF contig length must be positive")
    return tuple(
        QueryChunk(chromosome, start, min(chromosome_length, start + chunk_size_bp - 1))
        for start in range(1, chromosome_length + 1, chunk_size_bp)
    )


def prepare_carrier_inputs(
    vcf_path: Path,
    annotation_path: Path,
    chromosomes: Sequence[str],
    vcf_index_provenance: str,
    annotation_index_provenance: str,
    schema_output: Path,
    qc_output: Path,
) -> None:
    if len(set(chromosomes)) != len(chromosomes):
        raise ValueError("Requested chromosomes must be unique")
    if not chromosomes:
        raise ValueError("At least one chromosome is required")
    for label, value in (
        ("VCF", vcf_index_provenance),
        ("transcript", annotation_index_provenance),
    ):
        if value not in {"supplied", "generated"}:
            raise ValueError(f"{label} index provenance must be supplied or generated")

    sample_ids, contig_lengths, vcf_header = _read_vcf_header(vcf_path)
    missing_contigs = [chromosome for chromosome in chromosomes if chromosome not in contig_lengths]
    if missing_contigs:
        raise ValueError("Requested chromosomes are absent from the VCF: " + ", ".join(missing_contigs))
    with open_text(annotation_path) as handle:
        annotation_header = _first_nonempty_fields(handle, "Transcript annotations are empty")
    schema = TranscriptCarrierSchema.from_header(annotation_header)
    schema.write_json(schema_output)
    write_json(
        qc_output,
        {
            "filtered_vcf": str(vcf_path),
            "transcript_annotations": str(annotation_path),
            "vcf_header": vcf_header,
            "transcript_header": list(schema.header),
            "transcript_schema": schema.as_dict(),
            "selected_chromosomes": list(chromosomes),
            "vcf_contig_lengths": {
                chromosome: contig_lengths[chromosome] for chromosome in chromosomes
            },
            "sample_count": len(sample_ids),
            "revel_available": True,
            "vcf_index_provenance": vcf_index_provenance,
            "transcript_index_provenance": annotation_index_provenance,
            "quality_or_frequency_filters_applied": False,
        },
    )


def extract_chromosome_carriers(
    vcf_path: Path,
    annotation_path: Path,
    schema_path: Path,
    chromosome: str,
    chunk_size_bp: int,
    audit_output: Path,
    qc_output: Path,
) -> None:
    sample_ids, contig_lengths, _ = _read_vcf_header(vcf_path)
    if chromosome not in contig_lengths:
        raise ValueError(f"Requested chromosome is absent from the VCF: {chromosome}")
    schema = TranscriptCarrierSchema.read_json(schema_path)
    chunks = build_chromosome_chunks(chromosome, contig_lengths[chromosome], chunk_size_bp)
    qc: dict[str, object] = {
        "chromosome": chromosome,
        "chunk_size_bp": chunk_size_bp,
        "annotation_chunk_count": len(chunks),
        "vcf_records": 0,
        "alt_alleles": 0,
        "vat_joined_alt_alleles": 0,
        "vat_unmatched_alt_alleles": 0,
        "carrier_audit_rows": 0,
        "transcript_rows": 0,
        "duplicate_transcript_rows": 0,
        "missing_gene_transcript_rows": 0,
        "unique_annotation_alleles": 0,
        "unique_annotation_allele_gene_pairs": 0,
        "quality_or_frequency_filters_applied": False,
    }

    with _open_deterministic_gzip_text(audit_output) as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for chunk in chunks:
            with CarrierAnnotationChunkStore(audit_output.parent, schema) as annotations:
                _stream_tabix_tsv(annotation_path, chunk.tabix_region, annotations.ingest)
                chunk_qc = annotations.finalize()
                for key, value in chunk_qc.items():
                    qc[key] = int(qc[key]) + value

                for fields in _stream_tabix_vcf(vcf_path, chunk.tabix_region):
                    qc["vcf_records"] = int(qc["vcf_records"]) + 1
                    alleles = parse_variant_alleles(
                        fields, sample_ids, set(sample_ids)
                    )
                    qc["alt_alleles"] = int(qc["alt_alleles"]) + len(alleles)
                    for allele in alleles:
                        if allele.chrom != chromosome:
                            raise ValueError(
                                f"Tabix returned chromosome {allele.chrom} for {chromosome}"
                            )
                        key = VariantKey(allele.chrom, allele.pos, allele.ref, allele.alt)
                        collapsed = annotations.annotations_for_allele(key)
                        if not collapsed:
                            qc["vat_unmatched_alt_alleles"] = (
                                int(qc["vat_unmatched_alt_alleles"]) + 1
                            )
                            continue
                        qc["vat_joined_alt_alleles"] = int(qc["vat_joined_alt_alleles"]) + 1
                        for annotation in collapsed:
                            classes = initial_variant_classes(annotation)
                            for sample_id, dosage in zip(
                                allele.carriers,
                                allele.carrier_alt_allele_counts,
                                strict=True,
                            ):
                                writer.writerow(
                                    {
                                        "sample_id": sample_id,
                                        "gene_id": annotation.gene_id,
                                        "gene_symbol": annotation.gene_symbol or "",
                                        "chrom": allele.chrom,
                                        "pos": allele.pos,
                                        "ref": allele.ref,
                                        "alt": allele.alt,
                                        "variant_id": (
                                            f"{allele.chrom}:{allele.pos}:{allele.ref}:{allele.alt}"
                                        ),
                                        "variant_ac": allele.ac,
                                        "variant_af": _format_optional_float(allele.af),
                                        "sample_alt_allele_count": dosage,
                                        "most_severe_consequence": (
                                            annotation.most_severe_consequence or ""
                                        ),
                                        "all_consequences": ",".join(annotation.all_consequences),
                                        "unknown_consequences": ",".join(
                                            annotation.unknown_consequences
                                        ),
                                        "loftee": annotation.loftee or "",
                                        "revel": _format_optional_float(annotation.revel),
                                        "gvs_max_af": _format_optional_float(annotation.gvs_max_af),
                                        "variant_classes": ",".join(classes),
                                    }
                                )
                                qc["carrier_audit_rows"] = int(qc["carrier_audit_rows"]) + 1
    write_json(qc_output, qc)


def _read_vcf_header(path: Path) -> tuple[list[str], dict[str, int], list[str]]:
    contig_lengths: dict[str, int] = {}
    header_lines: list[str] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.startswith("##"):
                header_lines.append(line)
            if line.startswith("##contig=<"):
                identifier = re.search(r"(?:^|,)ID=([^,>]+)", line.removeprefix("##contig=<"))
                length = re.search(r"(?:^|,)length=(\d+)", line.removeprefix("##contig=<"))
                if identifier and length:
                    contig_lengths[identifier.group(1)] = int(length.group(1))
            elif line.startswith("#CHROM"):
                fields = line.split("\t")
                if len(fields) < 9:
                    raise ValueError("VCF header must contain fixed columns and FORMAT")
                header_lines.append("\t".join(fields[:9]))
                if not contig_lengths:
                    raise ValueError("VCF contigs must declare positive lengths")
                return fields[9:], contig_lengths, header_lines
    raise ValueError("VCF does not contain a #CHROM header")


def _stream_tabix_tsv(
    path: Path, region: str, consume: Callable[[Sequence[str]], None]
) -> None:
    for fields in _stream_tabix(path, region):
        consume(fields)


def _stream_tabix_vcf(path: Path, region: str) -> Iterator[list[str]]:
    yield from _stream_tabix(path, region)


def _stream_tabix(path: Path, region: str) -> Iterator[list[str]]:
    command = ["tabix", str(path), region]
    with subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    ) as process:
        assert process.stdout is not None and process.stderr is not None
        for raw_line in process.stdout:
            line = raw_line.rstrip("\r\n")
            if line and not line.startswith("#"):
                yield line.split("\t")
        error_text = process.stderr.read()
        return_code = process.wait()
        if return_code != 0:
            raise subprocess.CalledProcessError(
                return_code, command, stderr=error_text
            )


def _first_nonempty_fields(handle, empty_message: str) -> list[str]:
    for raw_line in handle:
        if raw_line.strip():
            return raw_line.rstrip("\r\n").split("\t")
    raise ValueError(empty_message)


def _open_deterministic_gzip_text(path: Path) -> TextIOWrapper:
    raw = path.open("wb")
    compressed = gzip.GzipFile(fileobj=raw, mode="wb", mtime=0)
    return TextIOWrapper(compressed, encoding="utf-8", newline="")


def _format_optional_float(value: float | None) -> str:
    return "" if value is None else format(value, ".15g")
