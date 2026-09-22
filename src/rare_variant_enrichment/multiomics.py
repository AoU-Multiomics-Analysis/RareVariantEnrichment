"""Match gene/sample IDs and export inclusive multi-omics outlier matrices."""

from contextlib import ExitStack, closing
import csv
import gzip
from itertools import combinations
import logging
import math
import os
from pathlib import Path
import re
import sqlite3
import tempfile
from typing import Sequence

import numpy as np

from rare_variant_enrichment.io import open_text, write_json
from rare_variant_enrichment.haplo_matrix import OUTLIER_Z_THRESHOLD


LOGGER = logging.getLogger(__name__)
MANIFEST_COLUMNS = (
    "matrix_file", "datasets", "z_threshold", "gene_count", "sample_count",
    "outlier_count", "nonoutlier_count", "missing_count", "expression_haplo_required",
)


def require_local_file(path: Path) -> None:
    if re.match(r"^[A-Za-z][A-Za-z0-9+.-]*:/", str(path)):
        raise ValueError(f"Input localization error: unresolved cloud URI: {path}")
    if not path.is_file() or not os.access(path, os.R_OK):
        raise ValueError(f"Input localization error: file is not readable: {path}")


def read_list_file(path: Path) -> list[str]:
    require_local_file(path)
    with open_text(path) as handle:
        values = [line.rstrip("\r\n") for line in handle]
    if any(not value.strip() for value in values):
        raise ValueError(f"List file contains an empty entry: {path}")
    return values


def validate_multiomics_inputs(names: Sequence[str]) -> None:
    if len(names) < 2:
        raise ValueError("At least two datasets are required")
    if len(set(names)) != len(names):
        raise ValueError("Dataset names must be unique")
    if any(not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]{0,63}", name) for name in names):
        raise ValueError("Dataset names must start with a letter and contain only letters, digits, or underscores (up to 64 characters)")


def _store_matrix(connection: sqlite3.Connection, index: int, path: Path, *, binary: bool = False) -> tuple[list[str], int]:
    # Table names come only from enumerated input positions, never file contents.
    table = f"matrix_{index}"
    connection.execute(f"CREATE TABLE {table} (gene TEXT UNIQUE NOT NULL, position INTEGER PRIMARY KEY, scores BLOB)")
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if len(header) < 2 or header[0] != "gene_id":
            raise ValueError(f"Matrix must start with gene_id and sample columns: {path}")
        samples = [sample.strip() for sample in header[1:]]
        if any(not sample for sample in samples):
            raise ValueError(f"Matrix contains an empty sample ID: {path}")
        if len(set(samples)) != len(samples):
            raise ValueError(f"Duplicate sample ID in matrix: {path}")
        count = 0
        for line, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Matrix {path} line {line} has {len(row)} columns; expected {len(header)}")
            gene = row[0].strip()
            if not gene:
                raise ValueError(f"Matrix {path} line {line} has an empty gene ID")
            values = np.empty(len(samples), dtype=np.float64)
            for column, text in enumerate(row[1:]):
                if text.strip().lower() in {"", ".", "na", "nan"}:
                    values[column] = np.nan
                    continue
                try:
                    value = float(text)
                except ValueError as error:
                    raise ValueError(f"Matrix {path} line {line} has a nonnumeric Z score") from error
                if not math.isfinite(value):
                    raise ValueError(f"Matrix {path} line {line} requires finite Z scores or NA")
                if binary and value not in (0.0, 1.0):
                    raise ValueError(f"Haplo matrix {path} line {line} requires 0, 1, or NA")
                values[column] = value
            try:
                connection.execute(f"INSERT INTO {table} VALUES (?, ?, ?)", (gene, count, values.tobytes()))
            except sqlite3.IntegrityError as error:
                raise ValueError(f"Duplicate gene ID {gene!r} in matrix: {path}") from error
            count += 1
            if count % 5000 == 0:
                LOGGER.info("Loaded %d genes from %s", count, path)
    connection.commit()
    return samples, count


def build_outlier_intersections(
    matrix_paths: Sequence[Path], names: Sequence[str],
    output_directory: Path, manifest_output: Path, summary_output: Path,
    *, expression_haplo_path: Path | None = None,
) -> None:
    validate_multiomics_inputs(names)
    if len(matrix_paths) != len(names):
        raise ValueError("Matrix file count must match dataset name count")
    expression_index = names.index("expression") if "expression" in names else None
    if expression_index is not None and expression_haplo_path is None:
        raise ValueError("An expression dataset requires an expression haplo matrix")
    if expression_index is None and expression_haplo_path is not None:
        raise ValueError("Cannot supply an expression haplo matrix without an expression dataset")
    for path in matrix_paths:
        require_local_file(path)
    if expression_haplo_path is not None:
        require_local_file(expression_haplo_path)
    if output_directory.exists() and any(output_directory.iterdir()):
        raise ValueError("Intersection output directory must be empty")
    output_directory.mkdir(parents=True, exist_ok=True)
    expected_count = 2 ** len(names) - len(names) - 1
    if expression_index is not None:
        expected_count += 2 ** (len(names) - 1) - 1
    LOGGER.info("Starting intersections: datasets=%d Z<=%g output_matrices=%d", len(names), OUTLIER_Z_THRESHOLD, expected_count)
    entries = []
    dataset_qc = []
    with tempfile.TemporaryDirectory(prefix="multiomics-", dir=output_directory.parent) as temporary:
        with closing(sqlite3.connect(str(Path(temporary) / "matrices.sqlite"))) as connection:
            sample_lists = []
            for index, path in enumerate(matrix_paths):
                samples, genes = _store_matrix(connection, index, path)
                sample_lists.append(samples)
                dataset_qc.append({"name": names[index], "gene_count": genes, "sample_count": len(samples)})
            haplo_table = f"matrix_{len(names)}"
            haplo_sample_indexes = {}
            if expression_haplo_path is not None:
                haplo_samples, haplo_genes = _store_matrix(connection, len(names), expression_haplo_path, binary=True)
                if set(haplo_samples) != set(sample_lists[expression_index]):
                    raise ValueError("Expression haplo sample IDs must match expression Z-score sample IDs")
                expression_table = f"matrix_{expression_index}"
                missing_gene = connection.execute(
                    f"SELECT 1 FROM {expression_table} e WHERE NOT EXISTS "
                    f"(SELECT 1 FROM {haplo_table} h WHERE h.gene = e.gene) LIMIT 1"
                ).fetchone()
                if haplo_genes != dataset_qc[expression_index]["gene_count"] or missing_gene is not None:
                    raise ValueError("Expression haplo gene IDs must match expression Z-score gene IDs")
                haplo_sample_indexes = {sample: index for index, sample in enumerate(haplo_samples)}
            sample_indexes = [{name: index for index, name in enumerate(samples)} for samples in sample_lists]
            for group_size in range(2, len(names) + 1):
                for group in combinations(range(len(names)), group_size):
                    group_name = "-".join(names[index] for index in group)
                    first = group[0]
                    samples = [sample for sample in sample_lists[first] if all(sample in sample_indexes[index] for index in group[1:])]
                    columns = [np.asarray([sample_indexes[index][sample] for sample in samples], dtype=int) for index in group]
                    has_expression = expression_index in group
                    haplo_columns = np.asarray([haplo_sample_indexes[sample] for sample in samples], dtype=int) if has_expression else None
                    tables = [f"matrix_{index}" for index in group]
                    if has_expression:
                        tables.append(haplo_table)
                    query = "SELECT " + tables[0] + ".gene, " + ", ".join(table + ".scores" for table in tables)
                    query += " FROM " + tables[0] + " " + " ".join("INNER JOIN " + table + " USING (gene)" for table in tables[1:])
                    query += " ORDER BY " + tables[0] + ".position"
                    group_entries = []
                    with ExitStack() as stack:
                        writers = []
                        modes = [False, True] if has_expression else [False]
                        for haplo_required in modes:
                            rule = "z_le_minus3.expression_haplo" if haplo_required else "z_le_minus3"
                            filename = f"{group_name}.{rule}.tsv.gz"
                            handle = stack.enter_context(gzip.open(output_directory / filename, "wt", encoding="utf-8", newline=""))
                            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                            writer.writerow(["gene_id", *samples])
                            writers.append(writer)
                            group_entries.append(dict(zip(MANIFEST_COLUMNS, [filename, ",".join(names[index] for index in group), OUTLIER_Z_THRESHOLD, 0, len(samples), 0, 0, 0, haplo_required])))
                        for row in connection.execute(query):
                            values = np.stack([np.frombuffer(blob, dtype=np.float64)[indexes] for blob, indexes in zip(row[1:1 + len(group)], columns)])
                            z_complete = np.all(np.isfinite(values), axis=0)
                            z_pass = np.max(values, axis=0) <= OUTLIER_Z_THRESHOLD
                            haplo = np.frombuffer(row[-1], dtype=np.float64)[haplo_columns] if has_expression else None
                            for haplo_required, writer, entry in zip(modes, writers, group_entries):
                                complete = z_complete.copy()
                                outlier = z_pass.copy()
                                if haplo_required:
                                    complete &= np.isfinite(haplo)
                                    outlier &= haplo == 1
                                outlier &= complete
                                writer.writerow([row[0], *np.where(complete, np.where(outlier, "1", "0"), "NA")])
                                entry["gene_count"] += 1
                                entry["outlier_count"] += int(np.count_nonzero(outlier))
                                entry["nonoutlier_count"] += int(np.count_nonzero(complete & ~outlier))
                                entry["missing_count"] += int(np.count_nonzero(~complete))
                    entries.extend(group_entries)
                    LOGGER.info("Completed intersection %s: genes=%d samples=%d", group_entries[0]["datasets"], group_entries[0]["gene_count"], len(samples))
    with manifest_output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(entries)
    write_json(summary_output, {
        "datasets": dataset_qc, "z_threshold": OUTLIER_Z_THRESHOLD, "z_comparison": "<=", "matrix_count": len(entries),
        "outlier_rule": "all participating datasets have finite Z <= -3; expression combinations have Z-only and Z-plus-haplo versions; only the latter requires expression haplo = 1",
        "intersection_type": "inclusive; datasets outside each combination are ignored",
        "missing_rule": "NA when any participating Z score or required expression haplo call is missing",
        "alignment": "shared gene and sample IDs per combination; first participating dataset order",
        "intersections": entries,
    })
    LOGGER.info("Finished intersections: wrote %d matrices", len(entries))
