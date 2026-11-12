#!/usr/bin/env python3
"""Run the complete Olink WDL against synthetic data in GitHub Actions."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile


def make_fixture(path: Path) -> None:
    rng = random.Random(23)
    assays = [f"P{i}" for i in range(1, 13)] + ["CONSTANT", "GBP1", "MAP2K1"]
    with gzip.open(path, "wt", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["Assay", "SampleID", "NPX", "PCNormalizedNPX"])
        for sample in range(1, 49):
            for assay in assays:
                regular = rng.gauss(0, 1)
                normalized = regular + rng.gauss(0, 0.5)
                if assay == "CONSTANT":
                    regular = normalized = 5.0
                if sample == 1 and assay == "P1":
                    regular = "NA"
                if sample == 2 and assay == "P1":
                    normalized = "NA"
                writer.writerow([assay, f"{sample:03d}_run", regular, normalized])


def read_table(path: str) -> list[dict[str, str]]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def check_outputs(outputs: dict) -> None:
    prefix = "OlinkPcOutliers."
    results = outputs[prefix + "results"]
    assert len(results) == 4, results
    expected = {(metric, count) for metric in ("NPX", "PCNormalizedNPX") for count in (0, 2)}
    assert {(item["npx_type"], item["pcs_removed"]) for item in results} == expected
    for output, member in {
        "z_score_matrices": "z_scores",
        "protein_outlier_files": "protein_outliers",
        "sample_outlier_files": "sample_outliers",
        "undefined_residual_assay_files": "undefined_residual_assays",
        "run_summaries": "summary",
    }.items():
        assert [Path(path).resolve() for path in outputs[prefix + output]] == [
            Path(item[member]).resolve() for item in results
        ]
    vectors = {}
    for item in results:
        metric, count = item["npx_type"], item["pcs_removed"]
        for name in ("z_scores", "protein_outliers", "sample_outliers", "undefined_residual_assays", "summary"):
            path = item[name]
            assert Path(path).is_file(), path
            assert f".{metric}.PCs{count}." in Path(path).name, path
        matrix = read_table(item["z_scores"])
        assert len(matrix) == 12, matrix
        sample_ids = [name for name in matrix[0] if name != "Assay"]
        assert ("001" not in sample_ids) if metric == "NPX" else ("002" not in sample_ids)
        assert len(sample_ids) >= 3
        vector = []
        for assay in matrix:
            values = [float(assay[sample]) for sample in sample_ids if assay[sample] != "NA"]
            assert values
            assert all(math.isfinite(value) for value in values)
            assert abs(sum(values) / len(values)) < 1e-8
            assert abs(sum(value * value for value in values) / len(values) - 1) < 1e-8
            vector.extend(values)
        vectors[metric, count] = vector
        for call in read_table(item["protein_outliers"]):
            assert call["NPX_type"] == metric
            assert int(call["PCs_removed"]) == count
            assert abs(float(call["Z_score"])) > 3
    for metric in ("NPX", "PCNormalizedNPX"):
        assert vectors[metric, 0] != vectors[metric, 2]
    summary = read_table(outputs[prefix + "summary"])
    assert len(summary) == 4, summary
    assert {(row["NPX_type"], int(row["PCs_removed"])) for row in summary} == expected
    assert len(outputs[prefix + "sample_qc_files"]) == 2
    # miniwdl can expose output symlinks in separate per-variable directories.
    # Resolve the manifest to locate the physical collection directory.
    manifest_path = Path(outputs[prefix + "manifest"]).resolve()
    collection_root = manifest_path.parent
    manifest = read_table(str(manifest_path))
    relative_paths = [row["Relative_path"] for row in manifest]
    assert len(relative_paths) == len(set(relative_paths))
    assert len(manifest) == 30, manifest
    for relative in relative_paths:
        path = Path(relative)
        assert not path.is_absolute() and ".." not in path.parts
        assert (collection_root / path).is_file(), relative
    for item in results:
        for member in ("z_scores", "protein_outliers", "sample_outliers",
                       "undefined_residual_assays", "summary"):
            relative = Path(item[member]).resolve().relative_to(collection_root).as_posix()
            assert relative.startswith(item["npx_type"] + "/")
            assert relative in relative_paths
    collected = {
        Path(path).resolve().relative_to(collection_root).as_posix()
        for path in outputs[prefix + "collected_files"]
    }
    assert collected == set(relative_paths) | {"manifest.tsv"}
    for row in summary:
        for name, relative in row.items():
            if name.endswith("_file"):
                assert relative.startswith(row["NPX_type"] + "/"), relative
                assert relative in relative_paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", required=True)
    args = parser.parse_args()
    runner = shutil.which("miniwdl")
    if runner is None:
        raise SystemExit("Install miniwdl before running this CI smoke test.")
    repo = Path(__file__).resolve().parents[1]
    # GitHub-hosted runners use this filesystem for Docker bind mounts.
    base = Path(os.environ.get("RUNNER_TEMP", tempfile.gettempdir()))
    folder = Path(tempfile.mkdtemp(prefix="olink-wdl-", dir=base))
    data_path = folder / "synthetic_npx.tsv.gz"
    make_fixture(data_path)
    inputs = {
        "OlinkPcOutliers.npx_file": str(data_path),
        "OlinkPcOutliers.pc_counts": [0, 2],
        "OlinkPcOutliers.docker_image": args.image,
        "OlinkPcOutliers.prepare_cpu": 1,
        "OlinkPcOutliers.prepare_memory_gb": 2,
        "OlinkPcOutliers.prepare_disk_gb": 1,
        "OlinkPcOutliers.adjust_cpu": 1,
        "OlinkPcOutliers.adjust_memory_gb": 2,
        "OlinkPcOutliers.adjust_disk_gb": 1,
        "OlinkPcOutliers.max_retries": 0,
    }
    inputs_path = folder / "inputs.json"
    outputs_path = folder / "outputs.json"
    inputs_path.write_text(json.dumps(inputs, indent=2) + "\n")
    command = [runner, "run", str(repo / "workflows/olink_pc_outliers.wdl"),
               "--input", str(inputs_path), "--dir", str(folder / "run"),
               "-o", str(outputs_path), "--as-me", "--no-cache", "--verbose"]
    print(f"[olink-ci] Synthetic workflow files: {folder}", flush=True)
    subprocess.run(command, cwd=repo, check=True)
    result = json.loads(outputs_path.read_text())
    check_outputs(result["outputs"])
    print("[olink-ci] Complete WDL smoke passed: both NPX types at PCs 0 and 2.", flush=True)


if __name__ == "__main__":
    main()
