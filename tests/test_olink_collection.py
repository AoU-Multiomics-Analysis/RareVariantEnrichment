"""Check the final Olink collection without a container or cloud job."""

from __future__ import annotations

import copy
import csv
import gzip
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess

import pytest
import WDL


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "workflows" / "olink_pc_outliers.wdl"
SCRIPT = ROOT / "scripts" / "olink_outliers.R"


class _FileStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return filename

    def _devirtualize_filename(self, filename):
        return self.cloud_files.get(filename, filename)


@pytest.fixture(scope="module")
def real_rscript():
    executable = shutil.which("Rscript")
    if executable is None:
        pytest.skip("Rscript is required for the real collection command")
    dependency = subprocess.run(
        [executable, "--vanilla", "-e",
         'quit(status = if (requireNamespace("tidyverse", quietly = TRUE)) 0L else 1L)'],
        text=True, capture_output=True, check=False,
    )
    if dependency.returncode:
        pytest.skip("R tidyverse is required for the real collection command")
    return executable


@pytest.fixture
def collection_fixture(tmp_path):
    inputs = tmp_path / "inputs ' \" $(touch INJECTION) `touch INJECTION` $HOME"
    inputs.mkdir()
    run = tmp_path / "run"
    run.mkdir()
    mapping = {}
    copies = {}
    files = []
    summaries = []
    adjustment_paths = []
    preparation_paths = []

    def source(relative_path, payload):
        index = len(mapping)
        # The source name must not be used as the collected destination name.
        local = inputs / f"opaque-source-{index:03d}.dat"
        local.write_bytes(payload)
        uri = f"gs://collection-fixture/opaque-source-{index:03d}.dat"
        mapping[uri] = str(local)
        copies[relative_path] = payload
        return uri

    fields = [
        "NPX_type", "PCs_removed", "Protein_z_scores_file",
        "Protein_outliers_file", "Sample_outliers_file",
        "Undefined_residual_assays_file", "Summary_file",
    ]
    for metric in ["NPX", "PCNormalizedNPX"]:
        for count in [10, 2]:
            names = [
                f"protein_z_scores.{metric}.PCs{count}.tsv.gz",
                f"protein_outliers.{metric}.PCs{count}.tsv.gz",
                f"sample_outliers.{metric}.PCs{count}.tsv",
                f"undefined_residual_assays.{metric}.PCs{count}.tsv",
                f"summary.{metric}.PCs{count}.tsv",
            ]
            relative = [f"{metric}/{name}" for name in names]
            adjustment_paths.append(relative)
            payloads = [
                gzip.compress(f"Assay\t001\nP1\t{count}.0\n".encode(), mtime=0),
                gzip.compress(b"Assay\tSampleID\tZ_score\nP1\t001\t4\n", mtime=0),
                b"SampleID\tStatus\n001\tconnectivity_outlier\n",
                b"Assay\tReason\nP2\tundefined\n",
                ("\t".join(fields) + "\n" +
                 "\t".join([metric, str(count), *names]) + "\n").encode(),
            ]
            for path, payload in zip(relative, payloads):
                uri = source(path, payload)
                files.append(uri)
                if path.endswith(names[-1]):
                    summaries.append(uri)
    for metric in ["NPX", "PCNormalizedNPX"]:
        names = [
            f"sample_qc.{metric}.tsv", f"excluded_assays.{metric}.tsv",
            f"prepare_summary.{metric}.tsv", f"sessionInfo.{metric}.txt",
        ]
        relative = [f"{metric}/{name}" for name in names]
        preparation_paths.append(relative)
        for path in relative:
            files.append(source(path, f"fixture for {path}\n".encode()))
    pc_plan = source("pc_plan.tsv", b"PCs_removed\n10\n2\n")
    return {
        "run": run,
        "mapping": mapping,
        "copies": copies,
        "values": {
            "files_to_collect": files,
            "adjustment_paths": adjustment_paths,
            "preparation_paths": preparation_paths,
            "summary_files": summaries,
            "pc_plan_file": pc_plan,
            "docker_image": "unused",
            "max_retries": 0,
            "zones": "us-central1-a",
            "disk_gb": 1,
        },
    }


def _task_environment(fixture, values=None, *, unresolved=None):
    document = WDL.load(str(WORKFLOW))
    task = next(task for task in document.tasks if task.name == "CollectOlinkSummaries")
    values = fixture["values"] if values is None else values
    environment = WDL.Env.Bindings()
    stdlib = _FileStdLib("1.0", write_dir=str(fixture["run"]))
    # Cloud size queries work before localization, including in failure tests.
    stdlib.cloud_files = fixture["mapping"]
    for declaration in task.inputs:
        environment = environment.bind(
            declaration.name,
            WDL.Value.from_json(declaration.type, values[declaration.name]),
        )
    for declaration in task.postinputs:
        environment = environment.bind(
            declaration.name,
            declaration.expr.eval(environment, stdlib).coerce(declaration.type),
        )
    environment = environment.map(
        lambda binding: WDL.Env.Binding(
            binding.name,
            WDL.Value.rewrite_paths(
                binding.value,
                lambda file: (file.value if file.value == unresolved else
                              fixture["mapping"].get(file.value, file.value)),
            ),
        )
    )
    return task, environment, stdlib


def _run_collection(fixture, real_rscript, values=None, *, unresolved=None):
    task, environment, stdlib = _task_environment(fixture, values, unresolved=unresolved)
    command = task.command.eval(environment, stdlib).value
    # Support both the old CLI command and the embedded-R collection command.
    command = command.replace(
        "Rscript /opt/olink/olink_outliers.R",
        f"Rscript {shlex.quote(str(SCRIPT))}",
    )
    command = command.replace('"/opt/olink/olink_outliers.R"', json.dumps(str(SCRIPT)))
    command = command.replace("'/opt/olink/olink_outliers.R'", json.dumps(str(SCRIPT)))
    command = re.sub(r"(?m)^(\s*)Rscript\s", lambda m: m[1] + shlex.quote(real_rscript) + " ", command)
    result = subprocess.run(
        ["bash", "-c", command], cwd=fixture["run"], env=os.environ.copy(),
        text=True, capture_output=True, check=False,
    )
    return result


def _manifest(path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        assert reader.fieldnames == ["Measurement", "PCs_removed", "Output_type", "Relative_path"]
        return list(reader)


def _optional_pc(value):
    return None if value in ("", "NA", "NaN", "null") else int(value)


def test_collection_copies_all_final_files_and_builds_relative_manifest(collection_fixture, real_rscript):
    result = _run_collection(collection_fixture, real_rscript)
    assert result.returncode == 0, result.stderr
    root = collection_fixture["run"] / "collected"
    manifest_path = root / "manifest.tsv"
    assert manifest_path.is_file(), "The collector must create collected/manifest.tsv"
    rows = _manifest(manifest_path)
    expected = set(collection_fixture["copies"]) | {"summary.tsv"}
    assert {row["Relative_path"] for row in rows} == expected
    assert len(rows) == len(expected), "Manifest destination paths must be unique"
    assert "manifest.tsv" not in {row["Relative_path"] for row in rows}
    for relative, payload in collection_fixture["copies"].items():
        assert (root / relative).read_bytes() == payload
    for row in rows:
        relative = row["Relative_path"]
        assert not Path(relative).is_absolute()
        assert ".." not in Path(relative).parts
        assert (root / relative).is_file(), relative
        assert row["Output_type"]
        pc = re.search(r"\.PCs(\d+)\.", relative)
        assert _optional_pc(row["PCs_removed"]) == (int(pc[1]) if pc else None)
        if pc is None:
            assert row["PCs_removed"] == ""
        if "/" in relative:
            assert row["Measurement"] == relative.split("/")[0]
        else:
            assert row["Measurement"] == "shared"
    for metric in ["NPX", "PCNormalizedNPX"]:
        counts = [_optional_pc(row["PCs_removed"]) for row in rows
                  if row["Measurement"] == metric and _optional_pc(row["PCs_removed"]) is not None]
        assert counts == [2] * 5 + [10] * 5
    with (root / "summary.tsv").open(newline="") as stream:
        summaries = list(csv.DictReader(stream, delimiter="\t"))
    assert {(row["NPX_type"], int(row["PCs_removed"])) for row in summaries} == {
        (metric, count) for metric in ["NPX", "PCNormalizedNPX"] for count in [2, 10]
    }
    assert not (collection_fixture["run"] / "INJECTION").exists()


def test_collection_file_outputs_are_known_before_the_command_runs(collection_fixture):
    task, environment, stdlib = _task_environment(collection_fixture)
    assert any(isinstance(declaration.type, WDL.Type.Array) and
               isinstance(declaration.type.item_type, WDL.Type.File)
               for declaration in task.outputs), "Collected File arrays must be declared explicitly"
    output_files = set()
    for declaration in task.outputs:
        value = declaration.expr.eval(environment, stdlib).coerce(declaration.type)
        environment = environment.bind(declaration.name, value)
        if isinstance(value, WDL.Value.File):
            output_files.add(value.value)
        elif isinstance(value, WDL.Value.Array) and isinstance(value.type.item_type, WDL.Type.File):
            assert all(isinstance(item, WDL.Value.File) for item in value.value)
            output_files.update(item.value for item in value.value)
    expected = {"collected/" + name for name in collection_fixture["copies"]}
    expected |= {"collected/summary.tsv", "collected/manifest.tsv"}
    assert expected == output_files
    assert all(path.startswith("collected/") for path in output_files)


@pytest.mark.parametrize("problem", ["count", "duplicate", "absolute", "traversal", "nested-traversal"])
def test_collection_rejects_invalid_destination_paths(collection_fixture, real_rscript, problem):
    values = copy.deepcopy(collection_fixture["values"])
    escape = collection_fixture["run"] / "escaped.tsv.gz"
    if problem == "count":
        values["files_to_collect"].pop()
    elif problem == "duplicate":
        values["adjustment_paths"][0][1] = values["adjustment_paths"][0][0]
    elif problem == "absolute":
        values["adjustment_paths"][0][0] = str(escape)
    elif problem == "traversal":
        values["adjustment_paths"][0][0] = "../escaped.tsv.gz"
    else:
        values["adjustment_paths"][0][0] = "NPX/../../escaped.tsv.gz"
    result = _run_collection(collection_fixture, real_rscript, values)
    assert result.returncode != 0, "Invalid collection paths must stop the task"
    assert result.stderr
    assert not escape.exists()


def test_collection_reports_unresolved_source_with_localized_summaries(collection_fixture, real_rscript):
    unresolved = collection_fixture["values"]["files_to_collect"][0]
    assert unresolved not in collection_fixture["values"]["summary_files"]
    result = _run_collection(collection_fixture, real_rscript, unresolved=unresolved)
    assert result.returncode != 0
    assert "localization error" in result.stderr.lower()
    assert unresolved in result.stderr
