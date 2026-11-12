"""Check Terra file handling without starting a container or a cloud job."""

import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

import pytest
import WDL


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "workflows" / "olink_pc_outliers.wdl"
CHECKER = ROOT / "tools" / "check_workflow_file_writes.py"


def _check_file_writes(path: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(CHECKER), str(path)],
        text=True,
        capture_output=True,
        check=False,
    )


@pytest.mark.parametrize(
    "workflow_body",
    [
        'input { File source = write_lines(["x"]) }',
        'File source = write_lines(["x"])',
        'call Copy { input: source = write_lines(["x"]) }',
        'output { File result = write_lines(["x"]) }',
        'scatter (item in [1]) { if (item > 0) { '
        'call Copy { input: source = write_lines(["x"]) } } }',
        'scatter (item in read_lines(write_lines(["x"]))) { '
        'call Copy { input: source = "unused" } }',
        'if (size(write_lines(["x"])) > 0) { '
        'call Copy { input: source = "unused" } }',
    ],
    ids=["input", "declaration", "call", "output", "nested-call", "scatter", "conditional"],
)
def test_scope_check_rejects_workflow_file_writes(tmp_path, workflow_body):
    fixture = tmp_path / "invalid_terra.wdl"
    fixture.write_text(
        'version 1.0\nworkflow Example {\n' + workflow_body + '\n}\n'
        'task Copy { input { File source } command <<< cat ~{source} >>> }\n'
    )
    result = _check_file_writes(fixture)
    assert result.returncode == 1, result.stderr
    assert "write_lines" in result.stderr
    assert str(fixture) in result.stderr


@pytest.mark.parametrize(
    "expression",
    [
        'write_lines(["x"])',
        'write_tsv([["x"]])',
        'write_map({"x": "y"})',
        'write_json({"x": "y"})',
    ],
)
def test_scope_check_covers_all_file_serializers(tmp_path, expression):
    fixture = tmp_path / "invalid_terra.wdl"
    fixture.write_text(
        'version 1.0\nworkflow Example { File generated = ' + expression + ' }\n'
    )
    result = _check_file_writes(fixture)
    assert result.returncode == 1, result.stderr
    assert expression.split("(")[0] in result.stderr


def test_scope_check_allows_task_scope_and_literal_function_names(tmp_path):
    fixture = tmp_path / "valid_terra.wdl"
    fixture.write_text('''version 1.0
workflow Example {
  # A comment containing write_lines(["example"]) is not an expression.
  String note = "write_json is a function name"
  call Create
  output { File generated = Create.generated }
}
task Create {
  File generated_file = write_lines(["x"])
  command <<< cat "~{generated_file}" > generated.txt >>>
  output { File generated = "generated.txt" }
}
''')
    result = _check_file_writes(fixture)
    assert result.returncode == 0, result.stderr


def test_olink_workflow_has_no_workflow_scope_file_writes():
    result = _check_file_writes(WORKFLOW)
    assert result.returncode == 0, result.stderr


def test_olink_uses_wdl_1_0_and_typed_file_inputs():
    document = WDL.load(str(WORKFLOW))
    assert document.effective_wdl_version == "1.0"
    assert document.workflow is not None
    inputs = {declaration.name: declaration for declaration in document.workflow.inputs}
    assert str(inputs["npx_file"].type) == "File"
    assert {binding.name for binding in document.workflow.required_inputs} == {
        "npx_file", "pc_counts"
    }
    assert str(inputs["docker_image"].type) == "String"
    assert inputs["docker_image"].expr.eval(
        WDL.Env.Bindings(), WDL.StdLib.Base("1.0")
    ).value == "docker.io/epadhi/rarevariantenrichment-olink:root-no-entrypoint"
    assert str(inputs["pc_counts"].type) == "Array[Int]"
    assert inputs["pc_counts"].expr is None
    tasks = {task.name: task for task in document.tasks}
    for task_name, input_name, expected_type in [
        ("PrepareOlink", "npx_file", "File"),
        ("AdjustOlink", "prepared_file", "File"),
        ("CollectOlinkSummaries", "summary_files", "Array[File]"),
        ("CollectOlinkSummaries", "files_to_collect", "Array[File]"),
        ("CollectOlinkSummaries", "pc_plan_file", "File"),
    ]:
        declarations = {declaration.name: declaration for declaration in tasks[task_name].inputs}
        assert str(declarations[input_name].type) == expected_type
    assert set(tasks) == {
        "ValidateOlinkPcGrid", "PrepareOlink", "AdjustOlink", "CollectOlinkSummaries"
    }
    result_members = document.struct_typedefs["OlinkPcResult"].members
    for name in ("z_scores", "protein_outliers", "sample_outliers", "undefined_residual_assays", "summary"):
        assert str(result_members[name]) == "File"
    for task in tasks.values():
        assert {"docker", "cpu", "memory", "disks", "maxRetries", "zones"} <= set(task.runtime)
        assert "Starting" in str(task.command)
        assert "Completed" in str(task.command)


def test_olink_scatters_are_not_nested():
    document = WDL.load(str(WORKFLOW))

    def check_body(body, inside_scatter=False):
        for node in body:
            if isinstance(node, WDL.Tree.Scatter):
                assert not inside_scatter, "Olink must not contain a nested scatter"
                check_body(node.body, inside_scatter=True)
            elif isinstance(node, WDL.Tree.Conditional):
                check_body(node.body, inside_scatter=inside_scatter)

    check_body(document.workflow.body)


@pytest.mark.parametrize("pc_counts", [[5, 0, 2], [0]])
def test_olink_prepares_each_metric_once_then_reuses_it_for_the_pc_grid(pc_counts):
    document = WDL.load(str(WORKFLOW))
    workflow = document.workflow
    scatters = [
        node for node in workflow.body
        if isinstance(node, WDL.Tree.Scatter)
        and any(isinstance(child, WDL.Tree.Call) for child in node.body)
    ]
    assert len(scatters) == 2
    preparation_scatter = next(
        scatter for scatter in scatters
        if any(isinstance(node, WDL.Tree.Call) and node.callee.name == "PrepareOlink"
               for node in scatter.body)
    )
    adjustment_scatter = next(
        scatter for scatter in scatters
        if any(isinstance(node, WDL.Tree.Call) and node.callee.name == "AdjustOlink"
               for node in scatter.body)
    )
    preparation = [node for node in preparation_scatter.body if isinstance(node, WDL.Tree.Call)]
    assert [call.callee.name for call in preparation] == ["PrepareOlink"]
    assert str(preparation[0].inputs["max_pcs"]) == "ValidateOlinkPcGrid.max_pcs"
    adjustment = [node for node in adjustment_scatter.body if isinstance(node, WDL.Tree.Call)]
    assert [call.callee.name for call in adjustment] == ["AdjustOlink"]

    stdlib = WDL.StdLib.Base("1.0")
    environment = WDL.Env.Bindings().bind(
        "pc_counts", WDL.Value.Array(WDL.Type.Int(), [WDL.Value.Int(n) for n in pc_counts])
    )
    for node in workflow.body:
        if isinstance(node, WDL.Tree.Decl):
            environment = environment.bind(node.name, node.expr.eval(environment, stdlib))
    metric_values = preparation_scatter.expr.eval(environment, stdlib)
    prepared_metrics = [
        preparation[0].inputs["metric"].eval(
            environment.bind(preparation_scatter.variable, metric), stdlib
        ).value
        for metric in metric_values.value
    ]
    assert prepared_metrics == ["NPX", "PCNormalizedNPX"]
    state_paths = {
        "NPX": "gs://fixture/npx.rds",
        "PCNormalizedNPX": "gs://fixture/pcnormalized.rds",
    }
    environment = environment.bind(
        "PrepareOlink.prepared_file",
        WDL.Value.Array(WDL.Type.File(), [WDL.Value.File(state_paths[m]) for m in prepared_metrics]),
    )
    actual = []
    for setting in adjustment_scatter.expr.eval(environment, stdlib).value:
        scoped = environment.bind(adjustment_scatter.variable, setting)
        prepared_file = adjustment[0].inputs["prepared_file"].eval(scoped, stdlib)
        assert isinstance(prepared_file, WDL.Value.File)
        actual.append((
            adjustment[0].inputs["metric"].eval(scoped, stdlib).value,
            adjustment[0].inputs["n_pcs"].eval(scoped, stdlib).value,
            prepared_file.value,
        ))
    expected = [
        (metric, count, state_paths[metric])
        for metric in ["NPX", "PCNormalizedNPX"] for count in pc_counts
    ]
    assert actual == expected
    collection = next(
        node for node in workflow.body
        if isinstance(node, WDL.Tree.Call) and node.callee.name == "CollectOlinkSummaries"
    )
    assert str(collection.inputs["summary_files"].type) == "Array[File]"


class _FileStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return filename

    def _devirtualize_filename(self, filename):
        # Model cloud-object size queries without changing the File value.
        return self.cloud_files.get(filename, filename)


def _render_task(task_name, directory, mapping, *, localize):
    document = WDL.load(str(WORKFLOW))
    task = next(task for task in document.tasks if task.name == task_name)
    values = {
        "npx_file": "gs://olink-bucket/npx.tsv",
        "prepared_file": "gs://olink-bucket/prepared.rds",
        "summary_files": ["gs://olink-bucket/summary1.tsv", "gs://olink-bucket/summary2.tsv"],
        "files_to_collect": (
            ["gs://olink-bucket/summary1.tsv"] * 5 + ["gs://olink-bucket/summary2.tsv"] * 5
            + ["gs://olink-bucket/summary1.tsv"] * 4 + ["gs://olink-bucket/summary2.tsv"] * 4
        ),
        "adjustment_paths": [
            [f"{metric}/{stem}.{metric}.PCs0{extension}" for stem, extension in [
                ("protein_z_scores", ".tsv.gz"), ("protein_outliers", ".tsv.gz"),
                ("sample_outliers", ".tsv"), ("undefined_residual_assays", ".tsv"),
                ("summary", ".tsv"),
            ]]
            for metric in ["NPX", "PCNormalizedNPX"]
        ],
        "preparation_paths": [
            [f"{metric}/{stem}.{metric}{extension}" for stem, extension in [
                ("sample_qc", ".tsv"), ("excluded_assays", ".tsv"),
                ("prepare_summary", ".tsv"), ("sessionInfo", ".txt"),
            ]]
            for metric in ["NPX", "PCNormalizedNPX"]
        ],
        "pc_plan_file": "gs://olink-bucket/pc_plan.tsv",
        "pc_counts": [0, 1],
        "metric": "NPX",
        "max_pcs": 1,
        "n_pcs": 1,
        "excluded_assays": ["GBP1", "MAP2K1", "literal ' \" $(touch INJECTION) `touch INJECTION` $HOME"],
        "connectivity_cutoff": -3.0,
        "protein_cutoff": 3.0,
        "block_size": 500,
        "strip_sample_suffix": True,
        "docker_image": "unused",
        "cpu": 1,
        "memory_gb": 1,
        "disk_gb": 1,
        "max_retries": 0,
        "preemptible": 0,
        "zones": "us-central1-a",
    }
    environment = WDL.Env.Bindings()
    stdlib = _FileStdLib("1.0", write_dir=str(directory))
    stdlib.cloud_files = mapping
    for declaration in task.inputs:
        environment = environment.bind(
            declaration.name, WDL.Value.from_json(declaration.type, values[declaration.name])
        )
    # Cromwell evaluates task declarations before command-time localization.
    # If a declaration has changed a File into a String, rewrite_paths cannot
    # repair it. A list written too early retains cloud URIs inside the file.
    for declaration in task.postinputs:
        environment = environment.bind(
            declaration.name,
            declaration.expr.eval(environment, stdlib).coerce(declaration.type),
        )
    if localize:
        environment = environment.map(
            lambda binding: WDL.Env.Binding(
                binding.name,
                WDL.Value.rewrite_paths(
                    binding.value, lambda file: mapping.get(file.value, file.value)
                ),
            )
        )
    return task.command.eval(environment, stdlib).value


def _install_rscript_probe(directory):
    executable = directory / "Rscript"
    executable.write_text('''#!/usr/bin/env python3
import json
import os
from pathlib import Path
import sys

arguments = sys.argv[1:]
if arguments[0] == "--vanilla":
    arguments = arguments[1:]
if arguments[0] == "-":
    stage, arguments = "collect", arguments[1:]
    sys.stdin.read()
else:
    script, stage, *arguments = arguments
    assert script == "/opt/olink/olink_outliers.R", script
assert len(arguments) % 2 == 0, arguments
options = dict(zip(arguments[0::2], arguments[1::2]))
assert "--config" not in options
record = {"stage": stage, "options": options}

def readable(value):
    if value.startswith(("gs://", "s3://", "https://")):
        raise SystemExit("Localization error: unresolved cloud URI " + value)
    assert Path(value).is_file(), value
    with Path(value).open("rb") as handle:
        handle.read(1)

for option in ("--input", "--prepared-file", "--pc-plan-file"):
    if option in options:
        readable(options[option])
for option in ("--summary-files", "--excluded-assays-file", "--files-list", "--relative-paths-list"):
    if option in options:
        readable(options[option])
        lines = Path(options[option]).read_text().splitlines()
        record[option] = lines
        if option in ("--summary-files", "--files-list"):
            for line in lines:
                readable(line)
if "--output-dir" in options:
    Path(options["--output-dir"]).mkdir(exist_ok=True)
Path(os.environ["OLINK_INVOCATION"]).write_text(json.dumps(record))
''')
    executable.chmod(0o755)


def _cloud_mapping(directory):
    localized = directory / "localized ' \" $(touch INJECTION) `touch INJECTION` $HOME"
    localized.mkdir()
    mapping = {}
    for filename in ("npx.tsv", "prepared.rds", "summary1.tsv", "summary2.tsv", "pc_plan.tsv"):
        path = localized / filename
        path.write_text("fixture\n")
        mapping["gs://olink-bucket/" + filename] = str(path)
    return mapping


@pytest.mark.parametrize(
    "task_name,stage",
    [
        ("ValidateOlinkPcGrid", "validate"),
        ("PrepareOlink", "prepare"),
        ("AdjustOlink", "adjust"),
        ("CollectOlinkSummaries", "collect"),
    ],
)
def test_olink_task_commands_use_localized_files_and_preserve_literal_paths(tmp_path, task_name, stage):
    mapping = _cloud_mapping(tmp_path)
    _install_rscript_probe(tmp_path)
    command = _render_task(task_name, tmp_path, mapping, localize=True)
    assert "gs://" not in command
    record_path = tmp_path / "invocation.json"
    result = subprocess.run(
        ["bash", "-c", command],
        cwd=tmp_path,
        env={**os.environ, "PATH": str(tmp_path) + os.pathsep + os.environ["PATH"],
             "OLINK_INVOCATION": str(record_path)},
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "INJECTION").exists()
    record = json.loads(record_path.read_text())
    assert record["stage"] == stage
    if stage == "validate":
        assert record["options"]["--pc-counts"] == "0,1"
    elif stage == "prepare":
        assert record["options"]["--input"] == mapping["gs://olink-bucket/npx.tsv"]
        assert record["--excluded-assays-file"] == [
            "GBP1", "MAP2K1", "literal ' \" $(touch INJECTION) `touch INJECTION` $HOME"
        ]
    elif stage == "adjust":
        assert record["options"]["--prepared-file"] == mapping["gs://olink-bucket/prepared.rds"]
    else:
        assert record["--summary-files"] == [
            mapping["gs://olink-bucket/summary1.tsv"], mapping["gs://olink-bucket/summary2.tsv"]
        ]


@pytest.fixture(scope="module")
def real_rscript():
    executable = shutil.which("Rscript")
    if executable is None:
        pytest.skip("Rscript is required to check the real CLI's localization errors")
    dependency = subprocess.run(
        [executable, "-e", 'quit(status = if (requireNamespace("tidyverse", quietly = TRUE)) 0L else 1L)'],
        text=True, capture_output=True, check=False,
    )
    if dependency.returncode != 0:
        pytest.skip("R tidyverse is required to check the real CLI's localization errors")
    return executable


@pytest.mark.parametrize("task_name", ["PrepareOlink", "AdjustOlink", "CollectOlinkSummaries"])
def test_olink_task_reports_unresolved_cloud_inputs_before_reading(tmp_path, task_name, real_rscript):
    mapping = _cloud_mapping(tmp_path)
    command = _render_task(task_name, tmp_path, mapping, localize=False)
    # Replace only the container script path. Run the production CLI with each
    # rendered argument intact, so the error check does not test the probe.
    if task_name == "CollectOlinkSummaries":
        command = command.replace(
            'source("/opt/olink/olink_outliers.R")',
            "source(" + json.dumps(str(ROOT / "scripts" / "olink_outliers.R")) + ")",
        )
        command = 'Rscript() { command ' + shlex.quote(real_rscript) + ' "$@"; }\n' + command
    else:
        command = (
            "Rscript() { command " + shlex.quote(real_rscript) + " "
            + shlex.quote(str(ROOT / "scripts" / "olink_outliers.R"))
            + ' "${@:2}"; }\n' + command
        )
    result = subprocess.run(
        ["bash", "-c", command],
        cwd=tmp_path,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode != 0
    assert "localization" in result.stderr.lower()
    assert "gs://olink-bucket/" in result.stderr
    assert "Completed" not in result.stderr
