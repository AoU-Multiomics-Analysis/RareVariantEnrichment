import csv
import json
import math
from pathlib import Path
import shutil
import subprocess

import pytest

from rare_variant_enrichment.pc_selection import (
    analyze_carrier_pc_enrichment,
    analyze_lof_pc_enrichment,
    median_log_or_by_pc,
    select_minimum_sufficient_pc_count,
)


def _write_dynamic_results(
    path: Path,
    definitions: list[str],
    *,
    all_unestimable: bool = False,
) -> None:
    rows = [
        [
            "pc_count",
            "z_threshold",
            "carrier_definition",
            "odds_ratio",
            "total_observations",
            "carrier_observations",
        ]
    ]
    for pc_count in (0, 10):
        for definition in definitions:
            for threshold in (-3, -4):
                if all_unestimable or definition == "splice_any":
                    odds_ratio = "NA"
                    carrier_observations = "0"
                else:
                    odds_ratio = str(2.0 + pc_count / 10.0)
                    carrier_observations = "4"
                rows.append(
                    [
                        str(pc_count),
                        str(threshold),
                        definition,
                        odds_ratio,
                        "100",
                        carrier_observations,
                    ]
                )
    path.write_text("\n".join("\t".join(row) for row in rows) + "\n")


def _write_results(path: Path) -> None:
    header = [
        "pc_count",
        "z_threshold",
        "carrier_definition",
        "odds_ratio",
    ]
    rows = [header]
    values = {
        0: {"HC": {"-2": 1.5, "-3": 2.0, "-4": 2.5, "-5": 3.0, "-6": 3.5},
            "any_lof": {"-2": 1.4, "-3": 1.8, "-4": 2.3, "-5": 2.8, "-6": 3.2}},
        10: {"HC": {"-2": 3.0, "-3": 8.0, "-4": 16.0, "-5": 32.0, "-6": 64.0},
             "any_lof": {"-2": 2.8, "-3": 7.0, "-4": 14.0, "-5": 28.0, "-6": 56.0}},
        20: {"HC": {"-2": 3.2, "-3": 8.2, "-4": 16.2, "-5": 32.2, "-6": 64.2},
             "any_lof": {"-2": 3.0, "-3": 7.2, "-4": 14.2, "-5": 28.2, "-6": 56.2}},
    }
    for pc_count in values:
        for definition in ("HC", "any_lof"):
            for threshold in ("-2", "-3", "-4", "-5", "-6"):
                rows.append(
                    [
                        str(pc_count),
                        threshold,
                        definition,
                        str(values[pc_count][definition][threshold]),
                    ]
                )
    path.write_text("\n".join("\t".join(row) for row in rows) + "\n")


def test_median_log_or_excludes_z_minus_two_and_uses_log_scale():
    rows = [
        {"pc_count": "0", "z_threshold": str(z), "carrier_definition": "HC", "odds_ratio": str(or_value)}
        for z, or_value in ((-2, 1.0), (-3, 1.0), (-4, 2.0), (-5, 8.0), (-6, 64.0))
    ]

    summary = median_log_or_by_pc(rows, [-3.0, -4.0, -5.0, -6.0])

    assert summary == [{
        "pc_count": 0,
        "carrier_definition": "HC",
        "included_z_thresholds": [-6.0, -5.0, -4.0, -3.0],
        "median_log_odds_ratio": pytest.approx(math.log(4.0)),
        "median_odds_ratio": pytest.approx(4.0),
    }]


def test_median_log_or_skips_unestimable_zero_observation_rows():
    rows = [
        {"pc_count": "0", "z_threshold": str(z), "carrier_definition": definition, "odds_ratio": "2.0"}
        for definition in ("HC", "any_lof")
        for z in (-3, -4, -5, -6)
    ]
    rows.append(
        {
            "pc_count": "8890",
            "z_threshold": "-3",
            "carrier_definition": "HC",
            "odds_ratio": "NA",
            "total_observations": "0",
        }
    )

    summary = median_log_or_by_pc(rows, [-3.0, -4.0, -5.0, -6.0])

    assert [point["pc_count"] for point in summary] == [0, 0]


def test_median_log_or_rejects_missing_odds_ratio_with_observations():
    rows = [
        {
            "pc_count": "0",
            "z_threshold": str(z),
            "carrier_definition": "HC",
            "odds_ratio": "NA" if z == -3 else "2.0",
            "total_observations": "100",
        }
        for z in (-3, -4, -5, -6)
    ]

    with pytest.raises(ValueError, match="invalid selection fields"):
        median_log_or_by_pc(rows, [-3.0, -4.0, -5.0, -6.0])


def test_selection_uses_earliest_common_pc_within_plateau_fraction():
    rows = [
        {"pc_count": 0, "carrier_definition": "HC", "median_log_odds_ratio": math.log(2.0)},
        {"pc_count": 10, "carrier_definition": "HC", "median_log_odds_ratio": math.log(8.0)},
        {"pc_count": 20, "carrier_definition": "HC", "median_log_odds_ratio": math.log(8.1)},
        {"pc_count": 0, "carrier_definition": "any_lof", "median_log_odds_ratio": math.log(1.8)},
        {"pc_count": 10, "carrier_definition": "any_lof", "median_log_odds_ratio": math.log(7.0)},
        {"pc_count": 20, "carrier_definition": "any_lof", "median_log_odds_ratio": math.log(7.1)},
    ]

    selection = select_minimum_sufficient_pc_count(rows, plateau_fraction=0.95)

    assert selection["selected_pc_count"] == 10
    assert selection["minimum_pc_count_within_plateau_by_definition"] == {
        "HC": 10,
        "any_lof": 10,
    }


def test_analysis_writes_selection_json_and_svg_with_reference_lines(tmp_path: Path):
    results = tmp_path / "results.tsv"
    selection_json = tmp_path / "selection.json"
    plot_svg = tmp_path / "plot.svg"
    _write_results(results)

    analyze_lof_pc_enrichment(
        results,
        selection_json,
        plot_svg,
        selection_z_thresholds=[-3.0, -4.0, -5.0, -6.0],
        plateau_fraction=0.95,
    )

    summary = json.loads(selection_json.read_text())
    assert summary["selection"]["selected_pc_count"] == 10
    assert summary["selection"]["excluded_z_thresholds"] == [-2.0]
    svg = plot_svg.read_text()
    assert 'data-selection-pc="10"' in svg
    assert 'data-median-log-or="HC"' in svg
    assert 'data-median-log-or="any_lof"' in svg
    assert 'height="418"' in svg
    assert 'height="-418"' not in svg


def test_dynamic_selection_excludes_zero_carrier_definition_without_dropping_panel(
    tmp_path: Path,
):
    results = tmp_path / "results.tsv"
    selection_json = tmp_path / "selection.json"
    plot_svg = tmp_path / "plot.svg"
    definitions = ["lof_hc", "missense", "splice_any"]
    _write_dynamic_results(results, definitions)

    analyze_carrier_pc_enrichment(
        results,
        selection_json,
        plot_svg,
        carrier_definitions=definitions,
        selection_z_thresholds=[-3.0, -4.0],
        plateau_fraction=0.95,
    )

    payload = json.loads(selection_json.read_text())
    assert payload["selection"]["excluded_definitions"] == {
        "splice_any": "zero_carriers"
    }
    assert payload["selection"]["estimable_carrier_definitions"] == [
        "lof_hc",
        "missense",
    ]
    assert payload["selection"]["selected_pc_count"] == 10
    svg = plot_svg.read_text()
    assert 'data-carrier-definition="splice_any"' in svg
    assert 'data-exclusion-reason="zero_carriers"' in svg


def test_dynamic_selection_writes_valid_outputs_when_every_definition_is_unestimable(
    tmp_path: Path,
):
    results = tmp_path / "results.tsv"
    selection_json = tmp_path / "selection.json"
    plot_svg = tmp_path / "plot.svg"
    definitions = ["missense", "splice_any"]
    _write_dynamic_results(results, definitions, all_unestimable=True)

    analyze_carrier_pc_enrichment(
        results,
        selection_json,
        plot_svg,
        carrier_definitions=definitions,
        selection_z_thresholds=[-3.0, -4.0],
    )

    payload = json.loads(selection_json.read_text())
    assert payload["selection"]["selected_pc_count"] is None
    assert payload["selection"]["excluded_definitions"] == {
        "missense": "zero_carriers",
        "splice_any": "zero_carriers",
    }
    assert plot_svg.read_text().endswith("</svg>\n")


def test_r_pc_sweep_accepts_ordered_dynamic_definitions(tmp_path: Path):
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is unavailable")
    package_check = subprocess.run(
        [
            rscript,
            "-e",
            'quit(status=ifelse(requireNamespace("tidyverse", quietly=TRUE) && '
            'requireNamespace("ggrepel", quietly=TRUE), 0, 1))',
        ],
        check=False,
    )
    if package_check.returncode != 0:
        pytest.skip("tidyverse or ggrepel is unavailable")

    results = tmp_path / "results.tsv"
    definitions = ["splice_any", "lof_hc", "missense"]
    rows = [["pc_count", "z_threshold", "carrier_definition", "odds_ratio"]]
    for definition_index, definition in enumerate(definitions):
        for pc_count in (0, 10):
            for threshold in (-3, -4):
                rows.append(
                    [
                        str(pc_count),
                        str(threshold),
                        definition,
                        str(2 + definition_index + pc_count / 10),
                    ]
                )
    results.write_text("\n".join("\t".join(row) for row in rows) + "\n")
    summary = tmp_path / "summary.tsv"
    plot = tmp_path / "plot.png"

    completed = subprocess.run(
        [
            rscript,
            "scripts/pc_sweep_qc.R",
            "--results-input",
            str(results),
            "--summary-output",
            str(summary),
            "--plot-output",
            str(plot),
            "--selection-z-thresholds",
            "-3,-4",
            "--carrier-definitions",
            ",".join(definitions),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    with summary.open(newline="") as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert list(dict.fromkeys(row["carrier_definition"] for row in summary_rows)) == (
        definitions
    )
    assert plot.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
