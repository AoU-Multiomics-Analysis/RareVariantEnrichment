import json
import math
from pathlib import Path

import pytest

from rare_variant_enrichment.pc_selection import (
    analyze_lof_pc_enrichment,
    median_log_or_by_pc,
    select_minimum_sufficient_pc_count,
)


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
