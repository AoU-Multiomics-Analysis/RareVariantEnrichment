"""PC-count selection summaries and dependency-free SVG plots."""

from collections import defaultdict
import csv
from dataclasses import dataclass
import html
import math
from pathlib import Path
import statistics
from typing import Mapping, Sequence

from rare_variant_enrichment.io import open_text, write_json


DEFAULT_SELECTION_Z_THRESHOLDS = (-3.0, -4.0, -5.0, -6.0)
DEFAULT_PLATEAU_FRACTION = 0.95
PLOT_DEFINITIONS = ("HC", "any_lof")
PLOT_COLORS = {
    -2.0: "#2563eb",
    -3.0: "#ea580c",
    -4.0: "#16a34a",
    -5.0: "#9333ea",
    -6.0: "#dc2626",
}
REQUIRED_RESULT_COLUMNS = ("pc_count", "z_threshold", "carrier_definition", "odds_ratio")


@dataclass(frozen=True)
class MedianLogOrPoint:
    pc_count: int
    carrier_definition: str
    included_z_thresholds: tuple[float, ...]
    median_log_odds_ratio: float
    median_odds_ratio: float

    def as_dict(self) -> dict[str, object]:
        return {
            "pc_count": self.pc_count,
            "carrier_definition": self.carrier_definition,
            "included_z_thresholds": list(self.included_z_thresholds),
            "median_log_odds_ratio": self.median_log_odds_ratio,
            "median_odds_ratio": self.median_odds_ratio,
        }


def read_enrichment_rows(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or any(
            column not in reader.fieldnames for column in REQUIRED_RESULT_COLUMNS
        ):
            raise ValueError(
                "Enrichment results must contain columns: "
                + ", ".join(REQUIRED_RESULT_COLUMNS)
            )
        return [row for row in reader if row]


def median_log_or_by_pc(
    rows: Sequence[Mapping[str, object]],
    selection_z_thresholds: Sequence[float],
    carrier_definitions: Sequence[str] = PLOT_DEFINITIONS,
) -> list[dict[str, object]]:
    thresholds = _validate_selection_thresholds(selection_z_thresholds)
    definitions = _validate_definitions(carrier_definitions)
    grouped: dict[tuple[int, str], dict[float, float]] = defaultdict(dict)
    for row in rows:
        try:
            pc_count = int(row["pc_count"])
            threshold = float(row["z_threshold"])
            definition = str(row["carrier_definition"])
            odds_ratio = float(row["odds_ratio"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Enrichment result row has invalid selection fields") from error
        if definition not in definitions or threshold not in thresholds:
            continue
        if pc_count < 0 or not math.isfinite(odds_ratio) or odds_ratio <= 0:
            continue
        grouped[(pc_count, definition)][threshold] = odds_ratio

    points: list[MedianLogOrPoint] = []
    for (pc_count, definition), values in grouped.items():
        if set(values) != set(thresholds):
            continue
        included = tuple(sorted(values))
        log_odds = [math.log(values[threshold]) for threshold in included]
        median_log_odds_ratio = float(statistics.median(log_odds))
        points.append(
            MedianLogOrPoint(
                pc_count,
                definition,
                included,
                median_log_odds_ratio,
                math.exp(median_log_odds_ratio),
            )
        )
    points.sort(key=lambda point: (point.pc_count, definitions.index(point.carrier_definition)))
    return [point.as_dict() for point in points]


def select_minimum_sufficient_pc_count(
    median_points: Sequence[Mapping[str, object]],
    plateau_fraction: float = DEFAULT_PLATEAU_FRACTION,
    carrier_definitions: Sequence[str] = PLOT_DEFINITIONS,
) -> dict[str, object]:
    if not math.isfinite(plateau_fraction) or not 0 < plateau_fraction <= 1:
        raise ValueError("plateau_fraction must be finite and in (0, 1]")
    definitions = _validate_definitions(carrier_definitions)
    points_by_definition: dict[str, list[Mapping[str, object]]] = {
        definition: [] for definition in definitions
    }
    for point in median_points:
        definition = str(point.get("carrier_definition", ""))
        if definition in points_by_definition:
            points_by_definition[definition].append(point)
    if any(not points for points in points_by_definition.values()):
        missing = [definition for definition, points in points_by_definition.items() if not points]
        raise ValueError("No median logOR points for carrier definition(s): " + ", ".join(missing))

    max_log_by_definition: dict[str, float] = {}
    threshold_by_definition: dict[str, float] = {}
    first_pc_by_definition: dict[str, int] = {}
    for definition, points in points_by_definition.items():
        ordered = sorted(points, key=lambda point: int(point["pc_count"]))
        max_log = max(float(point["median_log_odds_ratio"]) for point in ordered)
        threshold = max_log + math.log(plateau_fraction)
        qualifying = [
            int(point["pc_count"])
            for point in ordered
            if float(point["median_log_odds_ratio"]) >= threshold
        ]
        if not qualifying:
            raise ValueError(f"No PC count reaches the plateau for {definition}")
        max_log_by_definition[definition] = max_log
        threshold_by_definition[definition] = threshold
        first_pc_by_definition[definition] = min(qualifying)

    return {
        "carrier_definitions": list(definitions),
        "plateau_fraction": plateau_fraction,
        "max_median_log_odds_ratio_by_definition": max_log_by_definition,
        "plateau_log_odds_ratio_threshold_by_definition": threshold_by_definition,
        "minimum_pc_count_within_plateau_by_definition": first_pc_by_definition,
        "selected_pc_count": max(first_pc_by_definition.values()),
    }


def analyze_lof_pc_enrichment(
    results_input: Path,
    selection_output: Path,
    plot_output: Path,
    *,
    selection_z_thresholds: Sequence[float] = DEFAULT_SELECTION_Z_THRESHOLDS,
    plateau_fraction: float = DEFAULT_PLATEAU_FRACTION,
) -> None:
    rows = read_enrichment_rows(results_input)
    median_points = median_log_or_by_pc(rows, selection_z_thresholds)
    selection = select_minimum_sufficient_pc_count(median_points, plateau_fraction)
    observed_thresholds = sorted(
        {
            float(row["z_threshold"])
            for row in rows
            if row.get("z_threshold") not in (None, "")
        }
    )
    selection["included_z_thresholds"] = list(_validate_selection_thresholds(selection_z_thresholds))
    selection["excluded_z_thresholds"] = [
        threshold
        for threshold in observed_thresholds
        if threshold not in selection["included_z_thresholds"]
    ]
    payload = {
        "selection": selection,
        "median_log_or": median_points,
    }
    write_json(selection_output, payload)
    write_lof_pc_svg(plot_output, rows, median_points, selection)


def write_lof_pc_svg(
    path: Path,
    rows: Sequence[Mapping[str, object]],
    median_points: Sequence[Mapping[str, object]],
    selection: Mapping[str, object],
) -> None:
    definitions = list(selection["carrier_definitions"])
    thresholds = sorted(
        {
            float(row["z_threshold"])
            for row in rows
            if row.get("z_threshold") not in (None, "")
        },
        reverse=True,
    )
    plot_rows = [
        row
        for row in rows
        if row.get("carrier_definition") in definitions
        and _positive_finite(row.get("odds_ratio"))
    ]
    if not plot_rows:
        raise ValueError("No finite positive odds ratios available for the plot")
    x_values = [int(row["pc_count"]) for row in plot_rows]
    y_values = [float(row["odds_ratio"]) for row in plot_rows]
    y_values.extend(math.exp(float(point["median_log_odds_ratio"])) for point in median_points)
    x_min, x_max = min(x_values), max(x_values)
    y_min = max(0.75, min(y_values) / 1.18)
    y_max = max(y_values) * 1.22
    width, height = 1400, 620
    margin_left, margin_right, gap = 90, 35, 70
    panel_width = (width - margin_left - margin_right - gap) / 2
    plot_top, plot_bottom = 92, 510
    plot_width, plot_height = panel_width - 10, plot_top - plot_bottom
    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-label="LoF enrichment across principal components">',
        "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#111827} .grid{stroke:#d1d5db;stroke-width:1} .axis{stroke:#6b7280;stroke-width:1} .reference{stroke:#6b7280;stroke-width:1.3;stroke-dasharray:5 4} .median{stroke:#111827;stroke-width:3;fill:none} .threshold{fill:none;stroke-width:2} .selection{stroke:#111827;stroke-width:2;stroke-dasharray:7 4} .definition-selection{stroke:#6b7280;stroke-width:1.5;stroke-dasharray:4 4}</style>",
    ]
    legend_x = 460
    for threshold in thresholds:
        color = PLOT_COLORS.get(threshold, "#4b5563")
        svg.append(f'<line x1="{legend_x}" y1="31" x2="{legend_x + 22}" y2="31" stroke="{color}" stroke-width="3"/>')
        svg.append(f'<text x="{legend_x + 30}" y="36" font-size="13">z &lt;= {threshold:g}</text>')
        legend_x += 105
    svg.append('<line x1="995" y1="31" x2="1017" y2="31" stroke="#111827" stroke-width="4"/>')
    svg.append('<text x="1025" y="36" font-size="13">median logOR</text>')

    for index, definition in enumerate(definitions):
        x0 = margin_left + index * (panel_width + gap)
        svg.extend(_svg_panel(
            x0,
            definition,
            rows,
            median_points,
            selection,
            thresholds,
            x_min,
            x_max,
            y_min,
            y_max,
            plot_top,
            plot_bottom,
            plot_width,
            plot_height,
        ))
    svg.append("</svg>\n")
    path.write_text("\n".join(svg), encoding="utf-8")


def _svg_panel(
    x0: float,
    definition: str,
    rows: Sequence[Mapping[str, object]],
    median_points: Sequence[Mapping[str, object]],
    selection: Mapping[str, object],
    thresholds: Sequence[float],
    x_min: int,
    x_max: int,
    y_min: float,
    y_max: float,
    plot_top: float,
    plot_bottom: float,
    plot_width: float,
    plot_height: float,
) -> list[str]:
    def sx(value: float) -> float:
        return x0 + (value - x_min) / max(1, x_max - x_min) * plot_width

    log_min, log_max = math.log10(y_min), math.log10(y_max)

    def sy(value: float) -> float:
        return plot_bottom - (math.log10(value) - log_min) / (log_max - log_min) * plot_height

    output = [f'<text x="{x0:g}" y="70" font-size="16">{html.escape(_definition_label(definition))}</text>']
    for tick in (1, 2, 5, 10, 20, 50, 100, 200):
        if y_min <= tick <= y_max:
            y = sy(tick)
            output.append(f'<line class="grid" x1="{x0:g}" y1="{y:g}" x2="{x0 + plot_width:g}" y2="{y:g}"/>')
            output.append(f'<text x="{x0 - 12:g}" y="{y + 4:g}" text-anchor="end" font-size="12">{tick:g}</text>')
    output.append(f'<rect x="{x0:g}" y="{plot_top:g}" width="{plot_width:g}" height="{plot_height:g}" fill="none" class="axis"/>')
    output.append(f'<line class="reference" x1="{x0:g}" y1="{sy(1):g}" x2="{x0 + plot_width:g}" y2="{sy(1):g}"/>')
    for tick in _x_ticks(x_min, x_max):
        x = sx(tick)
        output.append(f'<line class="axis" x1="{x:g}" y1="{plot_bottom:g}" x2="{x:g}" y2="{plot_bottom + 5:g}"/>')
        output.append(f'<text x="{x:g}" y="{plot_bottom + 23:g}" text-anchor="middle" font-size="12">{tick:,}</text>')
    output.append(f'<text x="{x0 + plot_width / 2:g}" y="{plot_bottom + 53:g}" text-anchor="middle" font-size="13">Number of PCs included (k)</text>')
    output.append(f'<text x="{x0 + plot_width / 2:g}" y="{plot_top - 22:g}" text-anchor="middle" font-size="13">Odds ratio (log scale)</text>')

    for threshold in thresholds:
        values = sorted(
            (
                (int(row["pc_count"]), float(row["odds_ratio"]))
                for row in rows
                if row.get("carrier_definition") == definition
                and float(row["z_threshold"]) == threshold
                and _positive_finite(row.get("odds_ratio"))
            ),
            key=lambda value: value[0],
        )
        if values:
            path_data = " ".join(
                ("M" if index == 0 else "L") + f" {sx(pc):g},{sy(or_value):g}"
                for index, (pc, or_value) in enumerate(values)
            )
            color = PLOT_COLORS.get(threshold, "#4b5563")
            output.append(f'<path class="threshold" data-z-threshold="{threshold:g}" d="{path_data}" stroke="{color}"/>')

    median_values = sorted(
        (
            int(point["pc_count"]),
            math.exp(float(point["median_log_odds_ratio"])),
        )
        for point in median_points
        if point.get("carrier_definition") == definition
    )
    if median_values:
        path_data = " ".join(
            ("M" if index == 0 else "L") + f" {sx(pc):g},{sy(or_value):g}"
            for index, (pc, or_value) in enumerate(median_values)
        )
        output.append(f'<path class="median" data-median-log-or="{html.escape(definition)}" d="{path_data}"/>')

    per_definition = selection["minimum_pc_count_within_plateau_by_definition"]
    definition_pc = int(per_definition[definition])
    common_pc = int(selection["selected_pc_count"])
    x = sx(definition_pc)
    output.append(f'<line class="definition-selection" data-definition-selection-pc="{definition_pc}" x1="{x:g}" y1="{plot_top:g}" x2="{x:g}" y2="{plot_bottom:g}"/>')
    output.append(f'<text x="{x + 5:g}" y="{plot_top + 15:g}" font-size="11">plateau K={definition_pc:,}</text>')
    x = sx(common_pc)
    output.append(f'<line class="selection" data-selection-pc="{common_pc}" x1="{x:g}" y1="{plot_top:g}" x2="{x:g}" y2="{plot_bottom:g}"/>')
    output.append(f'<text x="{x + 5:g}" y="{plot_top + 31:g}" font-size="11">selected K={common_pc:,}</text>')
    return output


def _validate_selection_thresholds(thresholds: Sequence[float]) -> tuple[float, ...]:
    values = tuple(float(value) for value in thresholds)
    if not values or any(not math.isfinite(value) or value >= 0 for value in values):
        raise ValueError("selection_z_thresholds must be finite, unique, and negative")
    if len(values) != len(set(values)):
        raise ValueError("selection_z_thresholds must be finite, unique, and negative")
    return values


def _validate_definitions(definitions: Sequence[str]) -> tuple[str, ...]:
    values = tuple(str(value) for value in definitions)
    if not values or len(values) != len(set(values)):
        raise ValueError("carrier_definitions must be non-empty and unique")
    return values


def _positive_finite(value: object) -> bool:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(parsed) and parsed > 0


def _definition_label(definition: str) -> str:
    return "Any LoF / HC_or_LC" if definition == "any_lof" else definition


def _x_ticks(x_min: int, x_max: int) -> list[int]:
    step = 1000
    return list(range((x_min // step) * step, x_max + step, step))
