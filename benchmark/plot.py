#!/usr/bin/env python3
"""Create runtime and speed-ratio plots from aggregated benchmark CSV files."""

from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def safe_name(value: str) -> str:
    return "".join(character if character.isalnum() else "-" for character in value).strip("-")


def plot_runtime(rows: list[dict[str, str]], plots: Path) -> int:
    grouped: dict[tuple[str, str, str, int], list[float]] = defaultdict(list)
    for row in rows:
        grouped[
            (row["family"], row["operation"], row["engine"], int(row["size"]))
        ].append(float(row["median_seconds"]))

    combinations = sorted({(family, operation) for family, operation, _, _ in grouped})
    for family, operation in combinations:
        fig, axis = plt.subplots(figsize=(7.2, 4.5), constrained_layout=True)
        drawn = False
        for engine, color in (("dtmc", "#3366cc"), ("pydtmc", "#dc3912")):
            points = sorted(
                (
                    size,
                    statistics.median(values),
                )
                for (candidate_family, candidate_operation, candidate_engine, size), values in grouped.items()
                if candidate_family == family
                and candidate_operation == operation
                and candidate_engine == engine
            )
            if points:
                axis.plot(
                    [point[0] for point in points],
                    [point[1] for point in points],
                    marker="o",
                    color=color,
                    label=engine,
                )
                drawn = True
        if not drawn:
            plt.close(fig)
            continue
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("states")
        axis.set_ylabel("median seconds")
        axis.set_title(f"{family}: {operation}")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend()
        fig.savefig(plots / f"{safe_name(family)}--{safe_name(operation)}.svg")
        plt.close(fig)
    return len(combinations)


def plot_ratios(rows: list[dict[str, str]], plots: Path) -> int:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["family"], row["operation"])].append(row)

    for (family, operation), group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: int(row["size"]))
        sizes = [int(row["size"]) for row in ordered]
        ratios = [float(row["geometric_median_pydtmc_over_dtmc"]) for row in ordered]
        lows = [float(row["ci95_low"]) for row in ordered]
        highs = [float(row["ci95_high"]) for row in ordered]

        fig, axis = plt.subplots(figsize=(7.2, 4.5), constrained_layout=True)
        axis.plot(sizes, ratios, marker="o", color="#6a3d9a", label="PyDTMC / dtmc")
        axis.fill_between(sizes, lows, highs, color="#6a3d9a", alpha=0.18, label="95% CI")
        axis.axhline(1.0, color="#555555", linestyle="--", linewidth=1.0)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel("states")
        axis.set_ylabel("median speed ratio")
        axis.set_title(f"{family}: {operation} speed ratio")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend()
        fig.savefig(
            plots / f"{safe_name(family)}--{safe_name(operation)}--speed-ratio.svg"
        )
        plt.close(fig)
    return len(grouped)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    summary_path = RESULTS / "summary.csv"
    ratio_path = RESULTS / "ratio-summary.csv"
    if not summary_path.exists() or not ratio_path.exists():
        raise SystemExit("run benchmark/analyse.py before plotting")

    plots = RESULTS / "plots"
    plots.mkdir(parents=True, exist_ok=True)
    for stale in plots.glob("*.svg"):
        stale.unlink()
    runtime_count = plot_runtime(read_csv(summary_path), plots)
    ratio_count = plot_ratios(read_csv(ratio_path), plots)
    print(f"wrote {runtime_count} runtime and {ratio_count} ratio plots to {plots}")


if __name__ == "__main__":
    main()
