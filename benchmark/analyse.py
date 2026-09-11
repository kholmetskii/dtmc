#!/usr/bin/env python3
"""Aggregate Criterion and pyperf JSON into comparable robust statistics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def criterion_samples(path: Path) -> Iterable[tuple[str, float]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if (
        isinstance(payload, list)
        and len(payload) == 3
        and payload[0] == "criterion"
        and isinstance(payload[2], list)
    ):
        reports = payload[2]
    else:
        reports = payload if isinstance(payload, list) else payload.get("reports", [])
    for report in reports:
        name = report.get("reportName") or report.get("name")
        measured = report.get("reportMeasured", [])
        samples: list[float] = []
        keys = report.get("reportKeys", [])
        for measurement in measured:
            if isinstance(measurement, dict):
                elapsed = measurement.get("measTime", measurement.get("time"))
                iterations = measurement.get("measIters", measurement.get("iters", 1))
                if elapsed is not None:
                    samples.append(float(elapsed) / max(1.0, float(iterations)))
            elif isinstance(measurement, list) and keys:
                values = dict(zip(keys, measurement))
                elapsed = values.get("time")
                iterations = values.get("iters", 1)
                if elapsed is not None:
                    samples.append(float(elapsed) / max(1.0, float(iterations)))
        if not samples:
            analysis = report.get("reportAnalysis", report.get("analysis", {}))
            estimate = analysis.get("anMean", analysis.get("mean", {}))
            point = estimate.get("estPoint", estimate.get("point")) if isinstance(estimate, dict) else None
            if point is not None:
                samples = [float(point)]
        for sample in samples:
            yield str(name), sample


def pyperf_samples(path: Path) -> Iterable[tuple[str, float]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    for benchmark in payload.get("benchmarks", []):
        metadata = benchmark.get("metadata", {})
        name = metadata.get("name") or benchmark.get("name")
        for run in benchmark.get("runs", []):
            for value in run.get("values", []):
                yield str(name), float(value)


def quartiles(values: list[float]) -> tuple[float, float]:
    if len(values) < 2:
        return values[0], values[0]
    cuts = statistics.quantiles(values, n=4, method="inclusive")
    return cuts[0], cuts[2]


def percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def bootstrap_median_interval(values: list[float], seed: int) -> tuple[float, float]:
    if len(values) == 1:
        return values[0], values[0]
    rng = random.Random(seed)
    estimates = [
        statistics.median(rng.choices(values, k=len(values)))
        for _ in range(2000)
    ]
    return percentile(estimates, 0.025), percentile(estimates, 0.975)


def stable_seed(value: str) -> int:
    return int.from_bytes(hashlib.sha256(value.encode("utf-8")).digest()[:8], "big")


def summarize(values: list[float], seed_key: str) -> dict[str, float | int]:
    median = statistics.median(values)
    deviations = [abs(value - median) for value in values]
    q1, q3 = quartiles(values)
    ci_low, ci_high = bootstrap_median_interval(values, stable_seed(seed_key))
    return {
        "median_seconds": median,
        "mad_seconds": statistics.median(deviations),
        "q1_seconds": q1,
        "q3_seconds": q3,
        "median_ci95_low_seconds": ci_low,
        "median_ci95_high_seconds": ci_high,
        "samples": len(values),
    }


def split_name(name: str, metadata: dict[str, dict[str, Any]]) -> tuple[str, str]:
    for dataset in sorted(metadata, key=len, reverse=True):
        prefix = dataset + "/"
        if name.startswith(prefix):
            return dataset, name[len(prefix) :]
    raise ValueError(f"benchmark name does not begin with a manifest dataset id: {name}")


def active_result_paths(
    results: Path, engine: str, datasets: Iterable[str]
) -> Iterable[Path]:
    directory = results / "raw" / engine
    for dataset in sorted(datasets):
        path = directory / f"{dataset}.json"
        if path.is_file():
            yield path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root", type=Path, default=ROOT / "data" / "generated"
    )
    parser.add_argument("--results-root", type=Path, default=RESULTS)
    args = parser.parse_args()
    data_root = args.data_root.resolve()
    results = args.results_root.resolve()

    manifest = json.loads(
        (data_root / "manifest.json").read_text(encoding="utf-8")
    )
    metadata = {entry["id"]: entry for entry in manifest["datasets"]}
    grouped: dict[tuple[str, str, str], list[float]] = defaultdict(list)

    for path in active_result_paths(results, "haskell", metadata):
        for name, value in criterion_samples(path):
            dataset, operation = split_name(name, metadata)
            grouped[("dtmc", dataset, operation)].append(value)
    for path in active_result_paths(results, "python", metadata):
        for name, value in pyperf_samples(path):
            dataset, operation = split_name(name, metadata)
            grouped[("pydtmc", dataset, operation)].append(value)

    rows: list[dict[str, Any]] = []
    for (engine, dataset, operation), values in sorted(grouped.items()):
        entry = metadata[dataset]
        rows.append(
            {
                "engine": engine,
                "dataset": dataset,
                "family": entry["family"],
                "size": entry["size"],
                "seed": entry["seed"],
                "operation": operation,
                **summarize(values, f"{engine}/{dataset}/{operation}"),
            }
        )

    results.mkdir(parents=True, exist_ok=True)
    fields = [
        "engine",
        "dataset",
        "family",
        "size",
        "seed",
        "operation",
        "median_seconds",
        "mad_seconds",
        "q1_seconds",
        "q3_seconds",
        "median_ci95_low_seconds",
        "median_ci95_high_seconds",
        "samples",
    ]
    with (results / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    by_key = {(row["engine"], row["dataset"], row["operation"]): row for row in rows}
    ratios: list[dict[str, Any]] = []
    for row in rows:
        if row["engine"] != "dtmc":
            continue
        peer = by_key.get(("pydtmc", row["dataset"], row["operation"]))
        if peer is None:
            continue
        ratios.append(
            {
                "dataset": row["dataset"],
                "family": row["family"],
                "size": row["size"],
                "seed": row["seed"],
                "operation": row["operation"],
                "pydtmc_over_dtmc": peer["median_seconds"] / row["median_seconds"],
            }
        )
    ratio_fields = ["dataset", "family", "size", "seed", "operation", "pydtmc_over_dtmc"]
    with (results / "ratios.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=ratio_fields)
        writer.writeheader()
        writer.writerows(ratios)

    ratio_groups: dict[tuple[str, int, str], list[float]] = defaultdict(list)
    for row in ratios:
        ratio_groups[(row["family"], int(row["size"]), row["operation"])].append(
            float(row["pydtmc_over_dtmc"])
        )
    ratio_summary: list[dict[str, Any]] = []
    for (family, size, operation), values in sorted(ratio_groups.items()):
        logs = [math.log(value) for value in values]
        log_low, log_high = bootstrap_median_interval(
            logs, stable_seed(f"ratio/{family}/{size}/{operation}")
        )
        ratio_summary.append(
            {
                "family": family,
                "size": size,
                "operation": operation,
                "geometric_median_pydtmc_over_dtmc": math.exp(statistics.median(logs)),
                "ci95_low": math.exp(log_low),
                "ci95_high": math.exp(log_high),
                "datasets": len(values),
            }
        )
    ratio_summary_fields = [
        "family",
        "size",
        "operation",
        "geometric_median_pydtmc_over_dtmc",
        "ci95_low",
        "ci95_high",
        "datasets",
    ]
    with (results / "ratio-summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=ratio_summary_fields)
        writer.writeheader()
        writer.writerows(ratio_summary)

    scaling_points: dict[tuple[str, str, str, int], list[float]] = defaultdict(list)
    for row in rows:
        scaling_points[
            (row["engine"], row["family"], row["operation"], int(row["size"]))
        ].append(float(row["median_seconds"]))
    scaling_groups: dict[tuple[str, str, str], list[tuple[int, float]]] = defaultdict(list)
    for (engine, family, operation, size), values in scaling_points.items():
        scaling_groups[(engine, family, operation)].append(
            (size, statistics.median(values))
        )
    scaling_rows: list[dict[str, Any]] = []
    for (engine, family, operation), points in sorted(scaling_groups.items()):
        if len(points) < 3:
            continue
        points.sort()
        xs = [math.log(float(size)) for size, _ in points]
        ys = [math.log(seconds) for _, seconds in points]
        x_mean = statistics.mean(xs)
        y_mean = statistics.mean(ys)
        denominator = sum((value - x_mean) ** 2 for value in xs)
        slope = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys)) / denominator
        intercept = y_mean - slope * x_mean
        fitted = [intercept + slope * value for value in xs]
        total = sum((value - y_mean) ** 2 for value in ys)
        residual = sum((actual - predicted) ** 2 for actual, predicted in zip(ys, fitted))
        scaling_rows.append(
            {
                "engine": engine,
                "family": family,
                "operation": operation,
                "slope": slope,
                "r_squared": 1.0 - residual / total if total > 0 else 1.0,
                "minimum_size": points[0][0],
                "maximum_size": points[-1][0],
                "points": len(points),
            }
        )
    scaling_fields = [
        "engine",
        "family",
        "operation",
        "slope",
        "r_squared",
        "minimum_size",
        "maximum_size",
        "points",
    ]
    with (results / "scaling.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=scaling_fields)
        writer.writeheader()
        writer.writerows(scaling_rows)

    print(
        f"wrote {len(rows)} summaries, {len(ratios)} paired ratios, "
        f"and {len(scaling_rows)} scaling fits"
    )


if __name__ == "__main__":
    main()
