#!/usr/bin/env python3
"""Compare untimed dtmc outputs with PyDTMC on the shared fixtures."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

import numpy as np
from pydtmc import MarkovChain

from bench import (
    Dataset,
    bounded_return_probability,
    bounded_visit_expectation,
    data_root,
    occupation_including_initial,
)


@dataclass
class Comparison:
    dataset: str
    operation: str
    max_absolute_error: float
    max_relative_error: float
    passed: bool


def load_dataset(root: Path, entry: dict[str, Any]) -> Dataset:
    n = int(entry["size"])
    values = np.fromfile(root / entry["file"], dtype="<f8")
    if values.size != n * n + n:
        raise RuntimeError(f"invalid payload size for {entry['id']}")
    return Dataset(
        entry,
        values[: n * n].reshape((n, n)).copy(),
        values[n * n :].copy(),
    )


def expectation_array(value: Any) -> np.ndarray:
    def decode(item: Any) -> float:
        if item is None:
            return math.inf
        if isinstance(item, dict) and "finite" in item:
            return float(item["finite"])
        if isinstance(item, list):
            raise TypeError("nested expectation list reached scalar decoder")
        return float(item)

    if isinstance(value, list) and value and isinstance(value[0], list):
        return np.asarray([[decode(item) for item in row] for row in value], dtype=float)
    return np.asarray([decode(item) for item in value], dtype=float)


def compare_numeric(
    dataset: str,
    operation: str,
    haskell: Any,
    python: Any,
    atol: float,
    rtol: float,
) -> Comparison:
    left = np.asarray(haskell, dtype=float)
    right = np.asarray(python, dtype=float)
    if left.shape != right.shape:
        raise AssertionError(
            f"{dataset}/{operation}: shape mismatch {left.shape} != {right.shape}"
        )
    if not np.array_equal(np.isinf(left), np.isinf(right)):
        raise AssertionError(f"{dataset}/{operation}: infinity masks differ")
    finite = np.isfinite(left) & np.isfinite(right)
    if not np.array_equal(np.isnan(left), np.isnan(right)):
        raise AssertionError(f"{dataset}/{operation}: NaN masks differ")
    if np.any(finite):
        absolute = np.abs(left[finite] - right[finite])
        scale = np.maximum(np.maximum(np.abs(left[finite]), np.abs(right[finite])), 1e-300)
        relative = absolute / scale
        max_absolute = float(np.max(absolute))
        max_relative = float(np.max(relative))
        passed = bool(np.all(absolute <= atol + rtol * np.abs(right[finite])))
    else:
        max_absolute = 0.0
        max_relative = 0.0
        passed = True
    return Comparison(dataset, operation, max_absolute, max_relative, passed)


def state_index(label: str) -> int:
    # Default PyDTMC labels are one-based decimal strings.
    return int(label) - 1


def normalized_partition(groups: list[list[int]]) -> list[list[int]]:
    return sorted(sorted(group) for group in groups)


def compare_dataset(
    record: dict[str, Any],
    dataset: Dataset,
    atol: float,
    rtol: float,
) -> list[Comparison]:
    name = dataset.name
    chain = MarkovChain(dataset.matrix)
    targets = dataset.targets
    rewards = np.zeros(dataset.size, dtype=np.float64)
    rewards[targets[0]] = 1.0
    comparisons = [
        compare_numeric(
            name,
            "evolve_1",
            record["evolve_1"],
            chain.redistribute(1, initial_status=dataset.initial, output_last=True),
            atol,
            rtol,
        ),
        compare_numeric(
            name,
            "evolve_10",
            record["evolve_10"],
            chain.redistribute(10, initial_status=dataset.initial, output_last=True),
            atol,
            rtol,
        ),
        compare_numeric(
            name,
            "power_10",
            record["power_10"],
            chain.to_nth_order(10).p,
            atol,
            rtol,
        ),
        compare_numeric(
            name,
            "hitting_probability",
            record["hitting_probability"],
            chain.hitting_probabilities(targets),
            atol,
            rtol,
        ),
        compare_numeric(
            name,
            "hitting_time",
            expectation_array(record["hitting_time"]),
            chain.hitting_times(targets),
            atol,
            rtol,
        ),
        compare_numeric(
            name,
            "return_bounded_10",
            [record["return_bounded_10"]],
            [bounded_return_probability(chain, 10)],
            atol,
            rtol,
        ),
        compare_numeric(
            name,
            "return_bounded_100",
            [record["return_bounded_100"]],
            [bounded_return_probability(chain, 100)],
            atol,
            rtol,
        ),
        compare_numeric(
            name,
            "visits_bounded_expectation_10",
            [record["visits_bounded_expectation_10"]],
            [bounded_visit_expectation(chain, 10, rewards)],
            atol,
            rtol,
        ),
        compare_numeric(
            name,
            "visits_bounded_expectation_100",
            [record["visits_bounded_expectation_100"]],
            [bounded_visit_expectation(chain, 100, rewards)],
            atol,
            rtol,
        ),
    ]

    if dataset.family in {"dense", "low-outdegree"}:
        py_race = chain.committor_probabilities(
            "forward", dataset.competing, targets
        )
        if py_race is None:
            raise AssertionError(f"{name}: PyDTMC returned no committor probabilities")
        comparisons.append(
            compare_numeric(
                name,
                "race_probability",
                record["race_probability"],
                py_race,
                atol,
                rtol,
            )
        )
        py_return_mean = chain.mean_recurrence_times()
        if py_return_mean is None:
            raise AssertionError(f"{name}: PyDTMC returned no mean recurrence times")
        comparisons.append(
            compare_numeric(
                name,
                "return_mean",
                expectation_array(record["return_mean"]),
                py_return_mean,
                atol,
                rtol,
            )
        )

    haskell_classes = sorted(sorted(group) for group in record["classes"])
    python_classes = sorted(sorted(state_index(state) for state in group) for group in chain.communicating_classes)
    if haskell_classes != python_classes:
        raise AssertionError(f"{name}: communicating classes differ")
    if bool(record["irreducible"]) != bool(chain.is_irreducible):
        raise AssertionError(f"{name}: irreducibility differs")

    if dataset.family == "periodic":
        if int(record["chain_period"]) != int(chain.period):
            raise AssertionError(f"{name}: chain periods differ")
        haskell_cyclic = normalized_partition(record["cyclic_classes"])
        python_cyclic = normalized_partition(
            [[state_index(state) for state in group] for group in chain.cyclic_classes]
        )
        if haskell_cyclic != python_cyclic:
            raise AssertionError(f"{name}: cyclic classes differ")

    haskell_stationary = {
        tuple(item["members"]): np.asarray(item["weights"], dtype=float)
        for item in record["stationary"]
    }
    python_stationary: dict[tuple[int, ...], np.ndarray] = {}
    for vector in chain.pi:
        array = np.asarray(vector, dtype=float)
        support = tuple(int(i) for i in np.flatnonzero(array > 0.0))
        python_stationary[support] = array
    if set(haskell_stationary) != set(python_stationary):
        raise AssertionError(f"{name}: stationary supports differ")
    for support, haskell_vector in haskell_stationary.items():
        comparisons.append(
            compare_numeric(
                name,
                f"stationary[{','.join(map(str, support))}]",
                haskell_vector,
                python_stationary[support],
                atol,
                rtol,
            )
        )

    if dataset.family == "absorbing":
        fundamental = record["fundamental"]
        py_fundamental = chain.fundamental_matrix
        if py_fundamental is None:
            raise AssertionError(f"{name}: PyDTMC did not classify fixture as absorbing")
        comparisons.append(
            compare_numeric(
                name,
                "fundamental",
                fundamental["values"],
                py_fundamental,
                atol,
                rtol,
            )
        )
        absorption_probability = record["absorption_probability"]
        py_absorbing = [state_index(state) for state in chain.absorbing_states]
        py_transient = [state_index(state) for state in chain.transient_states]
        if absorption_probability["absorbing"] != py_absorbing:
            raise AssertionError(f"{name}: absorbing-state order differs")
        if absorption_probability["transient"] != py_transient:
            raise AssertionError(f"{name}: transient-state order differs")
        py_absorption_probability = chain.absorption_probabilities()
        if py_absorption_probability is None:
            raise AssertionError(f"{name}: PyDTMC returned no absorption probabilities")
        comparisons.append(
            compare_numeric(
                name,
                "absorption_probability",
                absorption_probability["values"],
                py_absorption_probability,
                atol,
                rtol,
            )
        )
        transient = [int(index) for index in fundamental["states"]]
        py_absorption = chain.mean_absorption_times()
        if py_absorption is None:
            raise AssertionError(f"{name}: PyDTMC returned no absorption times")
        hs_absorption = expectation_array(record["absorption_time"])[transient]
        comparisons.append(
            compare_numeric(
                name,
                "absorption_time",
                hs_absorption,
                py_absorption,
                atol,
                rtol,
            )
        )

    if dataset.family in {"absorbing", "reducible"}:
        py_occupation = occupation_including_initial(chain)
        comparisons.append(
            compare_numeric(
                name,
                "occupation",
                expectation_array(record["occupation"]),
                py_occupation,
                atol,
                rtol,
            )
        )

    return comparisons


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("haskell_json", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--atol", type=float, default=1e-10)
    parser.add_argument("--rtol", type=float, default=1e-8)
    args = parser.parse_args()

    root = data_root()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    entries = {entry["id"]: entry for entry in manifest["datasets"]}
    haskell = json.loads(args.haskell_json.read_text(encoding="utf-8"))
    comparisons: list[Comparison] = []
    skipped: list[str] = []
    for record in haskell["datasets"]:
        if record.get("skipped"):
            skipped.append(record["id"])
            continue
        entry = entries[record["id"]]
        comparisons.extend(
            compare_dataset(record, load_dataset(root, entry), args.atol, args.rtol)
        )

    failed = [comparison for comparison in comparisons if not comparison.passed]
    output = {
        "schema_version": 1,
        "absolute_tolerance": args.atol,
        "relative_tolerance": args.rtol,
        "comparisons": [asdict(comparison) for comparison in comparisons],
        "skipped": skipped,
        "passed": not failed,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    if failed:
        names = ", ".join(f"{item.dataset}/{item.operation}" for item in failed)
        raise SystemExit(f"numerical verification failed: {names}")
    print(f"verified {len(comparisons)} operation results; maximum size {haskell['maximum_size']}")


if __name__ == "__main__":
    main()
