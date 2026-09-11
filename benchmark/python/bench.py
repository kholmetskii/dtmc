#!/usr/bin/env python3
"""Pyperf benchmarks matching the public dtmc benchmark operations."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any, Callable

import numpy as np
import pyperf
from pydtmc import MarkovChain


@dataclass(frozen=True)
class Dataset:
    metadata: dict[str, Any]
    matrix: np.ndarray
    initial: np.ndarray

    @property
    def name(self) -> str:
        return str(self.metadata["id"])

    @property
    def family(self) -> str:
        return str(self.metadata["family"])

    @property
    def size(self) -> int:
        return int(self.metadata["size"])

    @property
    def seed(self) -> int:
        return int(self.metadata["seed"])

    @property
    def targets(self) -> list[int]:
        return [int(value) for value in self.metadata["targets"]]

    @property
    def absorbing(self) -> list[int]:
        return [int(value) for value in self.metadata["absorbing"]]

    @property
    def competing(self) -> list[int]:
        targets = set(self.targets)
        return [state for state in range(self.size) if state not in targets][
            : max(1, len(targets))
        ]


def data_root() -> Path:
    return Path(os.environ.get("DTMC_BENCH_DATA", "benchmark/data/generated")).resolve()


def load_selected_dataset() -> Dataset:
    root = data_root()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    family = os.environ.get("DTMC_BENCH_FAMILY")
    size = os.environ.get("DTMC_BENCH_SIZE")
    seed = os.environ.get("DTMC_BENCH_SEED")
    matches = [
        entry
        for entry in manifest["datasets"]
        if (family is None or entry["family"] == family)
        and (size is None or entry["size"] == int(size))
        and (seed is None or entry["seed"] == int(seed))
    ]
    if len(matches) != 1:
        raise RuntimeError(
            "select exactly one dataset with DTMC_BENCH_FAMILY, "
            "DTMC_BENCH_SIZE, and DTMC_BENCH_SEED"
        )
    entry = matches[0]
    payload_path = root / entry["file"]
    digest = hashlib.sha256(payload_path.read_bytes()).hexdigest()
    if digest != entry["sha256"]:
        raise RuntimeError(f"SHA-256 mismatch for {payload_path}")
    n = int(entry["size"])
    values = np.fromfile(payload_path, dtype="<f8")
    expected = n * n + n
    if values.size != expected:
        raise RuntimeError(f"{payload_path} contains {values.size} values, expected {expected}")
    matrix = values[: n * n].reshape((n, n)).copy()
    initial = values[n * n :].copy()
    matrix.flags.writeable = False
    initial.flags.writeable = False
    return Dataset(entry, matrix, initial)


def checksum_array(value: Any) -> float:
    array = np.asarray(value, dtype=np.float64).ravel(order="C")
    return float(np.sum(array, dtype=np.float64))


def consume_chain(chain: MarkovChain) -> float:
    return checksum_array(chain.p)


def consume_classes(classes: list[list[str]]) -> int:
    return sum(int(state) - 1 for group in classes for state in group)


def consume_classes_access(classes: list[list[str]]) -> bool:
    return classes is not None


def consume_bool(value: bool) -> int:
    return int(value)


def consume_float(value: float) -> float:
    return float(value)


def consume_sequence(value: list[int]) -> int:
    return sum(map(int, value))


def consume_optional_array(value: Any) -> float:
    if value is None:
        raise RuntimeError("PyDTMC operation unexpectedly returned None")
    return checksum_array(value)


def occupation_including_initial(chain: MarkovChain) -> np.ndarray:
    """Adapt PyDTMC's positive-time visit counts to dtmc's time-zero convention."""
    values = chain.mean_number_visits()
    if values is None:
        raise RuntimeError("PyDTMC operation unexpectedly returned None")
    adjusted = np.asarray(values, dtype=np.float64).copy()
    adjusted.flat[:: adjusted.shape[0] + 1] += 1.0
    return adjusted


def bounded_return_probability(chain: MarkovChain, steps: int) -> float:
    probabilities = chain.first_passage_probabilities(steps, 0, [0])
    return float(np.sum(probabilities, dtype=np.float64))


def finite_time_probability(
    chain: MarkovChain,
    steps: int,
    initial: np.ndarray,
    destination: int,
) -> float:
    distribution = chain.redistribute(
        steps,
        initial_status=initial,
        output_last=True,
    )
    return float(distribution[destination])


def bounded_hitting_masses(
    chain: MarkovChain,
    steps: int,
    initial_state: int,
    targets: list[int],
) -> np.ndarray:
    return np.asarray(
        chain.first_passage_probabilities(steps, initial_state, targets),
        dtype=np.float64,
    )


def target_absorbing_matrix(
    matrix: np.ndarray,
    targets: list[int],
) -> np.ndarray:
    """Preserve first hits while making PyDTMC's target events disjoint."""
    adjusted = matrix.copy()
    adjusted[targets, :] = 0.0
    adjusted[targets, targets] = 1.0
    adjusted.flags.writeable = False
    return adjusted


def exact_hitting_probability(
    chain: MarkovChain,
    steps: int,
    initial_state: int,
    targets: list[int],
) -> float:
    masses = bounded_hitting_masses(chain, steps, initial_state, targets)
    return float(masses[-1])


def bounded_hitting_probability(
    chain: MarkovChain,
    steps: int,
    initial_state: int,
    targets: list[int],
) -> float:
    masses = bounded_hitting_masses(chain, steps, initial_state, targets)
    return float(np.sum(masses, dtype=np.float64))


def upper_hitting_probability(
    chain: MarkovChain,
    steps: int,
    initial_state: int,
    targets: list[int],
) -> float:
    return 1.0 - bounded_hitting_probability(
        chain,
        steps,
        initial_state,
        targets,
    )


def bounded_visit_expectation(
    chain: MarkovChain, steps: int, rewards: np.ndarray
) -> float:
    return float(chain.expected_rewards(steps - 1, rewards)[0])


def fresh_chain(dataset: Dataset) -> MarkovChain:
    return MarkovChain(dataset.matrix)


def operation_selected(name: str) -> bool:
    selected_many = os.environ.get("DTMC_BENCH_OPERATIONS")
    if selected_many is not None:
        operations = [item for item in selected_many.split(",") if item]
        return any(name.endswith(f"/{operation}") for operation in operations)

    selected = os.environ.get("DTMC_BENCH_OPERATION")
    return selected is None or selected in name


def register_fresh(
    runner: pyperf.Runner,
    name: str,
    setup: Callable[[], Any],
    operation: Callable[[Any], Any],
    consume: Callable[[Any], Any],
) -> None:
    """Exclude fixture setup while still using a fresh cache for every sample."""
    if not operation_selected(name):
        return

    def timed(loops: int) -> float:
        elapsed = 0.0
        for _ in range(loops):
            fixture = setup()
            started = perf_counter()
            result = operation(fixture)
            consume(result)
            elapsed += perf_counter() - started
        return elapsed

    runner.bench_time_func(name, timed)


def register_warm(
    runner: pyperf.Runner,
    name: str,
    fixture: Any,
    operation: Callable[[Any], Any],
    consume: Callable[[Any], Any],
) -> None:
    if not operation_selected(name):
        return
    inner_loops = 10_000

    def timed(loops: int) -> float:
        started = perf_counter()
        for _ in range(loops * inner_loops):
            consume(operation(fixture))
        return perf_counter() - started

    runner.bench_time_func(name, timed, inner_loops=inner_loops)


def register_benchmarks(runner: pyperf.Runner, dataset: Dataset) -> None:
    prefix = dataset.name
    matrix = dataset.matrix
    initial = dataset.initial
    targets = dataset.targets
    competing = dataset.competing
    target = targets[0]
    hitting_matrix = target_absorbing_matrix(matrix, targets)
    point_initial = np.zeros(dataset.size, dtype=np.float64)
    point_initial[0] = 1.0
    point_initial.flags.writeable = False

    register_fresh(
        runner,
        f"{prefix}/construction/public-consumed",
        lambda: matrix,
        MarkovChain,
        consume_chain,
    )

    for steps in (1, 10, 100):
        label = "one-step" if steps == 1 else f"k-{steps}"
        register_fresh(
            runner,
            f"{prefix}/evolution/{label}",
            lambda: fresh_chain(dataset),
            lambda chain, count=steps: chain.redistribute(
                count, initial_status=initial, output_last=True
            ),
            checksum_array,
        )

    if dataset.size <= 100:
        register_fresh(
            runner,
            f"{prefix}/finite-time/step",
            lambda: fresh_chain(dataset),
            lambda chain: chain.p[0, target],
            consume_float,
        )
        for steps in (10, 100):
            register_fresh(
                runner,
                f"{prefix}/finite-time/n-step/k-{steps}",
                lambda: fresh_chain(dataset),
                lambda chain, count=steps: finite_time_probability(
                    chain,
                    count,
                    point_initial,
                    target,
                ),
                consume_float,
            )
            register_fresh(
                runner,
                f"{prefix}/finite-time/observation/k-{steps}",
                lambda: fresh_chain(dataset),
                lambda chain, count=steps: finite_time_probability(
                    chain,
                    count,
                    initial,
                    target,
                ),
                consume_float,
            )

    if dataset.size <= 500:
        for exponent in (2, 10, 100):
            register_fresh(
                runner,
                f"{prefix}/power/{exponent}",
                lambda: fresh_chain(dataset),
                lambda chain, value=exponent: chain.to_nth_order(value),
                consume_chain,
            )

    if dataset.family != "dense" or dataset.size <= 500:
        register_fresh(
            runner,
            f"{prefix}/structure/classes-lifecycle",
            lambda: matrix,
            lambda value: MarkovChain(value).communicating_classes,
            consume_classes,
        )
        register_fresh(
            runner,
            f"{prefix}/structure/classes-cold",
            lambda: fresh_chain(dataset),
            lambda chain: chain.communicating_classes,
            consume_classes,
        )

        warm_classes = fresh_chain(dataset)
        consume_classes(warm_classes.communicating_classes)
        register_warm(
            runner,
            f"{prefix}/structure/classes-warm/access",
            warm_classes,
            lambda chain: chain.communicating_classes,
            consume_classes_access,
        )
        register_warm(
            runner,
            f"{prefix}/structure/classes-warm/consumed",
            warm_classes,
            lambda chain: chain.communicating_classes,
            consume_classes,
        )
        register_fresh(
            runner,
            f"{prefix}/structure/irreducible-cold",
            lambda: fresh_chain(dataset),
            lambda chain: chain.is_irreducible,
            consume_bool,
        )

        if dataset.family == "periodic":
            register_fresh(
                runner,
                f"{prefix}/structure/period-cold",
                lambda: fresh_chain(dataset),
                lambda chain: chain.period,
                consume_float,
            )
            warm_period = fresh_chain(dataset)
            consume_float(warm_period.period)
            register_warm(
                runner,
                f"{prefix}/structure/period-warm",
                warm_period,
                lambda chain: chain.period,
                consume_float,
            )
            register_fresh(
                runner,
                f"{prefix}/structure/cyclic-classes-cold",
                lambda: fresh_chain(dataset),
                lambda chain: chain.cyclic_classes,
                consume_classes,
            )
            warm_cyclic_classes = fresh_chain(dataset)
            consume_classes(warm_cyclic_classes.cyclic_classes)
            register_warm(
                runner,
                f"{prefix}/structure/cyclic-classes-warm",
                warm_cyclic_classes,
                lambda chain: chain.cyclic_classes,
                consume_classes,
            )

        warm_irreducible = fresh_chain(dataset)
        consume_bool(warm_irreducible.is_irreducible)
        register_warm(
            runner,
            f"{prefix}/structure/irreducible-warm",
            warm_irreducible,
            lambda chain: chain.is_irreducible,
            consume_bool,
        )

    if dataset.size <= 500:
        register_fresh(
            runner,
            f"{prefix}/stationary",
            lambda: fresh_chain(dataset),
            lambda chain: chain.pi,
            checksum_array,
        )
        register_fresh(
            runner,
            f"{prefix}/hitting-probability/cold-all-states",
            lambda: fresh_chain(dataset),
            lambda chain: chain.hitting_probabilities(targets),
            checksum_array,
        )
        register_fresh(
            runner,
            f"{prefix}/hitting-time/cold-all-states",
            lambda: fresh_chain(dataset),
            lambda chain: chain.hitting_times(targets),
            checksum_array,
        )

        if dataset.family in {"dense", "low-outdegree"}:
            register_fresh(
                runner,
                f"{prefix}/race/forward-committor",
                lambda: fresh_chain(dataset),
                lambda chain: chain.committor_probabilities(
                    "forward", competing, targets
                ),
                consume_optional_array,
            )

    if dataset.family in {"dense", "low-outdegree"} and dataset.size <= 500:
        register_fresh(
            runner,
            f"{prefix}/return/mean-recurrence",
            lambda: fresh_chain(dataset),
            lambda chain: chain.mean_recurrence_times(),
            consume_optional_array,
        )

    if dataset.size <= 100:
        rewards = np.zeros(dataset.size, dtype=np.float64)
        rewards[targets[0]] = 1.0
        rewards.flags.writeable = False
        for steps in (10, 100):
            for label, operation in (
                ("exact", exact_hitting_probability),
                ("at-most", bounded_hitting_probability),
                ("greater-than", upper_hitting_probability),
            ):
                register_fresh(
                    runner,
                    f"{prefix}/hitting/bounded/{label}/k-{steps}",
                    lambda: MarkovChain(hitting_matrix),
                    lambda chain, count=steps, query=operation: query(
                        chain,
                        count,
                        0,
                        targets,
                    ),
                    consume_float,
                )
            register_fresh(
                runner,
                f"{prefix}/return/bounded/k-{steps}",
                lambda: fresh_chain(dataset),
                lambda chain, count=steps: bounded_return_probability(chain, count),
                consume_float,
            )
            register_fresh(
                runner,
                f"{prefix}/visits/bounded-expectation/k-{steps}",
                lambda: fresh_chain(dataset),
                lambda chain, count=steps: bounded_visit_expectation(
                    chain, count, rewards
                ),
                consume_float,
            )

    if dataset.family == "absorbing" and dataset.size <= 500:
        register_fresh(
            runner,
            f"{prefix}/fundamental-matrix",
            lambda: fresh_chain(dataset),
            lambda chain: chain.fundamental_matrix,
            consume_optional_array,
        )
        register_fresh(
            runner,
            f"{prefix}/absorption/probabilities",
            lambda: fresh_chain(dataset),
            lambda chain: chain.absorption_probabilities(),
            consume_optional_array,
        )
        register_fresh(
            runner,
            f"{prefix}/absorption-time",
            lambda: fresh_chain(dataset),
            lambda chain: chain.mean_absorption_times(),
            consume_optional_array,
        )

    if dataset.family in {"absorbing", "reducible"} and dataset.size <= 250:
        register_fresh(
            runner,
            f"{prefix}/occupation-matrix",
            lambda: fresh_chain(dataset),
            occupation_including_initial,
            consume_optional_array,
        )

    if dataset.family in {"dense", "low-outdegree"}:
        for transitions in (10_000, 100_000):
            register_fresh(
                runner,
                f"{prefix}/simulation/{transitions}",
                lambda: fresh_chain(dataset),
                lambda chain, count=transitions: chain.simulate(
                    count,
                    initial_state=0,
                    output_indices=True,
                    seed=dataset.seed,
                ),
                consume_sequence,
            )


def main() -> None:
    dataset = load_selected_dataset()
    # One cold operation per value prevents pyperf calibration from recreating
    # an expensive untimed fixture thousands of times. Three processes with ten
    # values give 30 independent full-operation samples in the normal suite.
    runner = pyperf.Runner(processes=3, values=10, warmups=3, loops=1)
    runner.metadata["dataset"] = dataset.name
    runner.metadata["dataset_sha256"] = dataset.metadata["sha256"]
    runner.metadata["pydtmc_version"] = "9.0.0"
    register_benchmarks(runner, dataset)


if __name__ == "__main__":
    main()
