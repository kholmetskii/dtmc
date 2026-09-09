#!/usr/bin/env python3
"""Generate byte-identical benchmark inputs for dtmc and PyDTMC.

Each payload is little-endian float64 data in row-major order.  The n*n
transition probabilities are followed by the n coordinates of the initial
distribution.  Metadata and SHA-256 digests are stored in manifest.json.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Callable

import numpy as np


DEFAULT_SIZES = (10, 25, 50, 100, 250, 500, 1000)
DEFAULT_SEEDS = (1729, 2718, 31415)
FAMILY_OFFSETS = {
    "dense": 10_000_000,
    "low-outdegree": 20_000_000,
    "absorbing": 30_000_000,
    "reducible": 40_000_000,
}


def normalized_integer_weights(rng: np.random.Generator, size: int) -> np.ndarray:
    values = rng.integers(1, 1001, size=size, dtype=np.int64).astype(np.float64)
    return values / np.sum(values, dtype=np.float64)


def make_dense(n: int, rng: np.random.Generator) -> tuple[np.ndarray, list[int], list[int]]:
    matrix = np.vstack([normalized_integer_weights(rng, n) for _ in range(n)])
    width = max(1, n // 20)
    targets = list(range(n // 2, min(n, n // 2 + width)))
    return matrix, targets, []


def make_low_outdegree(
    n: int, rng: np.random.Generator
) -> tuple[np.ndarray, list[int], list[int]]:
    matrix = np.zeros((n, n), dtype=np.float64)
    degree = min(4, n)
    for source in range(n):
        destinations = {source, (source + 1) % n}
        while len(destinations) < degree:
            destinations.add(int(rng.integers(0, n)))
        ordered = sorted(destinations)
        matrix[source, ordered] = normalized_integer_weights(rng, len(ordered))
    width = max(1, n // 20)
    targets = list(range(n // 2, min(n, n // 2 + width)))
    return matrix, targets, []


def make_absorbing(
    n: int, rng: np.random.Generator
) -> tuple[np.ndarray, list[int], list[int]]:
    if n < 4:
        raise ValueError("absorbing fixtures require at least four states")
    absorbing = [n - 2, n - 1]
    matrix = np.zeros((n, n), dtype=np.float64)
    for source in range(n - 2):
        # Every transient row has a positive edge to each absorber, so the
        # transient block has spectral radius strictly below one.
        matrix[source, :] = normalized_integer_weights(rng, n)
    matrix[n - 2, n - 2] = 1.0
    matrix[n - 1, n - 1] = 1.0
    return matrix, absorbing.copy(), absorbing


def _closed_ring(
    matrix: np.ndarray,
    members: list[int],
    rng: np.random.Generator,
) -> None:
    for offset, source in enumerate(members):
        destinations = sorted({source, members[(offset + 1) % len(members)]})
        matrix[source, destinations] = normalized_integer_weights(rng, len(destinations))


def make_reducible(
    n: int, rng: np.random.Generator
) -> tuple[np.ndarray, list[int], list[int]]:
    if n < 6:
        raise ValueError("reducible fixtures require at least six states")
    class_size = max(2, n // 10)
    transient_count = n - 2 * class_size
    class_a = list(range(transient_count, transient_count + class_size))
    class_b = list(range(transient_count + class_size, n))
    matrix = np.zeros((n, n), dtype=np.float64)

    for source in range(transient_count):
        # A dense transient SCC avoids exponentially tiny reachability values
        # that fall on PyDTMC's internal np.isclose threshold. Threshold-stress
        # behavior is intentionally outside this ordinary performance suite.
        ordered = list(range(transient_count)) + [class_a[0], class_b[0]]
        matrix[source, ordered] = normalized_integer_weights(rng, len(ordered))

    _closed_ring(matrix, class_a, rng)
    _closed_ring(matrix, class_b, rng)
    return matrix, class_a, []


GENERATORS: dict[
    str, Callable[[int, np.random.Generator], tuple[np.ndarray, list[int], list[int]]]
] = {
    "dense": make_dense,
    "low-outdegree": make_low_outdegree,
    "absorbing": make_absorbing,
    "reducible": make_reducible,
}


def make_initial(n: int, rng: np.random.Generator) -> np.ndarray:
    return normalized_integer_weights(rng, n)


def dataset_id(family: str, size: int, seed: int) -> str:
    return f"{family}-n{size:04d}-s{seed:05d}"


def generate(output_dir: Path, sizes: list[int], seeds: list[int]) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, object]] = []

    for family, family_generator in GENERATORS.items():
        for size in sizes:
            for seed in seeds:
                stream_seed = FAMILY_OFFSETS[family] + size * 1009 + seed
                rng = np.random.Generator(np.random.PCG64(stream_seed))
                matrix, targets, absorbing = family_generator(size, rng)
                initial = make_initial(size, rng)

                if matrix.shape != (size, size) or initial.shape != (size,):
                    raise AssertionError("generator returned an invalid shape")
                if not np.all(np.isfinite(matrix)) or not np.all(matrix >= 0.0):
                    raise AssertionError("generator returned invalid probabilities")
                if not np.allclose(matrix.sum(axis=1), 1.0, atol=1e-14, rtol=0.0):
                    raise AssertionError("matrix rows are not stochastic")
                if not np.isclose(initial.sum(), 1.0, atol=1e-14, rtol=0.0):
                    raise AssertionError("initial distribution is not stochastic")

                identifier = dataset_id(family, size, seed)
                file_name = f"{identifier}.f64"
                path = output_dir / file_name
                payload = np.concatenate((matrix.ravel(order="C"), initial)).astype("<f8")
                payload.tofile(path)
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                records.append(
                    {
                        "id": identifier,
                        "family": family,
                        "size": size,
                        "seed": seed,
                        "stream_seed": stream_seed,
                        "file": file_name,
                        "sha256": digest,
                        "targets": targets,
                        "absorbing": absorbing,
                    }
                )

    return {
        "schema_version": 1,
        "format": {
            "dtype": "float64",
            "byte_order": "little",
            "layout": "row-major",
            "contents": ["transition_matrix[size*size]", "initial_distribution[size]"],
        },
        "generator": {
            "name": "benchmark/data/generate.py",
            "numpy_version": np.__version__,
            "bit_generator": "PCG64",
        },
        "datasets": records,
    }


def parse_csv_ints(value: str) -> list[int]:
    return [int(item) for item in value.split(",") if item]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "generated",
    )
    parser.add_argument("--sizes", default=",".join(map(str, DEFAULT_SIZES)))
    parser.add_argument("--seeds", default=",".join(map(str, DEFAULT_SEEDS)))
    args = parser.parse_args()

    sizes = parse_csv_ints(args.sizes)
    seeds = parse_csv_ints(args.seeds)
    if not sizes or any(size < 6 for size in sizes):
        parser.error("all sizes must be integers greater than or equal to six")
    if not seeds:
        parser.error("at least one seed is required")

    manifest = generate(args.output_dir, sizes, seeds)
    manifest_path = args.output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"generated {len(manifest['datasets'])} datasets in {args.output_dir}")
    print(f"manifest: {manifest_path}")


if __name__ == "__main__":
    main()
