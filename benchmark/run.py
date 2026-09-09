#!/usr/bin/env python3
"""Reproducible orchestration for the dtmc versus PyDTMC benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any


BENCHMARK_ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = BENCHMARK_ROOT.parent
DATA_DIR = BENCHMARK_ROOT / "data"
DATA_ROOT = DATA_DIR / "generated"
RESULTS_ROOT = BENCHMARK_ROOT / "results"
VENV = BENCHMARK_ROOT / ".venv"
CABAL_PROJECT = "benchmark/cabal.project"
SPEC = json.loads((DATA_DIR / "spec.json").read_text(encoding="utf-8"))
SIZES = tuple(map(int, SPEC["sizes"]))
SEEDS = tuple(map(int, SPEC["seeds"]))
FAMILIES = tuple(map(str, SPEC["families"]))
THREAD_ENV = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "MPLCONFIGDIR": str(BENCHMARK_ROOT / ".cache" / "matplotlib"),
    "XDG_CACHE_HOME": str(BENCHMARK_ROOT / ".cache"),
    "MPLBACKEND": "Agg",
}


def run(command: list[str], *, env: dict[str, str] | None = None) -> None:
    printable = " ".join(command)
    print(f"+ {printable}", flush=True)
    subprocess.run(command, cwd=PROJECT_ROOT, env=env, check=True)


def capture(command: list[str], *, env: dict[str, str] | None = None) -> str:
    try:
        return subprocess.run(
            command,
            cwd=PROJECT_ROOT,
            env=env,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as problem:
        return f"unavailable: {problem}"


def python_executable() -> Path:
    name = "python.exe" if os.name == "nt" else "python"
    return VENV / ("Scripts" if os.name == "nt" else "bin") / name


def haskell_executable() -> Path:
    result = subprocess.run(
        ["cabal", "list-bin", f"--project-file={CABAL_PROJECT}", "dtmc-bench"],
        cwd=PROJECT_ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    path = Path(result.stdout.strip())
    if not path.exists():
        raise RuntimeError("dtmc-bench is not built; run the bootstrap command first")
    return path


def bootstrap() -> None:
    (BENCHMARK_ROOT / ".cache" / "matplotlib").mkdir(parents=True, exist_ok=True)
    if not python_executable().exists():
        run([sys.executable, "-m", "venv", str(VENV)])
    run(
        [
            str(python_executable()),
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            "-r",
            str(BENCHMARK_ROOT / "python" / "requirements.lock"),
        ]
    )
    run(["cabal", "build", f"--project-file={CABAL_PROJECT}", "dtmc-bench"])


def generate(mode: str) -> None:
    sizes = SIZES if mode == "full" else (10,)
    seeds = SEEDS if mode == "full" else (SEEDS[0],)
    python = str(python_executable()) if python_executable().exists() else sys.executable
    run(
        [
            python,
            str(DATA_DIR / "generate.py"),
            "--sizes",
            ",".join(map(str, sizes)),
            "--seeds",
            ",".join(map(str, seeds)),
        ]
    )


def benchmark_environment(entry: dict[str, Any]) -> dict[str, str]:
    environment = os.environ.copy()
    environment.update(THREAD_ENV)
    environment.update(
        {
            "DTMC_BENCH_DATA": str(DATA_ROOT),
            "DTMC_BENCH_FAMILY": str(entry["family"]),
            "DTMC_BENCH_SIZE": str(entry["size"]),
            "DTMC_BENCH_SEED": str(entry["seed"]),
        }
    )
    return environment


def manifest_entries(mode: str) -> list[dict[str, Any]]:
    manifest = json.loads((DATA_ROOT / "manifest.json").read_text(encoding="utf-8"))
    entries = list(manifest["datasets"])
    if mode == "smoke":
        entries = [
            entry
            for entry in entries
            if entry["size"] == 10 and entry["seed"] == SEEDS[0]
        ]
    return entries


def verify(mode: str) -> None:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    (BENCHMARK_ROOT / ".cache" / "matplotlib").mkdir(parents=True, exist_ok=True)
    haskell_output = RESULTS_ROOT / "haskell-verification.json"
    verification_output = RESULTS_ROOT / "verification.json"
    environment = os.environ.copy()
    environment.update(THREAD_ENV)
    environment["DTMC_BENCH_DATA"] = str(DATA_ROOT)
    environment["DTMC_VERIFY_MAX_SIZE"] = "100" if mode == "full" else "10"
    run(
        [
            str(haskell_executable()),
            "--verify-json",
            str(haskell_output),
        ],
        env=environment,
    )
    run(
        [
            str(python_executable()),
            str(BENCHMARK_ROOT / "python" / "verify.py"),
            str(haskell_output),
            "--output",
            str(verification_output),
        ],
        env=environment,
    )


def benchmark(mode: str) -> None:
    raw_haskell = RESULTS_ROOT / "raw" / "haskell"
    raw_python = RESULTS_ROOT / "raw" / "python"
    raw_haskell.mkdir(parents=True, exist_ok=True)
    raw_python.mkdir(parents=True, exist_ok=True)

    for entry in manifest_entries(mode):
        name = str(entry["id"])
        environment = benchmark_environment(entry)
        haskell_json = raw_haskell / f"{name}.json"
        haskell_csv = raw_haskell / f"{name}.csv"
        haskell_raw = raw_haskell / f"{name}.raw"
        python_json = raw_python / f"{name}.json"
        for output in (haskell_json, haskell_csv, haskell_raw, python_json):
            output.unlink(missing_ok=True)
        criterion_command = [
            str(haskell_executable()),
            "--json",
            str(haskell_json),
            "--csv",
            str(haskell_csv),
            "--raw",
            str(haskell_raw),
            "--time-limit",
            "0.1" if mode == "smoke" else "5",
            "--resamples",
            "100" if mode == "smoke" else "10000",
        ]
        run(criterion_command, env=environment)

        pyperf_command = [
            str(python_executable()),
            str(BENCHMARK_ROOT / "python" / "bench.py"),
            "--copy-env",
            "-o",
            str(python_json),
        ]
        if mode == "smoke":
            pyperf_command.append("--fast")
        run(pyperf_command, env=environment)


def record_environment() -> None:
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update(THREAD_ENV)
    python = str(python_executable()) if python_executable().exists() else sys.executable
    haskell_binary = haskell_executable()

    def digest(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    record = {
        "captured_at_utc": capture(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_brand": capture(["sysctl", "-n", "machdep.cpu.brand_string"]),
        "python": capture([python, "--version"]),
        "python_packages": capture([python, "-m", "pip", "freeze"]),
        "numpy_configuration": capture(
            [python, "-c", "import numpy as n; n.show_config()"], env=environment
        ),
        "ghc": capture(["ghc", "--version"]),
        "cabal": capture(["cabal", "--version"]),
        "haskell_benchmark_binary": str(haskell_binary),
        "haskell_benchmark_linkage": capture(["otool", "-L", str(haskell_binary)]),
        "haskell_options": ["-O2", "-threaded", "-rtsopts"],
        "cabal_freeze_sha256": digest(BENCHMARK_ROOT / "cabal.project.freeze"),
        "python_lock_sha256": digest(BENCHMARK_ROOT / "python" / "requirements.lock"),
        "dataset_spec_sha256": digest(DATA_DIR / "spec.json"),
        "git_commit": capture(["git", "rev-parse", "HEAD"]),
        "git_status": capture(["git", "status", "--short"]),
        "thread_environment": THREAD_ENV,
    }
    (RESULTS_ROOT / "environment.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8"
    )


def profile() -> None:
    entries = manifest_entries("full")
    candidates = [entry for entry in entries if entry["family"] == "dense"]
    entry = next((item for item in candidates if item["size"] == 100), candidates[0])
    environment = benchmark_environment(entry)
    profile_root = BENCHMARK_ROOT / "profiles"
    profile_root.mkdir(parents=True, exist_ok=True)
    for operation in ("stationary", "hitting-probability", "simulation/100000"):
        safe = operation.replace("/", "-")
        haskell_path = profile_root / f"haskell-{safe}.txt"
        with haskell_path.open("w", encoding="utf-8") as stream:
            subprocess.run(
                [
                    str(haskell_executable()),
                    "--time-limit",
                    "0.2",
                    "--resamples",
                    "100",
                    f"{entry['id']}/{operation}",
                    "+RTS",
                    "-s",
                    "-RTS",
                ],
                cwd=PROJECT_ROOT,
                env=environment,
                check=True,
                stdout=stream,
                stderr=subprocess.STDOUT,
            )
        python_environment = environment.copy()
        python_environment["DTMC_BENCH_OPERATION"] = operation
        run(
            [
                str(python_executable()),
                str(BENCHMARK_ROOT / "python" / "bench.py"),
                "--copy-env",
                "--debug-single-value",
                "--profile",
                str(profile_root / f"python-{safe}.prof"),
            ],
            env=python_environment,
        )


def analyse() -> None:
    run([str(python_executable()), str(BENCHMARK_ROOT / "analyse.py")])
    run([str(python_executable()), str(BENCHMARK_ROOT / "plot.py")])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        choices=(
            "bootstrap",
            "generate",
            "verify",
            "benchmark",
            "profile",
            "analyse",
            "all",
        ),
    )
    parser.add_argument("--mode", choices=("smoke", "full"), default="smoke")
    args = parser.parse_args()

    if args.command in {"bootstrap", "all"}:
        bootstrap()
    if args.command in {"generate", "all"}:
        generate(args.mode)
    if args.command in {"verify", "all"}:
        verify(args.mode)
    if args.command in {"benchmark", "all"}:
        benchmark(args.mode)
    if args.command == "profile":
        profile()
    if args.command in {"analyse", "all"}:
        record_environment()
        analyse()


if __name__ == "__main__":
    main()
