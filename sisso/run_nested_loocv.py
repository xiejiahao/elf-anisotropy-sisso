#!/usr/bin/env python3
"""Run nested leave-one-material-out SISSO++ searches for D=1 and D=2."""

from __future__ import annotations

import argparse
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

try:
    from .workflow import dump_json, load_json, read_csv, read_prediction, regression_metrics, write_csv
except ImportError:  # Direct execution: python sisso/run_nested_loocv.py
    from workflow import (  # type: ignore
        dump_json,
        load_json,
        read_csv,
        read_prediction,
        regression_metrics,
        write_csv,
    )


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = SCRIPT_DIR / "configs" / "search.json"
DEFAULT_INPUTS = SCRIPT_DIR / "work" / "inputs"
DEFAULT_WORK = SCRIPT_DIR / "work"


def parse_dimensions(value: str) -> list[int]:
    dimensions = [int(item.strip()) for item in value.split(",") if item.strip()]
    if not dimensions or any(item not in {1, 2} for item in dimensions):
        raise argparse.ArgumentTypeError("dimensions must be 1, 2, or 1,2")
    return dimensions


def make_sissopp_config(
    data_file: Path,
    dimension: int,
    n_sis: int,
    held_out: list[int],
    settings: dict,
) -> dict:
    if len(set(held_out)) != len(held_out):
        raise ValueError("held-out indices must be unique")
    return {
        "data_file": str(data_file.resolve()),
        "data_file_relative_to_json": False,
        "property_key": settings["property_key"],
        "opset": settings["opset"],
        "param_opset": [],
        "calc_type": "regression",
        "desc_dim": dimension,
        "n_sis_select": n_sis,
        "max_rung": settings["max_rung"],
        "max_leaves": 1,
        "n_residual": 1,
        "n_models_store": 1,
        "n_rung_store": 0,
        "n_rung_generate": 0,
        "min_abs_feat_val": 1e-12,
        "max_abs_feat_val": 1e12,
        "leave_out_inds": held_out,
        "leave_out_frac": 0.0,
        "fix_intercept": False,
        "max_feat_cross_correlation": 1.0,
        "nlopt_seed": settings["random_seed"],
        "global_param_opt": False,
        "reparam_residual": True,
    }


def run_one(
    binary: Path,
    data_file: Path,
    folder: Path,
    dimension: int,
    n_sis: int,
    held_out: list[int],
    sample_id: str,
    settings: dict,
    resume: bool,
) -> tuple[float, float, str]:
    folder.mkdir(parents=True, exist_ok=True)
    config_path = folder / "sisso.json"
    model_path = folder / "models" / f"test_dim_{dimension}_model_0.dat"
    if not (resume and model_path.exists()):
        dump_json(
            config_path,
            make_sissopp_config(data_file, dimension, n_sis, held_out, settings),
        )
        environment = os.environ.copy()
        environment["OMP_NUM_THREADS"] = "1"
        with (folder / "run.log").open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                [str(binary), "sisso.json"],
                cwd=folder,
                env=environment,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        if completed.returncode != 0 or not model_path.exists():
            raise RuntimeError(
                f"SISSO++ failed in {folder} (return code {completed.returncode}); see run.log"
            )
    return read_prediction(model_path, sample_id)


def run_nested(
    binary: Path,
    input_dir: Path,
    work_dir: Path,
    config_path: Path,
    dimensions: list[int],
    workers: int,
    resume: bool,
) -> None:
    config = load_json(config_path)
    settings = config["sissopp"]
    data_file = input_dir / "data.csv"
    samples = read_csv(input_dir / "sample_map.csv")
    if len(samples) < 4:
        raise ValueError("nested leave-one-out validation requires at least four materials")
    if workers < 1:
        raise ValueError("workers must be positive")
    n_sis_grid = [int(value) for value in settings["n_sis_grid"]]
    results_dir = work_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    for dimension in dimensions:
        prediction_rows: list[dict] = []
        selection_rows: list[dict] = []
        for outer in range(len(samples)):
            inner_results: dict[int, list[tuple[float, float]]] = {
                n_sis: [] for n_sis in n_sis_grid
            }
            futures = {}
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for n_sis in n_sis_grid:
                    for inner in range(len(samples)):
                        if inner == outer:
                            continue
                        folder = (
                            work_dir
                            / "runs"
                            / f"D{dimension}"
                            / f"outer_{outer:02d}"
                            / "inner"
                            / f"nsis_{n_sis}"
                            / f"holdout_{inner:02d}"
                        )
                        future = pool.submit(
                            run_one,
                            binary,
                            data_file,
                            folder,
                            dimension,
                            n_sis,
                            [outer, inner],
                            samples[inner]["sample_id"],
                            settings,
                            resume,
                        )
                        futures[future] = (n_sis, inner)
                for future in as_completed(futures):
                    n_sis, _ = futures[future]
                    observed, predicted, _ = future.result()
                    inner_results[n_sis].append((observed, predicted))

            scored: list[tuple[float, int]] = []
            for n_sis in n_sis_grid:
                values = inner_results[n_sis]
                if len(values) != len(samples) - 1:
                    raise RuntimeError("an inner fold is missing from model selection")
                inner_rmse = (
                    sum((predicted - observed) ** 2 for observed, predicted in values) / len(values)
                ) ** 0.5
                selection_rows.append(
                    {"outer_fold": outer, "n_sis": n_sis, "inner_RMSE_meV": inner_rmse}
                )
                scored.append((inner_rmse, n_sis))
            selected_n_sis = min(scored)[1]

            final_folder = (
                work_dir
                / "runs"
                / f"D{dimension}"
                / f"outer_{outer:02d}"
                / "final"
                / f"nsis_{selected_n_sis}"
            )
            observed, predicted, expression = run_one(
                binary,
                data_file,
                final_folder,
                dimension,
                selected_n_sis,
                [outer],
                samples[outer]["sample_id"],
                settings,
                resume,
            )
            expected_observed = float(samples[outer]["FWHM_meV"])
            if abs(observed - expected_observed) > 1e-7:
                raise RuntimeError(
                    f"outer fold {outer} prediction does not match its mapped target; check row ordering"
                )
            prediction_rows.append(
                {
                    "outer_fold": outer,
                    "sample_id": samples[outer]["sample_id"],
                    "material": samples[outer]["material"],
                    "observed_FWHM_meV": observed,
                    "predicted_FWHM_meV": predicted,
                    "error_meV": predicted - observed,
                    "selected_n_sis": selected_n_sis,
                    "model_expression": expression,
                }
            )
            print(
                f"D={dimension} outer={outer:02d} selected_n_sis={selected_n_sis} "
                f"prediction={predicted:.6f}",
                flush=True,
            )

        write_csv(results_dir / f"D{dimension}_inner_selection.csv", selection_rows)
        write_csv(results_dir / f"D{dimension}_outer_predictions.csv", prediction_rows)
        metrics = regression_metrics(
            [float(row["observed_FWHM_meV"]) for row in prediction_rows],
            [float(row["predicted_FWHM_meV"]) for row in prediction_rows],
        )
        metrics.update(
            {
                "dimension": dimension,
                "n_samples": len(samples),
                "n_sis_grid": n_sis_grid,
            }
        )
        dump_json(results_dir / f"D{dimension}_metrics.json", metrics)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary",
        type=Path,
        default=Path(os.environ["SISSOPP_BINARY"]) if "SISSOPP_BINARY" in os.environ else None,
        help="path to the pinned SISSO++ executable (or set SISSOPP_BINARY)",
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUTS)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--dimensions", type=parse_dimensions, default=[1, 2])
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--resume", action="store_true", help="reuse completed SISSO++ model files")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.binary is None:
        raise SystemExit("provide --binary or set SISSOPP_BINARY")
    binary = args.binary.expanduser().resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise SystemExit(f"SISSO++ executable is not available: {binary}")
    run_nested(
        binary,
        args.input_dir.resolve(),
        args.work_dir.resolve(),
        args.config.resolve(),
        args.dimensions,
        args.workers,
        args.resume,
    )
    print(f"Nested results written to {(args.work_dir / 'results').resolve()}")


if __name__ == "__main__":
    main()
