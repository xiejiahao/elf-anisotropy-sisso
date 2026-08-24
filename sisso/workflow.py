"""Shared validation, parsing, and metric helpers for the SISSO++ workflow."""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable, Sequence


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def read_header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return next(csv.reader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fieldnames: Sequence[str] | None = None) -> None:
    if not rows and fieldnames is None:
        raise ValueError(f"cannot infer columns for empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(fieldnames or rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def dump_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_bool(value: str, *, field: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"{field} must contain boolean values; found {value!r}")


def require_columns(path: Path, actual: Iterable[str], required: Iterable[str]) -> None:
    missing = sorted(set(required) - set(actual))
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")


def require_unique(values: Sequence[str], label: str) -> None:
    duplicates = sorted({value for value in values if values.count(value) > 1})
    if duplicates:
        preview = ", ".join(duplicates[:5])
        raise ValueError(f"{label} contains duplicate values: {preview}")


def finite_float(value: str, *, label: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"{label} is not numeric: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{label} is not finite: {value!r}")
    return parsed


def regression_metrics(observed: Sequence[float], predicted: Sequence[float]) -> dict[str, float]:
    if len(observed) != len(predicted) or len(observed) < 2:
        raise ValueError("observed and predicted must have the same length of at least two")
    errors = [prediction - observation for observation, prediction in zip(observed, predicted)]
    mean_observed = sum(observed) / len(observed)
    sse = sum(error * error for error in errors)
    sst = sum((value - mean_observed) ** 2 for value in observed)
    if sst == 0.0:
        raise ValueError("Q2 is undefined for a constant target")
    return {
        "RMSE_meV": math.sqrt(sse / len(errors)),
        "MAE_meV": sum(abs(error) for error in errors) / len(errors),
        "Q2_LOO": 1.0 - sse / sst,
    }


def read_prediction(model_path: Path, sample_id: str) -> tuple[float, float, str]:
    expression = ""
    for line in model_path.read_text(encoding="utf-8").splitlines():
        if line.startswith("# c0"):
            expression = line[2:].strip()
        if line.startswith("#") or not line.strip():
            continue
        fields = [field.strip() for field in line.split(",")]
        if fields[0] == sample_id:
            if len(fields) < 3:
                raise ValueError(f"malformed prediction row in {model_path}: {line}")
            return float(fields[1]), float(fields[2]), expression
    raise ValueError(f"sample {sample_id!r} was not found in {model_path}")


def expression_feature_ids(expression: str, known_ids: Iterable[str]) -> list[str]:
    known = set(known_ids)
    tokens = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expression)
    return [token for token in tokens if token in known]

