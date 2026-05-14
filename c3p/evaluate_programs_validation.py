"""Evaluate all program files against validation examples in a dataset.

This script loads each program from a directory (default: c3p/programs), maps it to a
CHEBI class (prefer metadata id, fallback to safe-name/file-stem match), evaluates on
validation examples, and writes per-class metrics.
"""

from __future__ import annotations

import argparse
import ast
import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from tqdm import tqdm

from c3p.datamodel import Dataset
from c3p.learn import (
    evaluate_program,
    get_positive_and_negative_validate_instances,
)

logger = logging.getLogger(__name__)


def _extract_metadata_from_code(code: str) -> Dict[str, Any]:
    """Parse __metadata__ dict from a program file without executing code."""
    try:
        module = ast.parse(code)
    except SyntaxError as e:
        raise ValueError(f"Failed to parse program for metadata extraction: {e}")

    if (
        len(
            [
                target
                for node in module.body
                if isinstance(node, ast.Assign)
                for target in node.targets
                if isinstance(target, ast.Name) and target.id == "__metadata__"
            ]
        )
        > 1
    ):
        raise ValueError("Multiple assignments to __metadata__ found in program")

    for node in module.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "__metadata__":
                    try:
                        value = ast.literal_eval(node.value)
                    except Exception as e:
                        raise ValueError(f"Could not literal-eval __metadata__: {e}")
                    if isinstance(value, dict):
                        return value
                    raise ValueError(f"__metadata__ is not a dict: {value}")
    raise ValueError(f"__metadata__ not found in program: {code}")


def _find_dataset_class(dataset: Dataset, program_path: Path, metadata: Dict[str, Any]):
    """Find dataset class for a program via metadata CHEBI id, then by safe name."""
    chebi_id = metadata.get("chemical_class", {}).get("id")
    if chebi_id:
        try:
            return dataset.get_chemical_class_by_id(chebi_id)
        except ValueError:
            raise ValueError(
                f"Class id {chebi_id} from {program_path.name} not found in dataset"
            )
    else:
        raise ValueError(
            f"No CHEBI id in metadata for {program_path.name}, falling back to filename match"
        )


def _evaluate_one_program(program_path: Path, dataset: Dataset) -> Dict[str, Any]:
    code = program_path.read_text()
    try:
        metadata = _extract_metadata_from_code(code)
    except ValueError as e:
        raise ValueError(
            {
                "program": program_path.name,
                "status": "error",
                "error": str(e),
            }
        )
    cls = _find_dataset_class(dataset, program_path, metadata)

    if cls is None:
        raise ValueError(
            f"No matching class found in dataset for program {program_path.name}"
        )

    pos, neg = get_positive_and_negative_validate_instances(cls, dataset)
    if not pos and not neg:
        raise ValueError(
            f"No validation examples for class {cls.id} ({cls.name}) in program {program_path.name}"
        )

    try:
        result = evaluate_program(
            code,
            cls.lite_copy(),
            pos,
            neg,
            threshold=-0.01,  # some dummy value, not used
        )
    except Exception as e:
        raise ValueError(
            f"Error occurred while evaluating program {program_path.name}: {e}"
        )

    return {
        "program": program_path.name,
        "chebi_id": cls.id,
        "class_name": cls.name,
        "num_validate_pos": len(pos),
        "num_validate_neg": len(neg),
        "status": "ok" if result.success else "runtime_error",
        "success": result.success,
        "precision": result.precision,
        "recall": result.recall,
        "f1": result.f1,
        "accuracy": result.accuracy,
        "negative_predictive_value": result.negative_predictive_value,
        "num_true_positives": result.num_true_positives,
        "num_false_positives": result.num_false_positives,
        "num_true_negatives": result.num_true_negatives,
        "num_false_negatives": result.num_false_negatives,
        "error": result.error,
    }


def evaluate_programs(
    dataset_path: Path,
    program_dir: Path,
    output_csv: Path,
    output_json: Optional[Path] = None,
) -> pd.DataFrame:
    with dataset_path.open("r") as f:
        dataset = Dataset.model_validate_json(f.read())

    programs = sorted(
        p for p in program_dir.glob("*.py") if not p.name.startswith("__")
    )
    rows: List[Dict[str, Any]] = []

    # Load existing results if any, to skip already-evaluated programs
    done_programs = set()
    if output_json and output_json.exists():
        try:
            rows = json.loads(output_json.read_text())
            done_programs = {r.get("program") for r in rows if r.get("program")}
            logger.info(
                "Loaded %d existing results, will skip %d programs",
                len(rows),
                len(done_programs),
            )
        except Exception as e:
            logger.warning("Failed to load existing output_json %s: %s", output_json, e)
            rows = []
            done_programs = set()

    # Ensure output directory exists before writing iteratively
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    if output_json:
        output_json.parent.mkdir(parents=True, exist_ok=True)

    for program_path in tqdm(programs, desc="Programs"):
        if program_path.name in done_programs:
            logger.info("Skipping %s (already evaluated)", program_path.name)
            continue

        logger.info("Evaluating %s", program_path.name)
        row = _evaluate_one_program(program_path, dataset)
        rows.append(row)
        done_programs.add(row.get("program"))

        # Persist results after each program so progress is saved.
        try:
            df_partial = pd.DataFrame(rows)
            df_partial.to_csv(output_csv, index=False)
            if output_json:
                output_json.write_text(json.dumps(rows, indent=2))
        except Exception as e:
            raise ValueError("Failed to write incremental results: %s", e)

    df = pd.DataFrame(rows)
    return df


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Evaluate all c3p program files against dataset validation examples."
    )
    parser.add_argument(
        "--dataset",
        required=True,
        type=Path,
        help="Path to dataset JSON file",
    )
    parser.add_argument(
        "--program-dir",
        type=Path,
        default=Path("c3p/programs"),
        help="Directory containing classifier program files",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/validation_program_eval.csv"),
        help="Where to write CSV results",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("results/validation_program_eval.json"),
        help="Where to write JSON results",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="count",
        default=0,
        help="Increase logging verbosity",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    level = logging.WARNING
    if args.verbose == 1:
        level = logging.INFO
    elif args.verbose >= 2:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(levelname)s:%(name)s:%(message)s")

    df = evaluate_programs(
        dataset_path=args.dataset,
        program_dir=args.program_dir,
        output_csv=args.output_csv,
        output_json=args.output_json,
    )

    print(f"Evaluated {len(df)} programs")
    if "status" in df.columns:
        print(df["status"].value_counts(dropna=False).to_string())
    metric_cols = [
        c for c in ["precision", "recall", "f1", "accuracy"] if c in df.columns
    ]
    if metric_cols:
        ok_df = df[df["status"] == "ok"]
        if not ok_df.empty:
            print("\nSummary of successful evaluations:")
            print(ok_df[metric_cols].describe().to_string())


if __name__ == "__main__":
    main()
