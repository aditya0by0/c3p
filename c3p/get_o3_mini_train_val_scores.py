"""Evaluate classifier program files under learned/o3-mini against validation examples.

Scans program files in a directory (default: learned/o3-mini), extracts __metadata__
safely without executing code, maps each program to a CHEBI class (prefer metadata id,
then metadata chemical_class name, then filename/safe-name match), evaluates on
validation examples, and writes per-class metrics incrementally to CSV/JSON.
Skips programs already present in an existing output JSON and persists progress after
each program so partial results are saved. Returns a pandas DataFrame of results.
"""

from __future__ import annotations

import argparse
import ast
import json
import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from tqdm import tqdm

from c3p.datamodel import Dataset
from c3p.learn import (
    evaluate_program,
    get_positive_and_negative_validate_instances,
    safe_name,
)

logger = logging.getLogger(__name__)
CHEBI_ID_RE = re.compile(r"CHEBI:\d+")


def _extract_chebi_ids_from_header(code: str) -> List[str]:
    """Extract all CHEBI ids from the file header comments or module docstring."""
    try:
        module = ast.parse(code)
    except SyntaxError:
        return []

    chebi_ids = []
    docstring = ast.get_docstring(module)
    if docstring:
        chebi_ids.extend(CHEBI_ID_RE.findall(docstring))

    for line in code.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            chebi_ids.extend(CHEBI_ID_RE.findall(stripped))
            continue
        break

    unique_chebi_ids = sorted(set(chebi_ids))
    if len(unique_chebi_ids) > 1:
        raise ValueError(
            f"Multiple distinct CHEBI ids found in program header: {', '.join(unique_chebi_ids)}"
        )
    return unique_chebi_ids


def _find_dataset_class(dataset: Dataset, program_path: Path):
    """Find dataset class for a program via metadata id, then name, then filename."""
    header_chebi_ids = _extract_chebi_ids_from_header(program_path.read_text())

    if header_chebi_ids:
        header_chebi_id = header_chebi_ids[0]
        try:
            return dataset.get_chemical_class_by_id(header_chebi_id)
        except ValueError as e:
            raise ValueError(
                f"Class id {header_chebi_id} from header of {program_path.name} not found in dataset"
            ) from e

    program_safe_name = safe_name(program_path.stem)
    for chemical_class in dataset.classes:
        if safe_name(chemical_class.name) == program_safe_name:
            return chemical_class

    raise ValueError(
        f"No matching class found in dataset for program {program_path.name}"
    )


def _evaluate_one_program(program_path: Path, dataset: Dataset) -> Dict[str, Any]:
    code = program_path.read_text()

    cls = _find_dataset_class(dataset, program_path)
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
            threshold=-0.01,  # dummy value, not used
        )
    except Exception as e:
        raise ValueError(
            {
                "program": program_path.name,
                "chebi_id": cls.id,
                "class_name": cls.name,
                "status": "error",
                "error": f"Error occurred while evaluating program {program_path.name}: {e}",
            }
        )

    val_metrics = {
        "num_pos": len(pos),
        "num_neg": len(neg),
        "status": "ok" if result.success else "runtime_error",
        "success": result.success,
        "precision": getattr(result, "precision", None),
        "recall": getattr(result, "recall", None),
        "f1": getattr(result, "f1", None),
        "accuracy": getattr(result, "accuracy", None),
        "negative_predictive_value": getattr(result, "negative_predictive_value", None),
        "num_true_positives": getattr(result, "num_true_positives", None),
        "num_false_positives": getattr(result, "num_false_positives", None),
        "num_true_negatives": getattr(result, "num_true_negatives", None),
        "num_false_negatives": getattr(result, "num_false_negatives", None),
        "error": getattr(result, "error", ""),
    }

    return {
        "program": program_path.name,
        "chebi_id": cls.id,
        "class_name": cls.name,
        "val": val_metrics,
        "train": {},
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

        try:
            df_partial = pd.DataFrame(rows)
            df_partial.to_csv(output_csv, index=False)
            if output_json:
                output_json.write_text(json.dumps(rows, indent=2))
        except Exception as e:
            raise ValueError(f"Failed to write incremental results: {e}") from e

    return pd.DataFrame(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Evaluate all learned/o3-mini program files against dataset validation examples."
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
        default=Path("learned/o3-mini"),
        help="Directory containing learned/o3-mini program files",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results/o3_mini_validation_program_eval.csv"),
        help="Where to write CSV results",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("results/o3_mini_validation_program_eval.json"),
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
