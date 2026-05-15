import pprint
from pathlib import Path
from typing import Union, List, Optional

import pandas as pd
import yaml
from pydantic import BaseModel

from c3p.datamodel import Config, Result, EvaluationResult, EvaluationExperiment, ResultSet
from c3p.learn import safe_name


def object_to_header(obj: Union[BaseModel, dict]):
    obj_yaml = yaml.dump(obj.model_dump() if isinstance(obj, BaseModel) else obj)
    obj_yaml = "\n".join([f"# {line}" for line in obj_yaml.split("\n")])
    return obj_yaml


def result_prog_path(r: Result, results_dir: Path):
    cn = safe_name(r.chemical_class.name)
    prog_dir = results_dir / "programs"
    prog_dir.mkdir(exist_ok=True, parents=True)
    return f"{prog_dir / cn}.py"

def write_program(result: Union[Result, EvaluationResult], result_dir: Path, prog_path=None, include_metadata=True, config_yaml: Union[str, Config]=None):
    """
    Write a program to disk.

    Args:
        r:
        result_dir:
        prog_path:
        include_metadata:
        config_yaml:

    Returns:
    """
    #if isinstance(config_yaml, Config):
    #    config_yaml = object_to_header(config_yaml)
    if isinstance(result, EvaluationResult):
        r = result.train_results.best_result
    else:
        r = result
    if prog_path is None:
        prog_path = result_prog_path(r, result_dir)
    with open(prog_path, "w") as f:
        f.write(f'"""\nClassifies: {r.chemical_class.id} {r.chemical_class.name}\n"""')
        f.write("\n")
        f.write(r.code)
        if include_metadata:
            f.write("\n\n\n__metadata__ = ")
            #f.write("\n\n# Classifier Metadata\n# ---\n")
            r_obj = r.model_dump()
            del r_obj["code"]
            cc = r_obj["chemical_class"]
            #del cc["instances"]
            #del cc["negative_instances"]
            del r_obj["true_positives"]
            del r_obj["false_positives"]
            del r_obj["true_negatives"]
            del r_obj["false_negatives"]
            pp = pprint.PrettyPrinter(indent=4, sort_dicts=False)
            f.write(pp.pformat(r_obj))


def write_programs(results: List[Union[Result, EvaluationResult]], result_dir: Path, config_yaml: Union[str, Config]=None):
    for r in results:
        write_program(r, result_dir, config_yaml)

def write_programs_for_commits(results: List[EvaluationResult],  result_dir: Path):
    """
    Write each iteration of the program to a separate file, to be staged for commit

    Args:
        results:
        result_dir:

    Returns:

    """
    staged_dir = result_dir / "staged"
    staged_dir.mkdir(exist_ok=True, parents=True)
    with open(staged_dir / "commit.sh", "w") as script_f:
        script_f.write("#!/bin/bash\n")

        for r in results:
            write_program_for_commits(r, result_dir, script_f)

def write_program_for_commits(result: EvaluationResult, result_dir: Path, script_f):
    """
    Write a program to disk, and add the file to a commit script.

    Args:
        result:
        result_dir:
        script_f: file object

    Returns:

    """
    expt = result_dir.name.replace("-undef", "")
    base_path = result_dir / "staged"
    learned_dir = Path("learned") / expt
    script_f.write(f"mkdir -p {learned_dir}\n")
    base_path.mkdir(exist_ok=True, parents=True)
    for r in result.train_results.results:
        cn = safe_name(r.chemical_class.name)
        path = base_path / f"{cn}_{r.attempt}.py"
        write_program(r, result_dir, prog_path=path, include_metadata=False)
        script_f.write(f"cp {path} {learned_dir}/{cn}.py\n")
        script_f.write(f"git add {learned_dir}/{cn}.py\n")
        msg_path = base_path / f"{cn}_{r.attempt}.txt"
        with open(msg_path, "w") as f:
            f.write(f"Attempt {r.attempt+1} at {r.chemical_class.name}, f1={r.f1}\n\n")
            if not r.success:
                f.write(f"(currently fails to compile, but committing anyway\n\n")
            f.write(f"{r.reasoning}\n")
            if r.message:
                f.write(f"\n# PREVIOUS ATTEMPT:\n\n{r.message}\n")
            r_obj = r.model_dump()
            tps = r_obj.get("true_positives")
            if not tps:
                tps = []
            del r_obj["code"]
            del r_obj["message"]
            del r_obj["reasoning"]
            r_obj["sample_true_positives"] = tps[:5]
            del r_obj["true_positives"]
            del r_obj["false_positives"]
            code_stats = r_obj.get("code_statistics", {})
            if code_stats and "indent_by_line" in code_stats:
                del code_stats["indent_by_line"]
            f.write("# INFORMATION:\n\n")
            f.write("\n\n---\n")
            f.write(yaml.dump(r_obj, sort_keys=False))
        script_f.write(f"git commit -F {msg_path} {learned_dir}/{cn}.py\n")
    br = result.train_results.best_result
    cn = safe_name(br.chemical_class.name)
    path = base_path / f"{cn}_best.py"
    write_program(br, result_dir, prog_path=path, include_metadata=False)
    msg_path = base_path / f"{cn}_best.txt"
    with open(msg_path, "w") as f:
        r = br
        f.write(f"Reverting to best solution for {br.chemical_class.name}, f1={br.f1} iteration {br.attempt}\n")
        r_obj = r.model_dump()
        del r_obj["code"]
        del r_obj["reasoning"]
        del r_obj["true_positives"]
        del r_obj["false_positives"]
        code_stats = r_obj.get("code_statistics", {})
        if code_stats and "indent_by_line" in code_stats:
            del code_stats["indent_by_line"]
        f.write("\n\n---\n")
        f.write(yaml.dump(r_obj, sort_keys=False))
    script_f.write(f"cp {path} {learned_dir}/{cn}.py\n")
    script_f.write(f"git commit -F {msg_path} {learned_dir}/{cn}.py\n")

def result_as_dict(e: Union[Result, EvaluationResult, ResultSet]) -> dict:
    if isinstance(e, ResultSet):
        e = e.best_result
    if isinstance(e, Result):
        d = e.model_dump()
        del d["code"]
        del d["config"]
        del d["message"]
        d["chemical_class_id"] = e.chemical_class.id
        d["chemical_class"] = e.chemical_class.name
        return d
    if not isinstance(e, EvaluationResult):
        raise AssertionError(f"Expected Result or EvaluationResult, got {type(e)}")
    r = e.test_result
    d = r.model_dump()
    del d["code"]
    d["chemical_class"] = r.chemical_class.name
    del d["config"]
    #d["attempt"] = e.train_results.best_result.attempt
    # d["num_instances"] = len(r.chemical_class.positive_instances)
    for k in ["f1", "precision", "accuracy", "attempt"]:
        d[f"train_{k}"] = getattr(e.train_results.best_result, k, None)
    return d

def xxxwrite_eval_results(expt: EvaluationExperiment, results_dir: Path, f1_threshold: Optional[float]=None) -> pd.DataFrame:
    """
    Write the results of an evaluation experiment to disk.

    Args:
        expt:
        results_dir:

    Returns: dataframe of filtered results

    """
    if f1_threshold is None:
        f1_threshold = expt.config.f1_threshold
    expt_path = results_dir / "experiment.yaml"
    with open(expt_path, "w") as f:
        f.write(yaml.dump(expt.model_dump()))
    expt_path_json = results_dir / "experiment.json"
    with open(expt_path_json, "w") as f:
        f.write(expt.model_dump_json(indent=2))
    eval_results = expt.evaluation_results
    for r in eval_results:
        write_program(r, results_dir, expt.config)
    config = expt.config

    df_filtered = pd.DataFrame(
        [result_as_dict(r) for r in eval_results if r.train_results.best_result.f1 > f1_threshold])
    df_filtered["num_test_instances"] = df_filtered["num_true_positives"] + df_filtered["num_false_negatives"]
    df_filtered = df_filtered[df_filtered["num_test_instances"] > 0]
    df_filtered.to_csv(results_dir / "filtered_results.csv", index=False)
    summary_stats = calculate_metrics_pandas(df_filtered)
    summary_stats.to_csv(results_dir / "outcomes_macro.csv", header=False)
    summary_stats = df_filtered.select_dtypes(include=['int64', 'float64']).mean()
    summary_stats.to_csv(results_dir / "outcomes_micro.csv", header=False)
    desc = df_filtered.describe().transpose()
    desc.to_csv(results_dir / "summary_stats_desc.csv")
    return df_filtered


