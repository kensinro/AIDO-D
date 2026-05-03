"""Run the AIDO-D discriminability pipeline from a YAML config file."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

import yaml


def run(cmd: list[str]) -> None:
    print("\n$ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    inp = cfg["inputs"]
    cols = cfg["columns"]
    settings = cfg.get("settings", {})
    outdir = cfg["outputs"]["output_dir"]
    os.makedirs(outdir, exist_ok=True)

    script_dir = os.path.dirname(__file__)
    py = sys.executable

    run([
        py,
        os.path.join(script_dir, "01_compute_pathway_discriminability.py"),
        "--expression", inp["expression_matrix"],
        "--survival", inp["survival_table"],
        "--gmt", inp["hallmark_gmt"],
        "--sample-col", cols["sample_id"],
        "--time-col", cols["survival_time"],
        "--event-col", cols["survival_event"],
        "--outdir", outdir,
        "--orientation", settings.get("expression_orientation", "genes_by_samples"),
        "--sample-type", str(settings.get("sample_type_filter", "01")),
        "--min-genes", str(settings.get("min_genes_per_pathway", 5)),
        "--make-km",
    ])

    ranking = os.path.join(outdir, "AIDO_D_pathway_D_ranking.csv")
    score_table = os.path.join(outdir, "AIDO_D_pathway_scores_with_survival.csv")

    run([
        py,
        os.path.join(script_dir, "04_discriminability_guided_selection.py"),
        "--ranking", ranking,
        "--outdir", outdir,
        "--top-k", "10",
        "--threshold-D", "2",
    ])

    run([
        py,
        os.path.join(script_dir, "03_low_high_D_figure_panels.py"),
        "--score-table", score_table,
        "--ranking", ranking,
        "--outdir", outdir,
    ])

    print("\nAIDO-D pipeline completed.")


if __name__ == "__main__":
    main()
