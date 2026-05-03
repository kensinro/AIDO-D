"""Generate size-matched random gene-set baseline for AIDO-D.

The script compares the discriminability of a selected pathway against random
sets of the same gene-set size. This is used to assess whether observed D is
above a stochastic baseline.
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List

import numpy as np
import pandas as pd

sys.path.append(os.path.dirname(__file__))
from aido_d_common import (
    SurvivalColumns,
    align_scores_survival,
    ensure_dir,
    evaluate_all_observables,
    filter_primary_tumor_samples,
    load_expression_matrix,
    load_survival_table,
    parse_gmt,
    zscore_genes,
)


def random_gene_set_scores(z_expr: pd.DataFrame, size: int, trials: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    genes: List[str] = list(z_expr.index)
    scores = {}
    for t in range(trials):
        sampled = rng.choice(genes, size=size, replace=False)
        scores[f"random_{t+1:04d}"] = z_expr.loc[sampled].mean(axis=0, skipna=True)
    return pd.DataFrame(scores)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True)
    parser.add_argument("--survival", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--target-pathway", required=True)
    parser.add_argument("--sample-col", required=True)
    parser.add_argument("--time-col", required=True)
    parser.add_argument("--event-col", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--orientation", default="genes_by_samples", choices=["genes_by_samples", "samples_by_genes"])
    parser.add_argument("--sample-type", default="01")
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    ensure_dir(args.outdir)
    expr = load_expression_matrix(args.expression, orientation=args.orientation)
    expr = filter_primary_tumor_samples(expr, sample_type=args.sample_type)
    z = zscore_genes(expr)
    survival = load_survival_table(args.survival, SurvivalColumns(args.sample_col, args.time_col, args.event_col))
    gene_sets = parse_gmt(args.gmt)
    if args.target_pathway not in gene_sets:
        raise ValueError(f"target pathway not found in GMT: {args.target_pathway}")
    overlap = [g for g in gene_sets[args.target_pathway] if g in set(z.index)]
    if len(overlap) < 5:
        raise ValueError(f"Not enough overlapping genes for {args.target_pathway}: {len(overlap)}")

    rand_scores = random_gene_set_scores(z, size=len(overlap), trials=args.trials, seed=args.seed)
    table = align_scores_survival(rand_scores, survival)
    ranking = evaluate_all_observables(table, rand_scores.columns)
    ranking.to_csv(os.path.join(args.outdir, f"random_baseline_{args.target_pathway}.csv"), index=False)

    summary = {
        "target_pathway": args.target_pathway,
        "gene_set_size": len(overlap),
        "random_trials": args.trials,
        "random_median_D": ranking["D"].median(),
        "random_95th_percentile_D": ranking["D"].quantile(0.95),
        "random_max_D": ranking["D"].max(),
    }
    pd.DataFrame([summary]).to_csv(os.path.join(args.outdir, f"random_baseline_summary_{args.target_pathway}.csv"), index=False)
    print(f"Saved random baseline outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
