"""Compute pathway-level discriminability for AIDO-D.

This script constructs Hallmark pathway observables from a gene-expression matrix,
stratifies patients by each pathway score, performs Kaplan-Meier/log-rank analysis,
and ranks observables by D = -log10(p).
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.append(os.path.dirname(__file__))
from aido_d_common import (
    SurvivalColumns,
    align_scores_survival,
    compute_pathway_scores,
    ensure_dir,
    evaluate_all_observables,
    filter_primary_tumor_samples,
    load_expression_matrix,
    load_survival_table,
    parse_gmt,
    save_km_plot,
    zscore_genes,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", required=True)
    parser.add_argument("--survival", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--sample-col", required=True)
    parser.add_argument("--time-col", required=True)
    parser.add_argument("--event-col", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--orientation", default="genes_by_samples", choices=["genes_by_samples", "samples_by_genes"])
    parser.add_argument("--sample-type", default="01")
    parser.add_argument("--min-genes", type=int, default=5)
    parser.add_argument("--make-km", action="store_true")
    args = parser.parse_args()

    ensure_dir(args.outdir)

    expr = load_expression_matrix(args.expression, orientation=args.orientation)
    expr = filter_primary_tumor_samples(expr, sample_type=args.sample_type)
    z = zscore_genes(expr)

    survival = load_survival_table(
        args.survival,
        SurvivalColumns(args.sample_col, args.time_col, args.event_col),
    )
    gene_sets = parse_gmt(args.gmt)
    scores = compute_pathway_scores(z, gene_sets, min_genes=args.min_genes)
    table = align_scores_survival(scores, survival)

    ranking = evaluate_all_observables(table, scores.columns)
    ranking.to_csv(os.path.join(args.outdir, "AIDO_D_pathway_D_ranking.csv"), index=False)
    table.to_csv(os.path.join(args.outdir, "AIDO_D_pathway_scores_with_survival.csv"))

    if args.make_km and not ranking.empty:
        top = ranking.iloc[0]["observable"]
        save_km_plot(table, top, os.path.join(args.outdir, "KM_top_D_pathway.png"), title=f"Top-D pathway: {top}")
        low = ranking.dropna().sort_values("D", ascending=True).iloc[0]["observable"]
        save_km_plot(table, low, os.path.join(args.outdir, "KM_low_D_pathway.png"), title=f"Low-D pathway: {low}")

    print(f"Saved outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
