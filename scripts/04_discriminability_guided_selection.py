"""Discriminability-guided observable selection.

This script ranks candidate observables by D and exports selected observables under
threshold or top-k rules. It supports the manuscript claim that observable design can
be treated as an outcome-oriented selection problem.
"""

from __future__ import annotations

import argparse
import os

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--threshold-D", type=float, default=None)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ranking = pd.read_csv(args.ranking).sort_values("D", ascending=False)

    top = ranking.head(args.top_k)
    top.to_csv(os.path.join(args.outdir, f"top_{args.top_k}_observables_by_D.csv"), index=False)

    if args.threshold_D is not None:
        selected = ranking[ranking["D"] >= args.threshold_D]
        selected.to_csv(os.path.join(args.outdir, f"observables_D_ge_{args.threshold_D}.csv"), index=False)

    summary = {
        "n_observables": int(ranking.shape[0]),
        "top_k": args.top_k,
        "max_D": float(ranking["D"].max()),
        "median_D": float(ranking["D"].median()),
        "n_D_ge_1": int((ranking["D"] >= 1).sum()),
        "n_D_ge_2": int((ranking["D"] >= 2).sum()),
        "n_D_ge_4": int((ranking["D"] >= 4).sum()),
    }
    pd.DataFrame([summary]).to_csv(os.path.join(args.outdir, "selection_summary.csv"), index=False)
    print(f"Saved selection outputs to: {args.outdir}")


if __name__ == "__main__":
    main()
