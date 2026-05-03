"""Create low-D vs high-D Kaplan-Meier figure panels from computed pathway scores.

Input:
    AIDO_D_pathway_scores_with_survival.csv from script 01.
    AIDO_D_pathway_D_ranking.csv from script 01.

Output:
    KM plots for representative low-D and high-D observables.
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.append(os.path.dirname(__file__))
from aido_d_common import ensure_dir, save_km_plot


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-table", required=True)
    parser.add_argument("--ranking", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--low-observable", default=None)
    parser.add_argument("--high-observable", default=None)
    args = parser.parse_args()

    ensure_dir(args.outdir)
    table = pd.read_csv(args.score_table, index_col=0)
    ranking = pd.read_csv(args.ranking)

    if args.high_observable:
        high = args.high_observable
    else:
        high = ranking.sort_values("D", ascending=False).iloc[0]["observable"]

    if args.low_observable:
        low = args.low_observable
    else:
        valid = ranking.dropna(subset=["D"])
        low = valid.sort_values("D", ascending=True).iloc[0]["observable"]

    save_km_plot(table, low, os.path.join(args.outdir, "KM_low_discriminability_observable.png"), title=f"Low-D observable: {low}")
    save_km_plot(table, high, os.path.join(args.outdir, "KM_high_discriminability_observable.png"), title=f"High-D observable: {high}")

    pd.DataFrame([{"low_observable": low, "high_observable": high}]).to_csv(
        os.path.join(args.outdir, "selected_low_high_observables.csv"), index=False
    )
    print(f"Selected low-D: {low}")
    print(f"Selected high-D: {high}")


if __name__ == "__main__":
    main()
