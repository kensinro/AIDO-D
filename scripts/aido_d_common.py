"""
Common utilities for the AIDO-D discriminability pipeline.

Core idea:
    D(O_k) = -log10(p_k)
where p_k is the log-rank p-value comparing survival curves after stratifying
patients by observable O_k.
"""

from __future__ import annotations

import os
import re
import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test


@dataclass
class SurvivalColumns:
    sample_id: str
    survival_time: str
    survival_event: str


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def read_table(path: str, index_col: Optional[int | str] = None) -> pd.DataFrame:
    """Read CSV/TSV robustly based on file extension or delimiter sniffing."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    lower = path.lower()
    if lower.endswith(".csv"):
        return pd.read_csv(path, index_col=index_col)
    return pd.read_csv(path, sep="\t", index_col=index_col)


def normalize_tcga_sample_id(sample_id: str, length: int = 15) -> str:
    """Normalize TCGA barcode to patient/sample-level prefix.

    length=15 keeps TCGA-XX-YYYY-01 style sample-type information.
    length=12 keeps patient-level barcode.
    """
    s = str(sample_id).replace(".", "-")
    return s[:length]


def load_expression_matrix(path: str, orientation: str = "genes_by_samples") -> pd.DataFrame:
    """Load expression matrix and return genes x samples."""
    df = read_table(path, index_col=0)
    if orientation == "samples_by_genes":
        df = df.T
    elif orientation != "genes_by_samples":
        raise ValueError("orientation must be genes_by_samples or samples_by_genes")

    df.index = df.index.astype(str)
    df.columns = [normalize_tcga_sample_id(c, length=15) for c in df.columns]
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def filter_primary_tumor_samples(expr: pd.DataFrame, sample_type: Optional[str] = "01") -> pd.DataFrame:
    """Keep TCGA primary tumor samples if sample type information is present."""
    if sample_type is None:
        return expr
    keep = []
    for c in expr.columns:
        parts = str(c).split("-")
        keep.append(len(parts) >= 4 and parts[3][:2] == sample_type)
    if any(keep):
        return expr.loc[:, keep]
    return expr


def zscore_genes(expr: pd.DataFrame) -> pd.DataFrame:
    """Gene-wise z-score across samples."""
    means = expr.mean(axis=1, skipna=True)
    stds = expr.std(axis=1, skipna=True).replace(0, np.nan)
    z = expr.sub(means, axis=0).div(stds, axis=0)
    return z.dropna(axis=0, how="all")


def load_survival_table(path: str, cols: SurvivalColumns) -> pd.DataFrame:
    df = read_table(path)
    missing = [c for c in [cols.sample_id, cols.survival_time, cols.survival_event] if c not in df.columns]
    if missing:
        raise ValueError(f"Survival table missing columns: {missing}. Available: {list(df.columns)}")
    out = df[[cols.sample_id, cols.survival_time, cols.survival_event]].copy()
    out.columns = ["sample", "time", "event"]
    out["sample"] = out["sample"].map(lambda x: normalize_tcga_sample_id(x, length=15))
    out["time"] = pd.to_numeric(out["time"], errors="coerce")
    out["event"] = pd.to_numeric(out["event"], errors="coerce")
    out = out.dropna(subset=["sample", "time", "event"])
    out = out[out["time"] > 0]
    out["event"] = out["event"].astype(int)
    return out.drop_duplicates("sample")


def parse_gmt(path: str) -> Dict[str, List[str]]:
    """Parse GMT gene-set file into {set_name: [genes]} dictionary."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    gene_sets: Dict[str, List[str]] = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                name = parts[0]
                genes = [g.strip() for g in parts[2:] if g.strip()]
                gene_sets[name] = genes
    return gene_sets


def compute_pathway_scores(z_expr: pd.DataFrame, gene_sets: Dict[str, List[str]], min_genes: int = 5) -> pd.DataFrame:
    """Average z-scored expression within each gene set; returns samples x pathways."""
    available = set(z_expr.index)
    scores = {}
    for name, genes in gene_sets.items():
        overlap = [g for g in genes if g in available]
        if len(overlap) >= min_genes:
            scores[name] = z_expr.loc[overlap].mean(axis=0, skipna=True)
    return pd.DataFrame(scores)


def align_scores_survival(scores: pd.DataFrame, survival: pd.DataFrame) -> pd.DataFrame:
    """Return one table with sample, time, event, and scores."""
    s = survival.set_index("sample")
    common = scores.index.intersection(s.index)
    merged = s.loc[common, ["time", "event"]].join(scores.loc[common], how="inner")
    merged = merged.dropna(subset=["time", "event"])
    return merged


def median_split(values: pd.Series) -> pd.Series:
    med = values.median(skipna=True)
    return pd.Series(np.where(values <= med, "Low", "High"), index=values.index)


def logrank_discriminability(table: pd.DataFrame, score_col: str) -> Tuple[float, float, int, int]:
    """Compute p-value and D after median split of score_col."""
    sub = table[["time", "event", score_col]].dropna()
    if sub.shape[0] < 10:
        return np.nan, np.nan, 0, 0
    groups = median_split(sub[score_col])
    low = sub[groups == "Low"]
    high = sub[groups == "High"]
    if low.empty or high.empty:
        return np.nan, np.nan, low.shape[0], high.shape[0]
    res = logrank_test(
        low["time"],
        high["time"],
        event_observed_A=low["event"],
        event_observed_B=high["event"],
    )
    p = float(res.p_value)
    if p <= 0:
        d = np.inf
    else:
        d = -math.log10(p)
    return p, d, low.shape[0], high.shape[0]


def evaluate_all_observables(table: pd.DataFrame, score_cols: Iterable[str]) -> pd.DataFrame:
    rows = []
    for col in score_cols:
        p, d, n_low, n_high = logrank_discriminability(table, col)
        rows.append({"observable": col, "D": d, "p_value": p, "n_low": n_low, "n_high": n_high})
    out = pd.DataFrame(rows)
    return out.sort_values("D", ascending=False, na_position="last").reset_index(drop=True)


def save_km_plot(table: pd.DataFrame, score_col: str, output_png: str, title: Optional[str] = None) -> None:
    import matplotlib.pyplot as plt

    sub = table[["time", "event", score_col]].dropna()
    groups = median_split(sub[score_col])
    p, d, n_low, n_high = logrank_discriminability(table, score_col)

    fig, ax = plt.subplots(figsize=(6.2, 4.8), dpi=300)
    kmf = KaplanMeierFitter()
    for label in ["Low", "High"]:
        grp = sub[groups == label]
        kmf.fit(grp["time"], grp["event"], label=f"{label} (n={grp.shape[0]})")
        kmf.plot_survival_function(ax=ax, ci_show=False)
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.text(0.05, 0.08, f"Log-rank p = {p:.2e}\nD = {d:.2f}", transform=ax.transAxes)
    ax.set_title(title or score_col)
    fig.tight_layout()
    fig.savefig(output_png)
    plt.close(fig)
