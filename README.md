# AIDO-D

This repository accompanies the manuscript:
"Discriminability reveals intrinsic limits of observable-based inference in cancer systems"
A framework for outcome-oriented evaluation of molecular observables based on discriminability (D), revealing intrinsic limits of observable-based inference in cancer systems.

---

## Overview

AIDO-D introduces **discriminability (D)** as a quantitative measure of how well molecular observables separate clinical outcomes.

The framework shows that:
- Observable variation does not necessarily translate into clinical relevance
- Only a small subset of observables aligns with outcome-relevant differences
- Observable-based inference is fundamentally constrained by projection from underlying biological states

---

## Key Concept

Discriminability is defined as:

D(O_k) = -log10(p)

where p is the log-rank p-value from survival analysis.

- High D → strong outcome separation  
- Low D → little or no outcome separation  

---

## Repository Contents

- `analysis/` : Core analysis scripts  
- `data_processing/` : Data preprocessing pipeline  
- `figures/` : Scripts for generating figures  
- `notebooks/` : Jupyter notebooks (if applicable)

---

## Data Source

Data used in this study are obtained from:
- TCGA via UCSC Xena (https://xenabrowser.net)

---

## Reproducibility

All analyses are based on:
- Gene expression (GE)
- Hallmark gene sets (MSigDB)
- Survival endpoints (OS, PFI)

Detailed procedures are described in the manuscript and Supplementary Methods.

---

## Citation

If you use this work, please cite:

(Your paper here after acceptance)

---

## Contact

For questions, please contact:
<your email>
