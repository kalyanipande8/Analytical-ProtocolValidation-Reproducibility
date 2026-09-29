# QC Anomaly Detection Report

> All figures and numbers below come directly from running `python run_all.py` on simulated QC data. Nothing is hardcoded.

## Overview

- **254 QC batch records** analyzed across **18 protein targets**.
- **12 batches flagged** as out-of-spec against target-specific, validation-derived acceptance criteria (mean +/- 3 SD control limits) - a 4.7% flag rate.

## Flag rate by target

| Target | QC batches | Flagged | Flag rate |
|---|---|---|---|
| TNF-alpha-mAb | 15 | 3 | 20.0% |
| Fc-Fusion-1 | 14 | 2 | 14.29% |
| Fc-Fusion-2 | 15 | 2 | 13.33% |
| PEG-Protein-1 | 14 | 1 | 7.14% |
| EGFR-mAb | 14 | 1 | 7.14% |
| PSMAxCD3-Bispecific | 14 | 1 | 7.14% |
| PD-1-mAb | 14 | 1 | 7.14% |
| HER2-mAb | 15 | 1 | 6.67% |
| ADC-Conjugate-1 | 14 | 0 | 0.0% |
| Recomb-Insulin-Analog | 13 | 0 | 0.0% |
| Recomb-EPO | 15 | 0 | 0.0% |
| Nanobody-VHH-1 | 13 | 0 | 0.0% |
| PD-L1-mAb | 14 | 0 | 0.0% |
| CD19-Bispecific | 14 | 0 | 0.0% |
| IL-6R-mAb | 13 | 0 | 0.0% |
| Fab-Fragment-1 | 14 | 0 | 0.0% |
| CD20-mAb | 15 | 0 | 0.0% |
| VEGF-mAb | 14 | 0 | 0.0% |

## Root-cause categorization of flagged batches

| Category | Count |
|---|---|
| purity_low_only | 6 |
| aggregation_high_only | 5 |
| purity_low_and_aggregation_high | 1 |

| Severity (max |z-score| vs. validation mean/SD) | Count |
|---|---|
| marginal (3<=|z|<4) | 6 |
| severe (|z|>=5) | 4 |
| major (4<=|z|<5) | 2 |

## Detector performance vs. injected ground truth

For this demo, a known subset of QC batches was deliberately drawn from a shifted (out-of-spec) distribution during simulation, tagged with a root cause (degraded reference standard, reagent lot drift, instrument fault, sample mishandling stress). This ground truth is **not** used by the detector - it's used here only to check the detector honestly, the way a method-validation exercise would check assay sensitivity against spiked samples.

- Injected out-of-spec batches: **9**
- Correctly flagged (true positives): **9**
- Missed (false negatives): **0**
- Flagged batches with no injected shift (false positives, expected at low frequency from the 3-sigma control limit): **3**
- Detector sensitivity: **100.0%**

## Figures

- `output/figures/flag_rate_by_target.png`
- `output/figures/control_chart_example.png`
- `output/figures/flag_category_breakdown.png`
