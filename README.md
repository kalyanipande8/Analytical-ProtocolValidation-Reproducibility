# Analytical Protocol Validation & Reproducibility

Method-validation and QC anomaly-detection pipeline: validates an analytical
method across **18 protein targets**, derives statistically-grounded
acceptance criteria from replicate precision data, and screens **250+ routine
QC batch records** against those criteria to flag out-of-spec batches with a
root-cause categorization.

> **Disclaimer:** This is a portfolio/demo project. All method-validation and
> QC data is **synthetically simulated** by `simulate_validation_data.py`. No
> real molecule, program, or employer data is used. Target names (mAbs,
> bispecifics, Fc-fusions, ADCs, etc.) and their purity/aggregation profiles
> are simulated from published, general biotherapeutic-analytics knowledge,
> not from any real dataset.

## The 18 protein targets

18 targets spanning common therapeutic-protein modalities (8 mAbs, 2
bispecifics, 2 Fc-fusions, 1 nanobody, 1 Fab fragment, 1 ADC, 1 PEGylated
protein, 1 insulin analog, 1 EPO), each with its own process-typical
%purity and %aggregation baseline and assay variability.

## What it does

1. **`simulate_validation_data.py`** simulates two linked datasets:
   - A **method validation study**: 20 replicate runs per target (3 analysts
     x 5 days each), the repeatability/intermediate-precision design a real
     ICH Q2(R2) validation would use, measuring %purity and %aggregation.
   - **250+ routine QC batch records** distributed across the 18 targets,
     drawn from the same process distribution, with a small deliberately
     injected subset (~4-5%) shifted out of range via one of four real
     root-cause profiles (degraded reference standard, reagent lot drift,
     instrument fault, sample mishandling stress). The injected ground truth
     is kept only to sanity-check the detector - it is not used by the
     detection logic itself.

2. **`standardize.py`** derives **acceptance criteria per target** directly
   from the validation data: lower control limit = validation mean - 3 SD
   for purity, upper control limit = validation mean + 3 SD for aggregation
   (a standard statistical-process-control convention). Every QC batch is
   then scored (z-score) and flagged against its own target's criteria - no
   single global spec, because each target has its own validated process
   range.

3. **`anomaly_analysis.py`** runs the anomaly-detection and QC triage:
   - Flag-rate summary per target
   - Root-cause categorization of every flagged batch (purity-low-only,
     aggregation-high-only, or both) and a severity bucket by |z-score|
   - An honest detector-performance check against the simulation's injected
     ground truth (sensitivity, false positives) - included as a method
     check, not something a real QC group would have
   - Control-chart-style figures

4. **`run_all.py`** runs the whole pipeline.

## How to run

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python run_all.py
pytest -q
```

## Results

_All numbers below are copied directly from an actual run of `python
run_all.py` - nothing is hardcoded. Fixed random seed makes them
reproducible._

- **360 method-validation records** (18 targets x 20 replicate runs) used to
  derive acceptance criteria.
- **254 QC batch records** analyzed across all 18 targets.
- **12 of 254 batches flagged (4.7%)** against target-specific, validation-
  derived acceptance criteria (mean +/- 3 SD control limits).
  - 6 purity-low-only, 5 aggregation-high-only, 1 both.
  - Severity: mostly marginal-to-major (3 <= |z| < 5), with the mishandling
    and instrument-fault cases (which shift both attributes at once) landing
    in the more severe bucket.
- **Detector performance vs. injected ground truth**: 9 batches were
  deliberately drawn from a shifted distribution during simulation; **9/9
  (100%) were caught**, plus 3 additional flags from the 3-sigma control
  limit's expected background false-positive rate on 254 x 2 checks. This is
  an honest sensitivity check the pipeline runs on itself - a real QC group
  would never have this ground truth, only the flag list.
- `pytest -q`: **9 passed**.

### Example acceptance criteria (mean +/- 3 SD, from `data/processed/acceptance_criteria.csv`)

| Target | Purity mean (SD) | Purity lower limit | Aggregation mean (SD) | Aggregation upper limit |
|---|---|---|---|---|
| TNF-alpha-mAb | 98.47 (0.25) | 97.72 | 1.07 (0.15) | 1.52 |
| ADC-Conjugate-1 | 95.45 (0.43) | 94.17 | 2.73 (0.46) | 4.11 |
| Nanobody-VHH-1 | 99.06 (0.24) | 98.33 | 0.52 (0.07) | 0.74 |
| PSMAxCD3-Bispecific | 95.80 (0.47) | 94.38 | 2.52 (0.41) | 3.77 |

(Full table in `data/processed/acceptance_criteria.csv` for all 18 targets -
this is a snapshot; re-running reproduces it exactly.) Note each target gets
its own criteria derived from its own validation data, rather than one
blanket spec - a bispecific and a nanobody have genuinely different validated
process ranges.

## Figures

- `output/figures/flag_rate_by_target.png` - QC flag rate by protein target
- `output/figures/control_chart_example.png` - control chart (purity +
  aggregation vs. acceptance limits) for the target with the most flags
- `output/figures/flag_category_breakdown.png` - flagged batches by root-
  cause category

## Project structure

```
simulate_validation_data.py   # 18 targets: validation study + 250+ QC batches -> data/raw/*.csv
standardize.py                 # acceptance criteria (mean +/- 3 SD) + QC flags -> data/processed/
anomaly_analysis.py            # flag summary, root-cause categorization, figures -> output/
run_all.py                     # orchestrates the full pipeline
tests/test_pipeline.py         # pytest checks incl. injected out-of-spec/in-spec sanity checks
```

## Tech stack

- Python 3
- pandas / numpy for data simulation and wrangling
- matplotlib for figures
- pytest
