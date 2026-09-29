"""
standardize.py

Derives real, statistically-grounded acceptance criteria per target from the
method validation study, then standardizes every routine QC batch record
against those criteria - the "validated method -> acceptance criteria -> QC
disposition" chain a real analytical/QC group runs.

Acceptance criteria (per target, per attribute):
    lower/upper control limit = validation mean +/- K_SIGMA * validation SD
K_SIGMA = 3 (a standard control-limit convention: ICH Q2(R2)-style precision
data defines the expected process spread, and +/-3 SD is the conventional
statistical-process-control action limit for flagging a batch as atypical).

Purity has only a LOWER limit (low purity is the failure mode that matters);
aggregation has only an UPPER limit (high aggregation is the failure mode).

Writes:
  data/processed/acceptance_criteria.csv   one row per target
  data/processed/qc_with_criteria.csv      every QC batch + its z-scores + flags
"""

import os

import numpy as np
import pandas as pd

from simulate_validation_data import RAW_DIR

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
PROCESSED_DIR = os.path.join(PROJECT_DIR, "data", "processed")

K_SIGMA = 3.0


def derive_acceptance_criteria(validation: pd.DataFrame) -> pd.DataFrame:
    g = validation.groupby("target_id").agg(
        purity_mean=("purity_pct", "mean"),
        purity_sd=("purity_pct", "std"),
        agg_mean=("aggregation_pct", "mean"),
        agg_sd=("aggregation_pct", "std"),
        n_validation_runs=("run_id", "count"),
    ).reset_index()
    g["purity_lower_limit"] = g["purity_mean"] - K_SIGMA * g["purity_sd"]
    g["aggregation_upper_limit"] = g["agg_mean"] + K_SIGMA * g["agg_sd"]
    return g.round(4)


def apply_criteria_to_qc(qc: pd.DataFrame, criteria: pd.DataFrame) -> pd.DataFrame:
    merged = qc.merge(criteria, on="target_id", how="left")

    merged["purity_z"] = (merged["purity_pct"] - merged["purity_mean"]) / merged["purity_sd"]
    merged["aggregation_z"] = (merged["aggregation_pct"] - merged["agg_mean"]) / merged["agg_sd"]

    merged["flag_purity_low"] = merged["purity_pct"] < merged["purity_lower_limit"]
    merged["flag_aggregation_high"] = merged["aggregation_pct"] > merged["aggregation_upper_limit"]
    merged["any_flag"] = merged["flag_purity_low"] | merged["flag_aggregation_high"]

    return merged.round(4)


def standardize():
    validation = pd.read_csv(os.path.join(RAW_DIR, "method_validation.csv"))
    qc = pd.read_csv(os.path.join(RAW_DIR, "qc_batches.csv"))

    criteria = derive_acceptance_criteria(validation)
    standardized_qc = apply_criteria_to_qc(qc, criteria)

    os.makedirs(PROCESSED_DIR, exist_ok=True)
    criteria.to_csv(os.path.join(PROCESSED_DIR, "acceptance_criteria.csv"), index=False)
    standardized_qc.to_csv(os.path.join(PROCESSED_DIR, "qc_with_criteria.csv"), index=False)

    return criteria, standardized_qc


def main():
    criteria, standardized_qc = standardize()
    print(f"Derived acceptance criteria for {len(criteria)} targets from "
          f"{criteria['n_validation_runs'].sum()} validation runs")
    print(f"Standardized {len(standardized_qc)} QC batch records against target-specific criteria")
    print(f"  Flagged: {int(standardized_qc['any_flag'].sum())} of {len(standardized_qc)}")
    print("  Wrote data/processed/acceptance_criteria.csv, data/processed/qc_with_criteria.csv")


if __name__ == "__main__":
    main()
