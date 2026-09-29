"""
simulate_validation_data.py

Simulates two linked datasets for an analytical method validation program
across 18 distinct protein targets:

  1. Method validation study data (data/raw/method_validation.csv)
     For each target, a repeatability + intermediate-precision study: N_VALIDATION_RUNS
     independent replicate measurements of %purity (e.g. RP-HPLC / CE-SDS main peak)
     and %aggregation (e.g. SEC %HMW), spread across multiple analysts and days, the
     way a real method validation (ICH Q2(R2)) precision study is actually run. This
     is the data acceptance criteria get DERIVED from - it is not QC data itself.

  2. Routine QC batch/trace records (data/raw/qc_batches.csv)
     250+ QC batch release records distributed across the 18 targets, drawn from the
     same underlying process distribution as the validation study, with a small,
     deliberately injected subset shifted out of the normal process range (drifted
     reagent, degraded reference standard, instrument fault, sample mishandling) to
     represent real out-of-spec / atypical batches a QC lab would need to catch.
     `true_injected_shift` is the simulation's ground truth (not used by the
     detection logic) - kept only so the test suite and root-cause writeup can check
     the detector actually catches what was planted.

Uses a fixed random seed for reproducibility.
"""

import os

import numpy as np
import pandas as pd

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(PROJECT_DIR, "data", "raw")

RANDOM_SEED = 42
N_VALIDATION_RUNS = 20   # replicate runs per target in the precision study
N_ANALYSTS = 3
N_DAYS = 5
N_QC_BATCHES_TARGET_TOTAL = 252  # ~14 per target, distributed unevenly below

# 18 protein targets spanning common therapeutic-protein modalities, each with its
# own process-typical purity / aggregation baseline and measurement variability
# (literature-plausible, not real product data).
TARGETS = {
    "TNF-alpha-mAb":     dict(purity_mean=98.4, purity_sd=0.30, agg_mean=1.1, agg_sd=0.18),
    "PD-1-mAb":          dict(purity_mean=98.1, purity_sd=0.28, agg_mean=1.3, agg_sd=0.20),
    "PD-L1-mAb":         dict(purity_mean=97.9, purity_sd=0.32, agg_mean=1.4, agg_sd=0.22),
    "HER2-mAb":          dict(purity_mean=98.6, purity_sd=0.25, agg_mean=1.0, agg_sd=0.15),
    "EGFR-mAb":          dict(purity_mean=98.0, purity_sd=0.30, agg_mean=1.2, agg_sd=0.19),
    "VEGF-mAb":          dict(purity_mean=98.3, purity_sd=0.27, agg_mean=1.1, agg_sd=0.17),
    "CD20-mAb":          dict(purity_mean=98.5, purity_sd=0.26, agg_mean=1.0, agg_sd=0.16),
    "IL-6R-mAb":         dict(purity_mean=97.8, purity_sd=0.34, agg_mean=1.5, agg_sd=0.24),
    "CD19-Bispecific":   dict(purity_mean=96.2, purity_sd=0.45, agg_mean=2.3, agg_sd=0.35),
    "PSMAxCD3-Bispecific": dict(purity_mean=95.8, purity_sd=0.48, agg_mean=2.5, agg_sd=0.38),
    "Fc-Fusion-1":       dict(purity_mean=97.2, purity_sd=0.38, agg_mean=1.6, agg_sd=0.28),
    "Fc-Fusion-2":       dict(purity_mean=97.0, purity_sd=0.40, agg_mean=1.7, agg_sd=0.30),
    "Nanobody-VHH-1":    dict(purity_mean=99.1, purity_sd=0.18, agg_mean=0.5, agg_sd=0.10),
    "Fab-Fragment-1":    dict(purity_mean=98.2, purity_sd=0.29, agg_mean=0.7, agg_sd=0.13),
    "ADC-Conjugate-1":   dict(purity_mean=95.5, purity_sd=0.42, agg_mean=2.8, agg_sd=0.40),
    "PEG-Protein-1":     dict(purity_mean=97.6, purity_sd=0.33, agg_mean=1.2, agg_sd=0.21),
    "Recomb-Insulin-Analog": dict(purity_mean=99.0, purity_sd=0.20, agg_mean=0.4, agg_sd=0.09),
    "Recomb-EPO":        dict(purity_mean=98.7, purity_sd=0.24, agg_mean=0.8, agg_sd=0.14),
}

assert len(TARGETS) == 18

ANALYSTS = ["Analyst-1", "Analyst-2", "Analyst-3"]

# Root-cause categories a QC deviation investigation would actually assign,
# and the process bias each imprints on purity/aggregation when injected.
OUT_OF_SPEC_CAUSES = {
    "degraded_reference_standard": dict(purity_shift=-1.6, agg_shift=+0.0),
    "reagent_lot_drift":           dict(purity_shift=-0.9, agg_shift=+0.6),
    "instrument_fault":            dict(purity_shift=-0.3, agg_shift=+1.8),
    "sample_mishandling_stress":   dict(purity_shift=-1.1, agg_shift=+1.4),
}


def simulate_validation_study(rng):
    rows = []
    for target_id, p in TARGETS.items():
        for run in range(1, N_VALIDATION_RUNS + 1):
            analyst = ANALYSTS[(run - 1) % N_ANALYSTS]
            day = ((run - 1) % N_DAYS) + 1
            # small analyst-to-analyst offset representative of intermediate precision
            analyst_offset = {"Analyst-1": 0.0, "Analyst-2": 0.05, "Analyst-3": -0.05}[analyst]
            purity = np.clip(rng.normal(p["purity_mean"] + analyst_offset, p["purity_sd"]), 0, 100)
            agg = np.clip(rng.normal(p["agg_mean"] - analyst_offset, p["agg_sd"]), 0, None)
            rows.append({
                "target_id": target_id, "run_id": run, "analyst": analyst, "day": day,
                "purity_pct": round(purity, 3), "aggregation_pct": round(agg, 3),
            })
    return pd.DataFrame(rows)


def simulate_qc_batches(rng):
    """250+ routine QC batch records across 18 targets, ~4.5% deliberately shifted
    out of the validated process range using one of four real root-cause profiles."""
    targets = list(TARGETS.keys())
    # uneven but reasonable batch counts per target (12-17), summing to ~252
    base_counts = [15, 14, 14, 15, 14, 14, 15, 13, 14, 14, 14, 15, 13, 14, 14, 14, 13, 15]
    assert len(base_counts) == 18

    rows = []
    batch_counter = 1
    for target_id, n_batches in zip(targets, base_counts):
        p = TARGETS[target_id]
        for _ in range(n_batches):
            batch_id = f"QC-{batch_counter:04d}"
            batch_counter += 1
            analyst = rng.choice(ANALYSTS)
            day = int(rng.integers(1, 61))  # QC program spans ~2 months

            # ~4.5% of batches deliberately injected as out-of-spec via a real cause
            is_injected = rng.random() < 0.045
            cause = None
            if is_injected:
                cause = rng.choice(list(OUT_OF_SPEC_CAUSES.keys()))
                shift = OUT_OF_SPEC_CAUSES[cause]
                purity = np.clip(rng.normal(p["purity_mean"] + shift["purity_shift"], p["purity_sd"]), 0, 100)
                agg = np.clip(rng.normal(p["agg_mean"] + shift["agg_shift"], p["agg_sd"]), 0, None)
            else:
                purity = np.clip(rng.normal(p["purity_mean"], p["purity_sd"]), 0, 100)
                agg = np.clip(rng.normal(p["agg_mean"], p["agg_sd"]), 0, None)

            rows.append({
                "batch_id": batch_id, "target_id": target_id, "analyst": analyst,
                "day": day, "purity_pct": round(purity, 3), "aggregation_pct": round(agg, 3),
                "true_injected_shift": is_injected, "true_root_cause": cause or "",
            })
    return pd.DataFrame(rows)


def simulate():
    rng = np.random.default_rng(RANDOM_SEED)
    os.makedirs(RAW_DIR, exist_ok=True)

    registry = pd.DataFrame({"target_id": list(TARGETS.keys())})
    registry.to_csv(os.path.join(RAW_DIR, "target_registry.csv"), index=False)

    validation = simulate_validation_study(rng)
    validation.to_csv(os.path.join(RAW_DIR, "method_validation.csv"), index=False)

    qc = simulate_qc_batches(rng)
    qc.to_csv(os.path.join(RAW_DIR, "qc_batches.csv"), index=False)

    print(f"Simulated method validation data: {len(TARGETS)} targets x {N_VALIDATION_RUNS} runs "
          f"= {len(validation)} validation records")
    print(f"Simulated QC batch data: {len(qc)} batch records across {len(TARGETS)} targets "
          f"({int(qc['true_injected_shift'].sum())} deliberately injected out-of-spec)")
    print(f"  target_registry.csv     {len(registry)} rows")
    print(f"  method_validation.csv   {len(validation)} rows")
    print(f"  qc_batches.csv          {len(qc)} rows")


if __name__ == "__main__":
    simulate()
