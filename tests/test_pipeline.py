import os
import sys

import numpy as np
import pandas as pd
import pytest

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_DIR)

import simulate_validation_data as sim  # noqa: E402
import standardize  # noqa: E402
import anomaly_analysis as aa  # noqa: E402


@pytest.fixture(scope="module", autouse=True)
def ensure_pipeline_ran():
    sim.simulate()
    standardize.standardize()
    yield


def test_18_targets_present():
    registry = pd.read_csv(os.path.join(sim.RAW_DIR, "target_registry.csv"))
    assert len(sim.TARGETS) == 18
    assert registry["target_id"].nunique() == 18


def test_validation_study_has_enough_replicates_per_target():
    validation = pd.read_csv(os.path.join(sim.RAW_DIR, "method_validation.csv"))
    counts = validation.groupby("target_id").size()
    assert len(counts) == 18
    assert (counts >= sim.N_VALIDATION_RUNS).all()


def test_at_least_250_qc_traces():
    qc = pd.read_csv(os.path.join(sim.RAW_DIR, "qc_batches.csv"))
    assert len(qc) >= 250
    assert qc["target_id"].nunique() == 18
    assert qc["batch_id"].is_unique


def test_acceptance_criteria_derived_for_every_target():
    criteria = pd.read_csv(os.path.join(standardize.PROCESSED_DIR, "acceptance_criteria.csv"))
    assert len(criteria) == 18
    assert (criteria["purity_lower_limit"] < criteria["purity_mean"]).all()
    assert (criteria["aggregation_upper_limit"] > criteria["agg_mean"]).all()


def test_flagged_count_in_reasonable_band():
    """Target is ~12 flagged batches; allow a reasonable statistical band."""
    qc = pd.read_csv(os.path.join(standardize.PROCESSED_DIR, "qc_with_criteria.csv"))
    n_flagged = int(qc["any_flag"].sum())
    assert 6 <= n_flagged <= 20, f"expected flagged count in [6, 20], got {n_flagged}"


def test_obviously_out_of_spec_record_is_flagged():
    """Injected sanity check: a batch far outside any target's validated range must
    be flagged by the same acceptance-criteria logic used in standardize.py."""
    criteria = pd.read_csv(os.path.join(standardize.PROCESSED_DIR, "acceptance_criteria.csv")).set_index("target_id")
    target_id = criteria.index[0]
    row = criteria.loc[target_id]
    bad_batch = pd.DataFrame([{
        "batch_id": "TEST-BAD", "target_id": target_id, "analyst": "Analyst-1", "day": 1,
        "purity_pct": row["purity_mean"] - 10 * row["purity_sd"],  # wildly low purity
        "aggregation_pct": row["agg_mean"] + 10 * row["agg_sd"],   # wildly high aggregation
        "true_injected_shift": True, "true_root_cause": "test_injected",
    }])
    checked = standardize.apply_criteria_to_qc(bad_batch, criteria.reset_index())
    assert bool(checked.loc[0, "any_flag"]) is True
    assert bool(checked.loc[0, "flag_purity_low"]) is True
    assert bool(checked.loc[0, "flag_aggregation_high"]) is True


def test_clearly_in_spec_record_is_not_flagged():
    """Injected sanity check: a batch sitting exactly at the validation mean for its
    target must NOT be flagged."""
    criteria = pd.read_csv(os.path.join(standardize.PROCESSED_DIR, "acceptance_criteria.csv")).set_index("target_id")
    target_id = criteria.index[0]
    row = criteria.loc[target_id]
    good_batch = pd.DataFrame([{
        "batch_id": "TEST-GOOD", "target_id": target_id, "analyst": "Analyst-1", "day": 1,
        "purity_pct": row["purity_mean"],
        "aggregation_pct": row["agg_mean"],
        "true_injected_shift": False, "true_root_cause": "",
    }])
    checked = standardize.apply_criteria_to_qc(good_batch, criteria.reset_index())
    assert bool(checked.loc[0, "any_flag"]) is False


def test_root_cause_categories_cover_all_flagged_batches():
    qc = pd.read_csv(os.path.join(standardize.PROCESSED_DIR, "qc_with_criteria.csv"))
    flagged = aa.flag_summary(qc)
    assert (flagged["flag_category"] != "not_flagged").all()
    assert set(flagged["flag_category"].unique()) <= {
        "purity_low_only", "aggregation_high_only", "purity_low_and_aggregation_high",
    }


def test_detector_catches_most_injected_out_of_spec_batches():
    qc = pd.read_csv(os.path.join(standardize.PROCESSED_DIR, "qc_with_criteria.csv"))
    perf = aa.detector_performance(qc)
    assert perf["n_injected"] > 0
    assert perf["sensitivity"] >= 0.5
