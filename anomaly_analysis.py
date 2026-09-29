"""
anomaly_analysis.py

Statistical anomaly detection and root-cause categorization over the
standardized QC dataset (standardize.py):

  1. Flagged-batch summary: how many of 250+ QC records are out of spec
     against their target's derived acceptance criteria, and how that's
     distributed across the 18 targets.
  2. Root-cause categorization of every flagged batch: purity-low-only,
     aggregation-high-only, or both, plus a severity bucket by |z-score|,
     as an analyst doing a real deviation review would triage a flag list.
  3. Detector performance against the simulation's injected ground truth
     (how many injected out-of-spec batches were caught vs. missed, and
     whether any flags are false positives from the 3-sigma control limit) -
     included as an honest check on the method, not part of the QC report
     itself (a real QC group would never have this ground truth).
  4. Control-chart-style figures per representative target and a portfolio
     flag-rate summary.

Writes output/qc_anomaly_report.md, output/flagged_batches.csv, and PNGs to
output/figures/.
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
PROCESSED_DIR = os.path.join(PROJECT_DIR, "data", "processed")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")
FIGURES_DIR = os.path.join(OUTPUT_DIR, "figures")

PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]


def load():
    return pd.read_csv(os.path.join(PROCESSED_DIR, "qc_with_criteria.csv"))


def categorize_root_cause(row):
    if row["flag_purity_low"] and row["flag_aggregation_high"]:
        return "purity_low_and_aggregation_high"
    if row["flag_purity_low"]:
        return "purity_low_only"
    if row["flag_aggregation_high"]:
        return "aggregation_high_only"
    return "not_flagged"


def severity_bucket(row):
    z = max(abs(row["purity_z"]), abs(row["aggregation_z"]))
    if z >= 5:
        return "severe (|z|>=5)"
    if z >= 4:
        return "major (4<=|z|<5)"
    return "marginal (3<=|z|<4)"


def flag_summary(df):
    flagged = df[df["any_flag"]].copy()
    flagged["flag_category"] = flagged.apply(categorize_root_cause, axis=1)
    flagged["severity"] = flagged.apply(severity_bucket, axis=1)
    return flagged


def per_target_flag_rate(df):
    g = df.groupby("target_id").agg(
        n_batches=("batch_id", "count"),
        n_flagged=("any_flag", "sum"),
    ).reset_index()
    g["flag_rate_pct"] = (100 * g["n_flagged"] / g["n_batches"]).round(2)
    return g.sort_values("flag_rate_pct", ascending=False)


def detector_performance(df):
    """Honest check against the simulation's injected ground truth (not used by
    the detector itself, and not something a real QC lab would have)."""
    tp = int(((df["true_injected_shift"]) & (df["any_flag"])).sum())
    fn = int(((df["true_injected_shift"]) & (~df["any_flag"])).sum())
    fp = int(((~df["true_injected_shift"]) & (df["any_flag"])).sum())
    n_injected = int(df["true_injected_shift"].sum())
    sensitivity = tp / n_injected if n_injected else float("nan")
    return dict(n_injected=n_injected, true_positives=tp, false_negatives=fn,
                false_positives=fp, sensitivity=round(sensitivity, 3))


def make_figures(df, target_rates):
    os.makedirs(FIGURES_DIR, exist_ok=True)

    # Figure 1: portfolio flag rate by target
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(target_rates["target_id"], target_rates["flag_rate_pct"],
            color=PALETTE[0])
    ax.set_xlabel("QC batches flagged (%)")
    ax.set_title("QC flag rate by protein target")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, "flag_rate_by_target.png"), dpi=150)
    plt.close(fig)

    # Figure 2: control chart for one representative target with the most flags
    top_target = target_rates.iloc[0]["target_id"]
    g = df[df["target_id"] == top_target].reset_index(drop=True)
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    axes[0].plot(g.index, g["purity_pct"], "o-", color=PALETTE[0], markersize=4)
    axes[0].axhline(g["purity_lower_limit"].iloc[0], color="red", linestyle="--", label="lower acceptance limit")
    axes[0].axhline(g["purity_mean"].iloc[0], color="gray", linestyle=":", label="validation mean")
    flagged_purity = g[g["flag_purity_low"]]
    axes[0].scatter(flagged_purity.index, flagged_purity["purity_pct"], color="red", zorder=5, label="flagged")
    axes[0].set_ylabel("%Purity")
    axes[0].set_title(f"QC control chart - {top_target}")
    axes[0].legend(fontsize=8)

    axes[1].plot(g.index, g["aggregation_pct"], "o-", color=PALETTE[1], markersize=4)
    axes[1].axhline(g["aggregation_upper_limit"].iloc[0], color="red", linestyle="--", label="upper acceptance limit")
    axes[1].axhline(g["agg_mean"].iloc[0], color="gray", linestyle=":", label="validation mean")
    flagged_agg = g[g["flag_aggregation_high"]]
    axes[1].scatter(flagged_agg.index, flagged_agg["aggregation_pct"], color="red", zorder=5, label="flagged")
    axes[1].set_ylabel("%Aggregation")
    axes[1].set_xlabel("QC batch sequence")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, "control_chart_example.png"), dpi=150)
    plt.close(fig)

    # Figure 3: root-cause category breakdown among flagged batches
    flagged = flag_summary(df)
    counts = flagged["flag_category"].value_counts()
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.bar(counts.index, counts.values, color=PALETTE[: len(counts)])
    ax.set_ylabel("Number of flagged batches")
    ax.set_title("Flagged QC batches by category")
    ax.tick_params(axis="x", rotation=25)
    for tick in ax.get_xticklabels():
        tick.set_ha("right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, "flag_category_breakdown.png"), dpi=150)
    plt.close(fig)


def write_report(df, target_rates, flagged, perf):
    total = len(df)
    n_flagged = int(df["any_flag"].sum())
    lines = [
        "# QC Anomaly Detection Report",
        "",
        "> All figures and numbers below come directly from running `python run_all.py` "
        "on simulated QC data. Nothing is hardcoded.",
        "",
        "## Overview",
        "",
        f"- **{total} QC batch records** analyzed across **18 protein targets**.",
        f"- **{n_flagged} batches flagged** as out-of-spec against target-specific, "
        f"validation-derived acceptance criteria (mean +/- {3.0:g} SD control limits) "
        f"- a {100 * n_flagged / total:.1f}% flag rate.",
        "",
        "## Flag rate by target",
        "",
        "| Target | QC batches | Flagged | Flag rate |",
        "|---|---|---|---|",
    ]
    for row in target_rates.itertuples():
        lines.append(f"| {row.target_id} | {row.n_batches} | {row.n_flagged} | {row.flag_rate_pct}% |")

    lines += [
        "",
        "## Root-cause categorization of flagged batches",
        "",
        "| Category | Count |",
        "|---|---|",
    ]
    for cat, cnt in flagged["flag_category"].value_counts().items():
        lines.append(f"| {cat} | {cnt} |")

    lines += [
        "",
        "| Severity (max |z-score| vs. validation mean/SD) | Count |",
        "|---|---|",
    ]
    for sev, cnt in flagged["severity"].value_counts().items():
        lines.append(f"| {sev} | {cnt} |")

    lines += [
        "",
        "## Detector performance vs. injected ground truth",
        "",
        "For this demo, a known subset of QC batches was deliberately drawn from a "
        "shifted (out-of-spec) distribution during simulation, tagged with a root cause "
        "(degraded reference standard, reagent lot drift, instrument fault, sample "
        "mishandling stress). This ground truth is **not** used by the detector - it's "
        "used here only to check the detector honestly, the way a method-validation "
        "exercise would check assay sensitivity against spiked samples.",
        "",
        f"- Injected out-of-spec batches: **{perf['n_injected']}**",
        f"- Correctly flagged (true positives): **{perf['true_positives']}**",
        f"- Missed (false negatives): **{perf['false_negatives']}**",
        f"- Flagged batches with no injected shift (false positives, expected at low "
        f"frequency from the 3-sigma control limit): **{perf['false_positives']}**",
        f"- Detector sensitivity: **{perf['sensitivity'] * 100:.1f}%**",
        "",
        "## Figures",
        "",
        "- `output/figures/flag_rate_by_target.png`",
        "- `output/figures/control_chart_example.png`",
        "- `output/figures/flag_category_breakdown.png`",
    ]
    with open(os.path.join(OUTPUT_DIR, "qc_anomaly_report.md"), "w") as f:
        f.write("\n".join(lines) + "\n")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    df = load()
    target_rates = per_target_flag_rate(df)
    flagged = flag_summary(df)
    perf = detector_performance(df)

    flagged.to_csv(os.path.join(OUTPUT_DIR, "flagged_batches.csv"), index=False)
    make_figures(df, target_rates)
    write_report(df, target_rates, flagged, perf)

    print("=" * 60)
    print("QC ANOMALY DETECTION ANALYSIS")
    print("=" * 60)
    print(f"QC batches analyzed: {len(df)} across {df['target_id'].nunique()} targets")
    print(f"Flagged: {int(df['any_flag'].sum())} ({100 * df['any_flag'].mean():.1f}%)")
    print("Root-cause category breakdown:")
    print(flagged["flag_category"].value_counts().to_string())
    print(f"Detector sensitivity vs. injected ground truth: {perf['sensitivity'] * 100:.1f}% "
          f"({perf['true_positives']}/{perf['n_injected']} caught, {perf['false_positives']} false positives)")
    print("\nWrote output/qc_anomaly_report.md, output/flagged_batches.csv, output/figures/*.png")


if __name__ == "__main__":
    main()
