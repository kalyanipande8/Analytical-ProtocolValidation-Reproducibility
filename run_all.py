"""
run_all.py

Runs the full pipeline end-to-end:
  1. Simulate method validation + QC batch data for 18 protein targets -> data/raw/
  2. Derive acceptance criteria and standardize QC records against them -> data/processed/
  3. Anomaly detection, root-cause categorization, and figures -> output/

Run with:  python run_all.py
"""

import simulate_validation_data
import standardize
import anomaly_analysis


def main():
    print("--- Simulating method validation + QC batch data ---")
    simulate_validation_data.simulate()

    print("\n--- Deriving acceptance criteria and standardizing QC data ---")
    standardize.main()

    print("\n--- Anomaly detection analysis ---")
    anomaly_analysis.main()


if __name__ == "__main__":
    main()
