from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from audit_dataset import audit_dataset


def test_audit_dataset_reports_groups_and_cohort_counts():
    counts = pd.DataFrame(
        {"GSM1": [10, 5], "GSM2": [12, 6], "GSM3": [14, 7]},
        index=["ENSG1", "ENSG2"],
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["GSM1", "GSM2", "GSM3"],
            "analysis_group": ["Control", "NAFLD_F0_F1", "NAFLD_F3_F4"],
            "disease": ["Control", "NAFLD", "NAFLD"],
        }
    )
    audit, groups = audit_dataset(counts, metadata, expected_samples=3)
    metrics = dict(zip(audit["metric"], audit["value"]))
    assert metrics["samples"] == 3
    assert metrics["genes"] == 2
    assert metrics["control_samples"] == 1
    assert metrics["nafld_samples"] == 2
    assert groups["sample_count"].sum() == 3


def test_audit_dataset_rejects_sample_mismatch():
    counts = pd.DataFrame({"GSM1": [1], "GSM2": [2]}, index=["ENSG1"])
    metadata = pd.DataFrame(
        {
            "sample_id": ["GSM1", "GSM3"],
            "analysis_group": ["Control", "NAFLD_F0_F1"],
        }
    )
    with pytest.raises(ValueError, match="sample mismatch"):
        audit_dataset(counts, metadata, expected_samples=None)


def test_audit_dataset_checks_expected_complete_series_size():
    counts = pd.DataFrame({"GSM1": [1]}, index=["ENSG1"])
    metadata = pd.DataFrame(
        {"sample_id": ["GSM1"], "analysis_group": ["Control"], "disease": ["Control"]}
    )
    with pytest.raises(ValueError, match="expected 216 samples"):
        audit_dataset(counts, metadata, expected_samples=216)
