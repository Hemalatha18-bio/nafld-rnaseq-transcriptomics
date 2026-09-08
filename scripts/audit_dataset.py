#!/usr/bin/env python3
"""Audit GSE135251 count/metadata consistency before downstream analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def audit_dataset(
    count_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    expected_samples: int | None = 216,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if count_matrix.empty:
        raise ValueError("count matrix is empty")
    if "sample_id" not in metadata.columns:
        raise ValueError("metadata must contain sample_id")
    if metadata["sample_id"].duplicated().any():
        raise ValueError("metadata contains duplicate sample_id values")
    if count_matrix.columns.duplicated().any():
        raise ValueError("count matrix contains duplicate sample columns")

    count_samples = set(map(str, count_matrix.columns))
    metadata_samples = set(metadata["sample_id"].astype(str))
    only_counts = sorted(count_samples - metadata_samples)
    only_metadata = sorted(metadata_samples - count_samples)
    if only_counts or only_metadata:
        raise ValueError(
            "count/metadata sample mismatch: "
            f"only_counts={len(only_counts)}, only_metadata={len(only_metadata)}"
        )

    if expected_samples is not None and len(count_samples) != expected_samples:
        raise ValueError(
            f"expected {expected_samples} samples for the complete GEO series, "
            f"found {len(count_samples)}"
        )

    if "analysis_group" not in metadata.columns:
        raise ValueError("metadata must contain analysis_group")

    aligned = metadata.set_index("sample_id").loc[count_matrix.columns]
    group_counts = (
        aligned["analysis_group"]
        .fillna("Missing")
        .astype(str)
        .value_counts()
        .rename_axis("analysis_group")
        .reset_index(name="sample_count")
    )

    disease_counts = (
        aligned.get("disease", pd.Series(index=aligned.index, dtype="object"))
        .fillna("Missing")
        .astype(str)
        .value_counts()
    )

    audit = pd.DataFrame(
        {
            "metric": [
                "samples",
                "genes",
                "metadata_rows",
                "unique_analysis_groups",
                "control_samples",
                "nafld_samples",
            ],
            "value": [
                count_matrix.shape[1],
                count_matrix.shape[0],
                len(metadata),
                group_counts.shape[0],
                int(disease_counts.get("Control", 0)),
                int(disease_counts.get("NAFLD", 0)),
            ],
        }
    )
    return audit, group_counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count-matrix", type=Path, default=Path("data/processed/count_matrix.tsv"))
    parser.add_argument("--metadata", type=Path, default=Path("metadata/sample_metadata.tsv"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/audit"))
    parser.add_argument("--expected-samples", type=int, default=216)
    args = parser.parse_args()

    counts = pd.read_csv(args.count_matrix, sep="\t", index_col=0)
    metadata = pd.read_csv(args.metadata, sep="\t")
    audit, groups = audit_dataset(counts, metadata, expected_samples=args.expected_samples)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output_dir / "dataset_audit.tsv", sep="\t", index=False)
    groups.to_csv(args.output_dir / "analysis_group_counts.tsv", sep="\t", index=False)
    print(audit.to_string(index=False))
    print(groups.to_string(index=False))


if __name__ == "__main__":
    main()
