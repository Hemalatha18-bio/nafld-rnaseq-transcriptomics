#!/usr/bin/env python3
"""Generate sample-level QC summaries and exploratory PCA from raw RNA-seq counts.

Differential expression should use the untransformed integer counts with a count-aware
method such as DESeq2. The log2-CPM transform here is only for exploratory QC/PCA.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def load_counts(path: Path) -> pd.DataFrame:
    counts = pd.read_csv(path, sep="\t", index_col=0)
    if counts.empty:
        raise ValueError("count matrix is empty")
    numeric = counts.apply(pd.to_numeric, errors="raise")
    if numeric.isna().any().any():
        raise ValueError("count matrix contains missing values")
    if (numeric < 0).any().any():
        raise ValueError("count matrix contains negative values")
    if ((numeric % 1) != 0).any().any():
        raise ValueError("count matrix contains non-integer values")
    return numeric.astype("int64")


def validate_metadata(counts: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    if "sample_id" not in metadata.columns:
        raise ValueError("metadata must contain sample_id")
    if metadata["sample_id"].duplicated().any():
        raise ValueError("metadata contains duplicate sample_id values")
    metadata = metadata.set_index("sample_id")
    missing = [sample for sample in counts.columns if sample not in metadata.index]
    if missing:
        raise ValueError(f"metadata missing {len(missing)} count-matrix samples")
    return metadata.loc[counts.columns]


def filter_genes(counts: pd.DataFrame, min_count: int = 10, min_samples: int = 10) -> pd.DataFrame:
    if min_count < 0 or min_samples < 1:
        raise ValueError("min_count must be >= 0 and min_samples must be >= 1")
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    filtered = counts.loc[keep]
    if filtered.empty:
        raise ValueError("gene filter removed every gene")
    return filtered


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    library_sizes = counts.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("all samples must have positive library sizes")
    cpm = counts.divide(library_sizes, axis=1) * 1_000_000.0
    return np.log2(cpm + 1.0)


def run_pca(log_expression: pd.DataFrame, top_variable_genes: int = 500) -> tuple[pd.DataFrame, np.ndarray]:
    if top_variable_genes < 2:
        raise ValueError("top_variable_genes must be >= 2")
    variances = log_expression.var(axis=1).sort_values(ascending=False)
    selected = variances.head(min(top_variable_genes, len(variances))).index
    sample_matrix = log_expression.loc[selected].T
    if sample_matrix.shape[0] < 3 or sample_matrix.shape[1] < 2:
        raise ValueError("PCA requires at least 3 samples and 2 genes")
    model = PCA(n_components=2)
    scores = model.fit_transform(sample_matrix)
    frame = pd.DataFrame(scores, index=sample_matrix.index, columns=["PC1", "PC2"])
    frame.index.name = "sample_id"
    return frame, model.explained_variance_ratio_


def generate_qc(
    count_matrix: Path,
    metadata_path: Path,
    results_dir: Path,
    figures_dir: Path,
    min_count: int = 10,
    min_samples: int = 10,
    top_variable_genes: int = 500,
) -> None:
    counts = load_counts(count_matrix)
    metadata = pd.read_csv(metadata_path, sep="\t")
    metadata_aligned = validate_metadata(counts, metadata)
    filtered = filter_genes(counts, min_count=min_count, min_samples=min_samples)

    library_sizes = counts.sum(axis=0)
    detected_genes = (counts > 0).sum(axis=0)
    qc = pd.DataFrame(
        {
            "sample_id": counts.columns,
            "library_size": library_sizes.to_numpy(),
            "detected_genes": detected_genes.to_numpy(),
        }
    )
    if "analysis_group" in metadata_aligned.columns:
        qc["analysis_group"] = metadata_aligned["analysis_group"].to_numpy()

    log_expression = log2_cpm(filtered)
    pca_scores, variance = run_pca(log_expression, top_variable_genes=top_variable_genes)
    if "analysis_group" in metadata_aligned.columns:
        pca_scores["analysis_group"] = metadata_aligned.loc[pca_scores.index, "analysis_group"]

    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    qc.to_csv(results_dir / "sample_qc.tsv", sep="\t", index=False)
    pca_scores.reset_index().to_csv(results_dir / "pca_scores.tsv", sep="\t", index=False)

    summary = pd.DataFrame(
        {
            "metric": [
                "samples",
                "genes_raw",
                "genes_after_filter",
                "median_library_size",
                "pca_pc1_variance_fraction",
                "pca_pc2_variance_fraction",
            ],
            "value": [
                counts.shape[1],
                counts.shape[0],
                filtered.shape[0],
                float(library_sizes.median()),
                float(variance[0]),
                float(variance[1]),
            ],
        }
    )
    summary.to_csv(results_dir / "qc_summary.tsv", sep="\t", index=False)

    plt.figure(figsize=(8, 5))
    plt.hist(library_sizes.to_numpy(), bins=25)
    plt.xlabel("Library size (total counts)")
    plt.ylabel("Samples")
    plt.title("GSE135251 library-size distribution")
    plt.tight_layout()
    plt.savefig(figures_dir / "library_size_distribution.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 6))
    groups = pca_scores.get("analysis_group")
    if groups is None:
        plt.scatter(pca_scores["PC1"], pca_scores["PC2"], alpha=0.75)
    else:
        for group in sorted(groups.astype(str).unique()):
            subset = pca_scores.loc[groups.astype(str) == group]
            plt.scatter(subset["PC1"], subset["PC2"], label=group, alpha=0.75)
        plt.legend(frameon=False, fontsize=8)
    plt.xlabel(f"PC1 ({variance[0] * 100:.1f}% variance)")
    plt.ylabel(f"PC2 ({variance[1] * 100:.1f}% variance)")
    plt.title("Exploratory PCA of filtered log2-CPM counts")
    plt.tight_layout()
    plt.savefig(figures_dir / "pca_log2cpm.png", dpi=160)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count-matrix", type=Path, default=Path("data/processed/count_matrix.tsv"))
    parser.add_argument("--metadata", type=Path, default=Path("metadata/sample_metadata.tsv"))
    parser.add_argument("--results-dir", type=Path, default=Path("results/qc"))
    parser.add_argument("--figures-dir", type=Path, default=Path("figures/qc"))
    parser.add_argument("--min-count", type=int, default=10)
    parser.add_argument("--min-samples", type=int, default=10)
    parser.add_argument("--top-variable-genes", type=int, default=500)
    args = parser.parse_args()

    generate_qc(
        args.count_matrix,
        args.metadata,
        args.results_dir,
        args.figures_dir,
        min_count=args.min_count,
        min_samples=args.min_samples,
        top_variable_genes=args.top_variable_genes,
    )


if __name__ == "__main__":
    main()
