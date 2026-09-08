from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from qc_counts import filter_genes, log2_cpm, run_pca, validate_metadata


def example_counts():
    return pd.DataFrame(
        {
            "GSM1": [100, 30, 0],
            "GSM2": [120, 20, 1],
            "GSM3": [90, 40, 0],
            "GSM4": [110, 25, 2],
        },
        index=["ENSG1", "ENSG2", "ENSG3"],
    )


def test_filter_genes_and_log2cpm():
    counts = example_counts()
    filtered = filter_genes(counts, min_count=10, min_samples=3)
    assert filtered.index.tolist() == ["ENSG1", "ENSG2"]
    transformed = log2_cpm(filtered)
    assert transformed.shape == (2, 4)
    assert (transformed >= 0).all().all()


def test_pca_returns_two_components():
    transformed = log2_cpm(example_counts().loc[["ENSG1", "ENSG2"]])
    scores, variance = run_pca(transformed, top_variable_genes=2)
    assert scores.shape == (4, 2)
    assert len(variance) == 2


def test_validate_metadata_requires_every_sample():
    counts = example_counts()
    metadata = pd.DataFrame({"sample_id": ["GSM1", "GSM2", "GSM3"]})
    with pytest.raises(ValueError, match="metadata missing"):
        validate_metadata(counts, metadata)
