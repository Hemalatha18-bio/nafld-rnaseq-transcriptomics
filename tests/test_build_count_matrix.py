import gzip
from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from build_count_matrix import build_matrix, read_count_file, sample_id_from_filename


def write_counts(path, rows):
    with gzip.open(path, "wt") as handle:
        for gene, count in rows:
            handle.write(f"{gene}\t{count}\n")


def test_sample_id_from_filename():
    assert sample_id_from_filename(Path("GSM3998336_example.counts.txt.gz")) == "GSM3998336"
    with pytest.raises(ValueError):
        sample_id_from_filename(Path("sample.counts.txt.gz"))


def test_read_count_file_removes_htseq_summary_rows(tmp_path):
    path = tmp_path / "GSM0000001_demo.counts.txt.gz"
    write_counts(path, [("ENSG1", 10), ("ENSG2", 0), ("__no_feature", 5)])
    counts = read_count_file(path)
    assert counts.to_dict() == {"ENSG1": 10, "ENSG2": 0}
    assert counts.name == "GSM0000001"


def test_build_matrix_combines_samples(tmp_path):
    write_counts(
        tmp_path / "GSM0000001_a.counts.txt.gz",
        [("ENSG1", 10), ("ENSG2", 2), ("__no_feature", 1)],
    )
    write_counts(
        tmp_path / "GSM0000002_b.counts.txt.gz",
        [("ENSG1", 12), ("ENSG2", 4), ("__no_feature", 2)],
    )
    matrix = build_matrix(tmp_path)
    expected = pd.DataFrame(
        {"GSM0000001": [10, 2], "GSM0000002": [12, 4]},
        index=pd.Index(["ENSG1", "ENSG2"], name="gene_id"),
    )
    pd.testing.assert_frame_equal(matrix, expected)


def test_build_matrix_rejects_mismatched_gene_rows(tmp_path):
    write_counts(tmp_path / "GSM0000001_a.counts.txt.gz", [("ENSG1", 10), ("ENSG2", 2)])
    write_counts(tmp_path / "GSM0000002_b.counts.txt.gz", [("ENSG1", 12), ("ENSG3", 4)])
    with pytest.raises(ValueError, match="gene rows differ"):
        build_matrix(tmp_path)
