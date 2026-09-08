#!/usr/bin/env python3
"""Combine GEO HTSeq-count files into a validated gene × sample matrix."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


SAMPLE_RE = re.compile(r"^(GSM\d+)_")


def sample_id_from_filename(path: Path) -> str:
    match = SAMPLE_RE.match(path.name)
    if not match:
        raise ValueError(f"cannot determine GEO sample accession from {path.name}")
    return match.group(1)


def read_count_file(path: Path) -> pd.Series:
    table = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["gene_id", "count"],
        usecols=[0, 1],
        compression="infer",
    )
    table["gene_id"] = table["gene_id"].astype(str)
    table = table.loc[~table["gene_id"].str.startswith("__")].copy()
    if table["gene_id"].duplicated().any():
        duplicated = table.loc[table["gene_id"].duplicated(), "gene_id"].iloc[0]
        raise ValueError(f"duplicate gene id {duplicated!r} in {path.name}")

    counts = pd.to_numeric(table["count"], errors="raise")
    if (counts < 0).any():
        raise ValueError(f"negative counts found in {path.name}")
    if ((counts % 1) != 0).any():
        raise ValueError(f"non-integer counts found in {path.name}")

    return pd.Series(counts.astype("int64").to_numpy(), index=table["gene_id"], name=sample_id_from_filename(path))


def build_matrix(count_dir: Path) -> pd.DataFrame:
    files = sorted(count_dir.glob("GSM*_*.counts.txt.gz"))
    if not files:
        files = sorted(count_dir.glob("GSM*.txt.gz"))
    if not files:
        raise ValueError(f"no GEO count files found in {count_dir}")

    series: list[pd.Series] = []
    reference_genes: pd.Index | None = None
    seen_samples: set[str] = set()

    for path in files:
        counts = read_count_file(path)
        if counts.name in seen_samples:
            raise ValueError(f"duplicate sample accession {counts.name}")
        seen_samples.add(str(counts.name))

        if reference_genes is None:
            reference_genes = counts.index
        elif not counts.index.equals(reference_genes):
            raise ValueError(f"gene rows differ between count files; first mismatch at {path.name}")
        series.append(counts)

    matrix = pd.concat(series, axis=1)
    matrix.index.name = "gene_id"
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count-dir", type=Path, default=Path("data/raw/counts"))
    parser.add_argument("--output", type=Path, default=Path("data/processed/count_matrix.tsv"))
    args = parser.parse_args()

    matrix = build_matrix(args.count_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(args.output, sep="\t")
    print(f"genes={matrix.shape[0]}")
    print(f"samples={matrix.shape[1]}")


if __name__ == "__main__":
    main()
