#!/usr/bin/env python3
"""Parse sample metadata from a GEO family SOFT file."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


KNOWN_COLUMNS = [
    "sample_id",
    "title",
    "source_name",
    "nas_score",
    "fibrosis_stage",
    "group_in_paper",
    "disease",
    "stage",
    "description",
    "analysis_group",
]


def normalize_key(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def derive_analysis_group(disease: str | None, fibrosis_stage: str | None) -> str:
    if str(disease).strip().lower() == "control":
        return "Control"
    try:
        fibrosis = int(float(str(fibrosis_stage).strip()))
    except (TypeError, ValueError):
        return "Unknown"
    if fibrosis <= 1:
        return "NAFLD_F0_F1"
    if fibrosis == 2:
        return "NAFLD_F2"
    return "NAFLD_F3_F4"


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_soft(path: Path) -> pd.DataFrame:
    samples: list[dict[str, str]] = []
    current: dict[str, str] | None = None

    with _open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)
                current = {"sample_id": line.split("=", 1)[1].strip()}
                continue
            if current is None:
                continue

            if line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 = "):
                current["source_name"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_description = "):
                current["description"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                payload = line.split("=", 1)[1].strip()
                if ":" in payload:
                    key, value = payload.split(":", 1)
                    current[normalize_key(key)] = value.strip()

    if current is not None:
        samples.append(current)
    if not samples:
        raise ValueError(f"no GEO sample blocks found in {path}")

    frame = pd.DataFrame(samples)
    frame["analysis_group"] = [
        derive_analysis_group(row.get("disease"), row.get("fibrosis_stage"))
        for row in samples
    ]

    for column in KNOWN_COLUMNS:
        if column not in frame.columns:
            frame[column] = pd.NA

    ordered = KNOWN_COLUMNS + sorted(c for c in frame.columns if c not in KNOWN_COLUMNS)
    frame = frame[ordered].drop_duplicates(subset=["sample_id"]).reset_index(drop=True)
    return frame


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--soft", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("metadata/sample_metadata.tsv"))
    args = parser.parse_args()

    metadata = parse_soft(args.soft)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.output, sep="\t", index=False)
    print(f"samples={len(metadata)}")
    print(metadata["analysis_group"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
