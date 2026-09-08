#!/usr/bin/env python3
"""Download processed GEO count files and family SOFT metadata for a GSE accession."""

from __future__ import annotations

import argparse
import shutil
import tarfile
from pathlib import Path
from urllib.request import urlopen


def geo_series_prefix(accession: str) -> str:
    if not accession.startswith("GSE") or not accession[3:].isdigit():
        raise ValueError("accession must look like GSE135251")
    number = accession[3:]
    if len(number) < 4:
        raise ValueError("GEO series accession is unexpectedly short")
    return f"GSE{number[:-3]}nnn"


def download(url: str, destination: Path, force: bool = False) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not force:
        return destination
    temporary = destination.with_suffix(destination.suffix + ".part")
    with urlopen(url) as response, temporary.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    temporary.replace(destination)
    return destination


def safe_extract_tar(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            target = (destination / member.name).resolve()
            if root not in target.parents and target != root:
                raise ValueError(f"unsafe tar member: {member.name}")
        tar.extractall(destination)


def download_geo(accession: str, raw_dir: Path, metadata_dir: Path, force: bool = False) -> dict[str, Path]:
    prefix = geo_series_prefix(accession)
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}"

    raw_tar = raw_dir / f"{accession}_RAW.tar"
    soft_gz = metadata_dir / f"{accession}_family.soft.gz"
    counts_dir = raw_dir / "counts"

    download(f"{base}/suppl/{accession}_RAW.tar", raw_tar, force=force)
    download(f"{base}/soft/{accession}_family.soft.gz", soft_gz, force=force)

    if force and counts_dir.exists():
        shutil.rmtree(counts_dir)
    if not counts_dir.exists():
        safe_extract_tar(raw_tar, counts_dir)

    return {"raw_tar": raw_tar, "soft": soft_gz, "counts_dir": counts_dir}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--accession", default="GSE135251")
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"))
    parser.add_argument("--metadata-dir", type=Path, default=Path("metadata/raw"))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    paths = download_geo(args.accession, args.raw_dir, args.metadata_dir, force=args.force)
    for key, value in paths.items():
        print(f"{key}={value}")


if __name__ == "__main__":
    main()
