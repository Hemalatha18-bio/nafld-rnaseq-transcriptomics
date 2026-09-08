# Data provenance

This project uses public processed RNA-seq count files from **NCBI GEO GSE135251**.

GEO accession page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135251

The series contains 216 human liver-biopsy RNA-seq samples (206 NAFLD cases and 10 controls). GEO provides the processed gene-level count files in `GSE135251_RAW.tar` and the original sequencing reads through SRA.

The downloader retrieves:

- `GSE135251_RAW.tar` from the GEO supplementary-file directory;
- `GSE135251_family.soft.gz` from the GEO family SOFT directory.

These downloaded files are excluded from Git because they are public, reproducible upstream inputs. The repository stores code, documentation, small result summaries, and selected figures rather than duplicating the complete public archive.

## Processing provenance reported by GEO

Individual sample records report:

- reads cropped to 75 bp with Trimmomatic;
- alignment to GRCh38 / Ensembl release 78 using STAR;
- gene-level counting with HTSeq-count;
- reverse-stranded counting (`-s reverse`).

This portfolio workflow begins from those GEO-provided processed integer counts. It does not claim to reproduce the upstream FASTQ alignment step unless a future raw-read workflow is explicitly added.
