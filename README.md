# NAFLD RNA-seq Transcriptomics — GSE135251

Reproducible analysis of **public human liver RNA-seq data** from NCBI GEO accession **GSE135251**, focused on transcriptomic changes across NAFLD fibrosis severity.

## Why this project

This repository is the real-public-data component of my bioinformatics portfolio. It demonstrates how I move from a published GEO accession to validated sample metadata, a gene-count matrix, RNA-seq QC, count-aware differential expression, FDR-controlled results, figures, tests, CI, and a reproducible workflow.

## Dataset

NCBI GEO describes GSE135251 as a multicenter human liver transcriptomics study containing **216 snap-frozen liver biopsies: 206 NAFLD cases spanning different fibrosis stages and 10 controls**. Samples were sequenced on the Illumina NextSeq 500 platform.

The public processed files are gene-level count tables. GEO sample records report that reads were cropped to 75 bp with Trimmomatic, aligned to **GRCh38 / Ensembl release 78** with STAR, and quantified with **HTSeq-count** using reverse-stranded counting.

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135251
- BioProject: PRJNA558102
- SRA: SRP217231

## Biological question

**How does the liver transcriptome change as fibrosis severity increases across NAFLD?**

The default differential-expression contrast is:

**NAFLD fibrosis F3–F4 vs NAFLD fibrosis F0–F1**

Controls and F2 samples are retained in the parsed metadata and can be used for additional comparisons. The repository preserves the original GEO fields and derives a separate `analysis_group` only for transparent downstream grouping.

## Workflow

```text
GEO GSE135251
   │
   ├── GSE135251_RAW.tar ──> 216 HTSeq gene-count files
   │                              │
   │                              └──> validated gene × sample count matrix
   │
   └── GSE135251_family.soft.gz ──> parsed sample metadata
                                      │
                                      ├── disease / NAS / fibrosis stage
                                      └── derived fibrosis analysis group

count matrix + metadata
   │
   ├── exploratory QC
   │   ├── library-size diagnostics
   │   ├── expression filtering
   │   └── PCA on filtered log2-CPM values
   │
   └── DESeq2 on raw integer counts
       ├── F3–F4 vs F0–F1
       ├── Benjamini-Hochberg FDR
       ├── normalized counts
       └── volcano plot + machine-readable results
```

### Important statistical boundary

The log2-CPM transform is used **only for exploratory QC/PCA**. Differential expression is performed on the original integer counts with **DESeq2**, so the inferential analysis remains count-aware.

## Repository structure

```text
nafld-rnaseq-transcriptomics/
├── .github/workflows/ci.yml
├── data/
│   └── README.md
├── metadata/
│   └── README.md
├── scripts/
│   ├── download_geo.py
│   ├── parse_geo_soft.py
│   ├── build_count_matrix.py
│   ├── qc_counts.py
│   └── differential_expression.R
├── tests/
│   ├── test_parse_geo_soft.py
│   ├── test_build_count_matrix.py
│   └── test_qc_counts.py
├── workflow/Snakefile
├── environment.yml
├── requirements-dev.txt
├── CITATION.cff
├── LICENSE
└── README.md
```

## Reproduce the analysis

Create the full Conda environment:

```bash
conda env create -f environment.yml
conda activate nafld-rnaseq
```

Run the complete workflow:

```bash
snakemake --snakefile workflow/Snakefile --cores 4
```

The workflow downloads the public GEO inputs automatically, builds the count matrix and metadata, generates QC outputs, and runs the default DESeq2 fibrosis contrast.

### Run preprocessing tests only

```bash
python -m pip install -r requirements-dev.txt
pytest -q
```

GitHub Actions runs the Python tests, compiles the scripts, and validates the Snakemake DAG on every pull request and push to `main`.

## Data handling and reproducibility

The large GEO downloads and full generated count matrix are intentionally excluded from Git. They are reproducible upstream inputs that can be regenerated from the accession. Small result summaries and selected figures can be versioned after they are produced and reviewed.

This repository begins from the **GEO-provided processed gene counts**, not from raw SRA FASTQ re-alignment. That boundary is explicit: upstream Trimmomatic/STAR/HTSeq processing is provenance reported by GEO, not code reproduced in this repository.

## Interpretation policy

No gene, pathway, or performance result is presented as a finding unless it is produced by the versioned workflow from the public data. Software checks, example tests, and biological conclusions are kept clearly separate.

## Author

Hemalatha Ponnam  
M.S. Bioinformatics & Computational Biology  
Bioinformatics & Research Computing
