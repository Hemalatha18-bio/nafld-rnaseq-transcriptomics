# NAFLD RNA-seq Transcriptomics — GSE135251

Reproducible analysis of **public human liver RNA-seq data** from NCBI GEO accession **GSE135251**, focused on transcriptomic changes across NAFLD fibrosis severity.

## Why this project

This repository is the real-public-data component of my bioinformatics portfolio. It demonstrates how I move from a published GEO accession to validated sample metadata, a gene-count matrix, RNA-seq QC, count-aware differential expression, FDR-controlled results, pathway analysis, figures, tests, CI, and a reproducible workflow.

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
   ├── dataset audit
   │   ├── expected 216-sample check
   │   ├── count/metadata sample consistency
   │   └── disease and fibrosis-group summaries
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
       ├── volcano plot
       └── GO Biological Process enrichment
```

### Important statistical boundary

The log2-CPM transform is used **only for exploratory QC/PCA**. Differential expression is performed on the original integer counts with **DESeq2**, so the inferential analysis remains count-aware.

GO enrichment is downstream of the DESeq2 FDR-filtered gene set. Ensembl gene IDs are mapped to Entrez IDs with `org.Hs.eg.db`, and the tested DESeq2 genes are used as the enrichment background rather than the whole genome.

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
│   ├── audit_dataset.py
│   ├── qc_counts.py
│   ├── differential_expression.R
│   └── pathway_enrichment.R
├── tests/
│   ├── test_parse_geo_soft.py
│   ├── test_build_count_matrix.py
│   ├── test_audit_dataset.py
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

The workflow downloads the public GEO inputs automatically, builds the count matrix and metadata, validates the complete dataset, generates QC outputs, runs the default DESeq2 fibrosis contrast, and performs GO Biological Process enrichment.

### Run preprocessing tests only

```bash
python -m pip install -r requirements-dev.txt
pytest -q
```

GitHub Actions runs the Python tests, compiles the scripts, and validates the Snakemake DAG on every pull request and push to `main`.

## Expected result layers

After a complete run, the versioned workflow produces:

- `results/audit/` — cohort/sample consistency and fibrosis-group counts;
- `results/qc/` — sample QC summaries and PCA coordinates;
- `figures/qc/` — library-size and PCA figures;
- `results/deseq2/` — all tested genes, FDR-significant genes, normalized counts, summary, and volcano plot;
- `results/enrichment/` — Ensembl-to-Entrez mapping, GO BP enrichment table, summary, and figure.

Large upstream data and the full generated count matrix stay excluded from Git because they are reproducible from GEO. Small reviewed summaries and selected figures can be committed after a full analysis run.

## Data handling and reproducibility

This repository begins from the **GEO-provided processed gene counts**, not from raw SRA FASTQ re-alignment. That boundary is explicit: upstream Trimmomatic/STAR/HTSeq processing is provenance reported by GEO, not code reproduced in this repository.

## Interpretation policy

No gene, pathway, or performance result is presented as a finding unless it is produced by the versioned workflow from the public data. Software checks, example tests, and biological conclusions are kept clearly separate.

## Author

Hemalatha Ponnam  
M.S. Bioinformatics & Computational Biology  
Bioinformatics & Research Computing
