# NAFLD RNA-seq Transcriptomics — GSE135251

Reproducible analysis of **public human liver RNA-seq data** from NCBI GEO accession **GSE135251**, focused on transcriptomic changes across NAFLD fibrosis severity.

## Dataset

NCBI GEO describes GSE135251 as a multicenter human liver transcriptomics study containing **216 snap-frozen liver biopsies: 206 NAFLD cases spanning different fibrosis stages and 10 controls**. Samples were sequenced on the Illumina NextSeq 500 platform.

The public processed files are gene-level count tables. GEO sample records report that reads were cropped to 75 bp with Trimmomatic, aligned to **GRCh38 / Ensembl release 78** with STAR, and quantified with **HTSeq-count**.

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135251
- BioProject: PRJNA558102
- SRA: SRP217231

## Biological Question

**How does the liver transcriptome change as fibrosis severity increases across NAFLD?**

The primary portfolio analysis compares public gene-count profiles across fibrosis groups while keeping the statistical and biological interpretation separate from software-demonstration claims.

## Planned Reproducible Workflow

1. Download the GEO processed count archive and family metadata.
2. Parse sample-level disease, NAS, fibrosis-stage, and paper-group metadata.
3. Build a gene × sample raw-count matrix.
4. Perform count-data QC, library-size diagnostics, filtering, and PCA.
5. Run count-aware differential-expression analysis with DESeq2.
6. Correct for multiple testing with Benjamini-Hochberg FDR.
7. Summarize fibrosis-associated genes and pathway-level patterns.
8. Export machine-readable results and publication-quality figures.
9. Validate preprocessing code with pytest and GitHub Actions.
10. Orchestrate the workflow with Snakemake.

## Reproducibility Boundary

This repository uses only public GEO/SRA data and contains no employer or restricted research data. The initial analysis begins from the GEO-provided processed gene-count files rather than re-aligning the full raw SRA FASTQ collection; that keeps the project reproducible on ordinary hardware while preserving a clear path to a future raw-read workflow.

## Status

Repository initialization complete. The data-ingestion, QC, differential-expression, workflow, tests, and results layers are being added as version-controlled analysis steps.

## Author

Hemalatha Ponnam  
M.S. Bioinformatics & Computational Biology  
Bioinformatics & Research Computing
