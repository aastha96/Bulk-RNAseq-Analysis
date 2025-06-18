# Bulk RNA-seq Analysis in R (Yeast: *S. cerevisiae*)

This repository contains a step-by-step pipeline for performing **Bulk RNA-seq analysis** in R using publicly available yeast data. It covers:

- Downloading paired-end FASTQ files
- Downloading reference genome and GTF annotation
- Index building and read alignment using `Rsubread`
- BAM file inspection
- Generation of a count matrix with `featureCounts`

---

## 🧬 Data Sources

- **FASTQ files**: NCBI SRA (e.g., SRR5924196, SRR5924198)
- **Reference genome**: Ensembl
- **Annotation (GTF)**: Ensembl

---

## 🧰 Tools and Packages Used

- R / RStudio
- [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html)
- Rsamtools
- FastQC (optional)
- 
---

## 🧪 Steps in the Workflow

### 1. Set Up & Download Files
- Set working directory
- Download FASTQ files using `download.file()`
- Decompress reference genome and GTF files

### 2. Index the Reference Genome
- Use `buildindex()` from Rsubread

### 3. Align Reads
- Use `align()` to generate BAM files from paired-end FASTQ

### 4. Inspect BAM Files
- Load and check contents using `scanBam()` and `Rsamtools`

### 5. Count Matrix Generation
- Run `featureCounts()` to count reads per gene
- Store output in a `data.frame` for downstream analysis

---
