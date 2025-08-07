# RNA-seq Analysis Pipeline (Nextflow)

This repository contains an RNA-seq analysis pipeline I developed using [Nextflow](https://www.nextflow.io/). It takes paired-end FASTQ files as input and processes them through quality control, trimming, alignment, and statistical analysis.

I created this as part of a project to better understand RNA-seq data analysis and workflow management with Nextflow. This can be used as a base or educational reference for similar projects.

---

## Pipeline Steps

The workflow consists of the following main steps:

1. **FastQC** – run before trimming and after trimming to check read quality.
2. **Trim Galore** – detects low-quality reads and adapter sequences to remove them.
3. **STAR Indexing** – indexes the reference genome.
4. **STAR Alignment** – aligns trimmed reads to the genome.
5. **R Analysis** – an R script (`counts_and_tests.R`) that:
   - Generates count tables
   - Performs differential expression with `edgeR`
   - Uses `clusterProfiler` for enrichment analysis
   - Works with GTF annotations using `rtracklayer` and `org.Hs.eg.db`
   - Generates volcano and heatmap plots

The workflow is written in Nextflow DSL2 and uses Conda for reproducible environments.

---

## Folder Structure
