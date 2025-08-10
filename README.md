# RNA-seq Analysis Pipeline (Nextflow)

This repository contains an RNA-seq analysis pipeline I developed using [Nextflow](https://www.nextflow.io/). It takes paired-end FASTQ files as input and processes them through quality control, trimming, alignment, and statistical analysis.

I created this as part of a project to better understand RNA-seq data analysis and workflow management with Nextflow. This can be used as a base or educational reference for similar projects.

---

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

## File Guideline

To use this Nextflow pipeline you must first run the `envsetup.slurm`.  
After this, your Conda environment will set up all necessary packages.  
Then run `nfrun.slurm` with your desired file names, and you are good to go!

---

## 📋 Parameters You Must Set

The pipeline has **four required parameters** that must be set before running.

| Parameter    | Description                                               | Where to Set                                         | Example                                                     |
|--------------|-----------------------------------------------------------|------------------------------------------------------|-------------------------------------------------------------|
| `--reads`    | Glob pattern for paired FASTQ files (`R1`/`R2`)            | `nfrun.slurm` or CLI when running Nextflow           | `/data/fastq/*_{1,2}.fastq.gz`                              |
| `--fasta`    | Reference genome FASTA                                     | `nfrun.slurm` or CLI when running Nextflow           | `/refs/GRCh38.primary_assembly.genome.fa`                   |
| `--gtf`      | Gene annotation GTF file                                   | `nfrun.slurm` or CLI when running Nextflow           | `/refs/gencode.v48.primary_assembly.basic.annotation.gtf`   |
| `--outdir`   | Output directory for all pipeline results                  | `nfrun.slurm` or CLI when running Nextflow           | `/project/results`                                          |

---

## 🚀 Running the Pipeline

### 1. Edit `nfrun.slurm`
Open the file and update the four parameters in the `nextflow run` command.

```bash
#!/bin/bash -l
#SBATCH --job-name=testnf
#SBATCH --output=testnf_%j.out
#SBATCH --error=testnf_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH -p boyoz

# --- USER INPUTS (EDIT BELOW) ---
nextflow run mainuzay.nf \
  --reads "/data/fastq/*_{1,2}.fastq.gz" \        # EDIT
  --fasta "/refs/GRCh38.primary_assembly.genome.fa" \    # EDIT
  --gtf "/refs/gencode.v48.primary_assembly.basic.annotation.gtf" \ # EDIT
  --outdir "/project/results" \                   # EDIT
  -with-report "/project/results/summary_$(date +%F_%H-%M-%S).html" \
  -resume
