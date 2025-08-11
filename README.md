# RNA-seq Analysis Pipeline (Nextflow)

This repository contains a RNA-seq analysis pipeline I developed using [Nextflow](https://www.nextflow.io/). It takes paired-end FASTQ files as input and processes them through quality control, trimming, alignment, and statistical analysis.

I created this as part of a project to better understand RNA-seq data analysis and workflow management with Nextflow. This can be used as a base or educational reference for similar projects. I also used some gtf and fa files from the GENCODE those can be changed the parts that I marked as #EDIT.

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

## Running the Pipeline

### 1. Edit `nfrun.slurm`
Open the file and update the four parameters in the `nextflow run` command.

```bash
#!/bin/bash -l
#SBATCH --job-name=testnf
#SBATCH --output=testnf_%j.out
#SBATCH --error=testnf_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH -p your_partition_name       # EDIT

# --- USER INPUTS (EDIT BELOW) ---
nextflow run mainuzay.nf \
  --reads "/PATH/TO/FASTQ/*_{1,2}.fastq.gz" \                               # EDIT
  --fasta "/PATH/TO/REFERENCE/GENOME.fa" \                                  # EDIT
  --gtf "/PATH/TO/ANNOTATION/GENCODE.gtf" \                                 # EDIT
  --outdir "/PATH/TO/OUTPUT/DIRECTORY" \                                    # EDIT
  -with-report "/PATH/TO/OUTPUT/DIRECTORY/summary_$(date +%F_%H-%M-%S).html" \
  -resume
```
---

### 2. Nextflow Script (`mainuzay.nf`)

The Nextflow workflow defines the full RNA-seq analysis process.  
You must provide the same **four required parameters** (`--reads`, `--fasta`, `--gtf`, `--outdir`) either:

- In the SLURM script (`nfrun.slurm`)
- Or directly on the Nextflow CLI

---

### Parameters in the Script
These are the defaults inside `mainuzay.nf` — you can edit them here **or** override via CLI.

```groovy
// --- USER INPUTS (EDIT THESE PATHS OR OVERRIDE VIA CLI) ---
params.reads  = "/PATH/TO/FASTQ/*_{1,2}.fastq.gz"                            // EDIT
params.fasta  = "/PATH/TO/REFERENCE/GENOME.fa"                               // EDIT
params.gtf    = "/PATH/TO/ANNOTATION/GENCODE.gtf"                            // EDIT
params.outdir = "/PATH/TO/OUTPUT/DIRECTORY"                                  // EDIT
```
---

### 3. R Analysis Script (`counts_and_tests.R`)

This R script performs the downstream RNA-seq analysis after alignment and counting.

#### User Inputs to Edit

| Section / Variable | Description | Example |
|--------------------|-------------|---------|
| `bams <- c(...)` | Paths to sorted BAM files | `/path/to/sample1.sorted.bam` |
| `annot.ext` in `featureCounts` | Path to GTF annotation | `/refs/gencode.v48.annotation.gtf` |
| Output file paths | Where counts and FPKM TSVs are written | `/project/results/counts.tsv` |
| `nthreads` | Number of CPU threads for counting | `20` |
| `group <- factor(...)` | Experimental groups for DE analysis | `c("Control", "Control", "Treatment", "Treatment")` |
| Filtering thresholds | Expression cutoffs and DE cutoffs | Adjust as needed |
| PCA color mapping | Colors assigned to samples in PCA plot | `"Sample1" = "red", "Sample2" = "blue"` |

#### Main Analysis Steps
1. **Counting & FPKM calculation**
2. **Annotation merge with GTF**
3. **Filtering low-expression genes**
4. **PCA plot**
5. **Differential expression analysis**
6. **Volcano plots**
7. **Heatmap**
8. **KEGG & GO enrichment analysis**

####  Output Files
| File | Description |
|------|-------------|
| `counts.tsv` | Raw gene counts |
| `fpkm_values.tsv` | FPKM-normalized values |
| `edgeR_glm_DEG.tsv` | DE results from glmLRT |
| `edgeR_exact_DEG.tsv` | DE results from exactTest |
| `kegg_enrichment.tsv` | KEGG enrichment results |
| `go_enrichment.tsv` | GO enrichment results |
| `kegg_plots.pdf` | KEGG enrichment plots |
| `go_plots.pdf` | GO enrichment plots |

---

## File Structure
- **scripts/**: Contains a Rscript for future usage in statistical analysis.
- **envsetup.slurm**: Contains the bash code for environment setup.
- **mainuzay.nf**: The pipeline.
- **nextflow.config**: Additional setting for the pipeline.
- **nfrun.slurm**: Contains the code and inputs for running the pipeline.
- **nfuzay.yml**: Contains needed packages and tools for pipeline.
