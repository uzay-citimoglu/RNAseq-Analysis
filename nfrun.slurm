#!/bin/bash -l
#SBATCH --job-name=testnf
#SBATCH --output=testnf_%j.out
#SBATCH --error=testnf_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH -p boyoz

# --- USER INPUTS (EDIT THE QUOTED PATHS BELOW) --------------------------------
# --reads : Glob for paired FASTQs (R1/R2), e.g. "/path/to/*_{1,2}.fastq.gz"
# --fasta : Reference genome FASTA, e.g. "/refs/GRCh38.primary_assembly.genome.fa"
# --gtf   : Gene annotation GTF, e.g. "/refs/gencode.v48.primary_assembly.basic.annotation.gtf"
# --outdir: Output directory for all results, e.g. "/project/results"
# -with-report: Path to Nextflow HTML report (kept timestamped by default)
# ------------------------------------------------------------------------------

nextflow run mainuzay.nf \
  --reads "/archive/binokayl/Uzay/nftest/data/*_{1,2}.fastq.gz" \   # EDIT
  --fasta "GRCh38.primary_assembly.genome.fa" \                      # EDIT
  --gtf "gencode.v48.primary_assembly.basic.annotation.gtf" \        # EDIT
  --outdir "/archive/binokayl/Uzay/nftest/results_2" \               # EDIT
  -with-report "/archive/binokayl/Uzay/nftest/summary_$(date +%F_%H-%M-%S).html" \
  -resume
