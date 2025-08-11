#!/bin/bash -l
#SBATCH --job-name=testnf
#SBATCH --output=testnf_%j.out
#SBATCH --error=testnf_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH -p boyoz

# --- USER INPUTS (EDIT BELOW) ---
nextflow run mainuzay.nf \
  --reads "/PATH/TO/FASTQ/*_{1,2}.fastq.gz" \                               # EDIT
  --fasta "/PATH/TO/REFERENCE/GENOME.fa" \                                  # EDIT
  --gtf "/PATH/TO/ANNOTATION/GENCODE.gtf" \                                 # EDIT
  --outdir "/PATH/TO/OUTPUT/DIRECTORY" \                                    # EDIT
  -with-report "/PATH/TO/OUTPUT/DIRECTORY/summary_$(date +%F_%H-%M-%S).html" \
  -resume
