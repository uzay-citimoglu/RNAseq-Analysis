#!/bin/bash -l
#SBATCH --job-name=testnf
#SBATCH --output=testnf_%j.out
#SBATCH --error=testnf_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH -p boyoz

nextflow run mainuzay.nf \
  --reads "/archive/binokayl/Uzay/nftest/data/*_{1,2}.fastq.gz" \
  --fasta "GRCh38.primary_assembly.genome.fa" \
  --gtf "gencode.v48.primary_assembly.basic.annotation.gtf" \
  --outdir "/archive/binokayl/Uzay/nftest/results_2" \
  -with-report "/archive/binokayl/Uzay/nftest/summary_$(date +%F_%H-%M-%S).html" \
  -resume
