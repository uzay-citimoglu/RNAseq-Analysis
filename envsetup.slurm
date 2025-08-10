#!/bin/bash -l
#SBATCH --job-name=envsetup
#SBATCH --output=testnf_%j.out
#SBATCH --error=testnf_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH -p boyoz

conda env create -f nfuzay.yml
