#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=3:00:00
#SBATCH --job-name=coloc-compare-MAFs
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=4G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

Rscript coloc_compare_MAFs.R
