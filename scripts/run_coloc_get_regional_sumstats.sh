#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=2:00:00
#SBATCH --job-name=coloc-get-regional-sumstats
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=4G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

Rscript coloc_get_regional_sumstats.R
