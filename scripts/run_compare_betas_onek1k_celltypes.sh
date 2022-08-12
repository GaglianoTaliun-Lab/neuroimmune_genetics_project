#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=3:00:00
#SBATCH --job-name=compare_celltypes_onek1k
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8G

module load StdEnv/2020
module load r/4.0.2

R_LIBS=~/.local/R/$EBVERSIONR/.

Rscript compare_betas_onek1k_celltypes.R
